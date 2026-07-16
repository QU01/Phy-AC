"""
QUASAR Phy-AC · neural_optimizer.py
===================================
Capas 2–4 del Computational Engineering Model: SURROGATE con incertidumbre,
BÚSQUEDA INVERSA multiobjetivo y ACTIVE LEARNING bayesiano.

Núcleo portado casi literal de Phy-CC (es agnóstico al dominio): deep
ensemble K=5 con residual learning sobre el L0, NSGA-II con dominancia
restringida de Deb, adquisición LCB + diversidad k-means, quality gate y
checkpointing. Solo cambian DesignSpec (espec axial), las restricciones de
espec y el acople a structures_core (peor etapa del tambor).

Dado únicamente un DesignSpec (PR objetivo, flujo másico, límites),
devuelve el frente de Pareto factible y el diseño óptimo VERIFICADO con la
física real — sin que el ingeniero proponga geometría.

Todo NumPy puro: corre en cualquier laptop, sin GPU.
"""

from __future__ import annotations

import json
import time

import numpy as np

from physics_core import (
    Fidelity, evaluate, evaluate_batch, physics_features,
    denormalize, normalize, VAR_NAMES, NDIM, N_FEATURES,
    N_CONSTRAINTS,
)
import structures_core

SEED = 71  # semilla LEAP 71 :)
rng = np.random.default_rng(SEED)

# U2/Mu son alias del record axial (U_tip y Mach de punta) — se conservan
# los nombres de Phy-CC para mantener el núcleo del optimizador intacto.
OUT_KEYS = ["PR", "eta_poly", "U2", "Mu", "power_W"]

#: claves del record L0 que la dominancia restringida necesita además de
#: las predicciones del surrogate (g, térmica y geometría para g_struct)
_REC0_KEYS = ("g", "T02", "stage_table", "AN2_m2_rpm2", "r_tip_in_mm",
              "n_stages", "T0_out")


# ==========================================================================
# 1. MLP base (NumPy, Adam, early stopping) — portado verbatim de Phy-CC
# ==========================================================================
class _Adam:
    """Adam sobre una lista de parámetros NumPy, actualizados in-place."""

    def __init__(self, params, lr=2e-3):
        self.params = list(params)
        self.lr = lr
        self.m = [np.zeros_like(p) for p in self.params]
        self.v = [np.zeros_like(p) for p in self.params]
        self.t = 0

    def step(self, grads):
        self.t += 1
        bias1 = 1.0 - 0.9 ** self.t
        bias2 = 1.0 - 0.999 ** self.t
        for p, g, m, v in zip(self.params, grads, self.m, self.v):
            m[:] = 0.9 * m + 0.1 * g
            v[:] = 0.999 * v + 0.001 * g ** 2
            p -= self.lr * (m / bias1) / (np.sqrt(v / bias2) + 1e-8)


class _MLP:
    def __init__(self, d_in, d_out, hidden=(96, 96, 64), seed=0):
        r = np.random.default_rng(seed)
        sizes = [d_in, *hidden, d_out]
        self.W, self.b = [], []
        for i in range(len(sizes) - 1):
            scale = np.sqrt(2.0 / (sizes[i] + sizes[i + 1]))
            self.W.append(r.normal(0, scale, (sizes[i], sizes[i + 1])))
            self.b.append(np.zeros(sizes[i + 1]))

    @staticmethod
    def _act(z):                       # SiLU
        s = 1.0 / (1.0 + np.exp(-np.clip(z, -30, 30)))
        return z * s

    @staticmethod
    def _act_grad(z):
        s = 1.0 / (1.0 + np.exp(-np.clip(z, -30, 30)))
        return s * (1 + z * (1 - s))

    def forward(self, X, cache=False):
        a, acts, zs = X, [X], []
        for i in range(len(self.W)):
            z = a @ self.W[i] + self.b[i]
            zs.append(z)
            a = self._act(z) if i < len(self.W) - 1 else z
            acts.append(a)
        return (a, acts, zs) if cache else a

    def _backward(self, acts, zs, out, yb, weight_decay):
        delta = 2.0 * (out - yb) / out.shape[0]
        gW = [None] * len(self.W)
        gb = [None] * len(self.b)
        for i in reversed(range(len(self.W))):
            gW[i] = acts[i].T @ delta + weight_decay * self.W[i]
            gb[i] = delta.sum(axis=0)
            if i > 0:
                delta = (delta @ self.W[i].T) * self._act_grad(zs[i - 1])
        return gW, gb

    def train(self, X, Y, Xv, Yv, epochs=600, lr=2e-3, batch=64,
              weight_decay=1e-5, patience=60, seed=0):
        r = np.random.default_rng(seed)
        opt = _Adam(self.W + self.b, lr=lr)
        best_val, best_state, since = np.inf, None, 0
        n = X.shape[0]
        for _ep in range(epochs):
            idx = r.permutation(n)
            for s in range(0, n, batch):
                bi = idx[s:s + batch]
                xb, yb = X[bi], Y[bi]
                out, acts, zs = self.forward(xb, cache=True)
                gW, gb = self._backward(acts, zs, out, yb, weight_decay)
                opt.step(gW + gb)
            val = float(np.mean((self.forward(Xv) - Yv) ** 2))
            if val < best_val - 1e-6:
                best_val, since = val, 0
                best_state = ([w.copy() for w in self.W],
                              [b.copy() for b in self.b])
            else:
                since += 1
                if since >= patience:
                    break
        if best_state:
            self.W, self.b = best_state
        return best_val


# ==========================================================================
# 2. Deep ensemble con residual learning sobre L0
# ==========================================================================
class EnsembleSurrogate:
    """K MLPs entrenados sobre bootstraps distintos. Predice el RESIDUAL
    entre la verdad física (L1/L2 calibrada) y el meanline L0 — el L0 ya
    viene embebido en los features de entrada (physics_features)."""

    def __init__(self, K=5):
        self.K = K
        self.members: list[_MLP] = []
        self.y_mean = self.y_std = None
        self.metrics: dict = {}

    @staticmethod
    def _targets(records: list[dict]) -> np.ndarray:
        return np.array([[r["PR"], r["eta_poly"], r["U2"], r["Mu"],
                          r["power_W"]] for r in records])

    def fit(self, X: np.ndarray, Y: np.ndarray, log=print, rand=None):
        r = rand if rand is not None else rng
        Yt = Y.copy()
        Yt[:, 0] = np.log(np.maximum(Yt[:, 0], 1.001))     # PR -> logPR
        Yt[:, 4] = np.log(np.maximum(Yt[:, 4], 1.0))       # power -> log
        self.y_mean = Yt.mean(axis=0)
        self.y_std = Yt.std(axis=0) + 1e-9
        Yn = (Yt - self.y_mean) / self.y_std

        n = X.shape[0]
        n_val = max(int(0.15 * n), 8)
        perm = r.permutation(n)
        vi, ti = perm[:n_val], perm[n_val:]
        Xv, Yv = X[vi], Yn[vi]

        self.members = []
        for k in range(self.K):
            boot = r.choice(ti, size=len(ti), replace=True)
            m = _MLP(X.shape[1], Yn.shape[1], seed=SEED + 13 * k)
            m.train(X[boot], Yn[boot], Xv, Yv, seed=SEED + 7 * k)
            self.members.append(m)

        # ---- Quality gate: R² en holdout + cobertura de incertidumbre ----
        mu, sd = self._predict_norm(Xv)
        r2 = {}
        for i, key in enumerate(OUT_KEYS):
            ss_res = np.sum((Yv[:, i] - mu[:, i]) ** 2)
            ss_tot = np.sum((Yv[:, i] - Yv[:, i].mean()) ** 2) + 1e-12
            r2[key] = float(1 - ss_res / ss_tot)
        inside = np.abs(Yv - mu) <= 2.0 * np.maximum(sd, 1e-6)
        coverage = float(inside.mean())
        self.metrics = dict(r2=r2, coverage_2sigma=coverage, n_train=n)
        log(f"      Surrogate R²: " +
            "  ".join(f"{k}={v:.3f}" for k, v in r2.items()) +
            f"  | coverage ±2σ = {coverage:.2f}")
        return self.metrics

    def gate_ok(self) -> bool:
        r2 = self.metrics.get("r2", {})
        cov = self.metrics.get("coverage_2sigma", 0.0)
        return (r2.get("PR", 0) >= 0.90 and r2.get("eta_poly", 0) >= 0.85
                and 0.80 <= cov <= 1.0)

    def _predict_norm(self, X):
        preds = np.stack([m.forward(X) for m in self.members])  # K x n x d
        return preds.mean(axis=0), preds.std(axis=0)

    def predict(self, X):
        """Devuelve (mu, sigma) en unidades físicas."""
        mu_n, sd_n = self._predict_norm(np.atleast_2d(X))
        mu = mu_n * self.y_std + self.y_mean
        sd = sd_n * self.y_std
        mu[:, 0], sd[:, 0] = np.exp(mu[:, 0]), np.exp(mu[:, 0]) * sd[:, 0]
        mu[:, 4], sd[:, 4] = np.exp(mu[:, 4]), np.exp(mu[:, 4]) * sd[:, 4]
        return mu, sd


# ==========================================================================
# 3. Especificación y objetivos
# ==========================================================================
class DesignSpec:
    """Lo único que el ingeniero escribe.

    `material`: aleación del rotor. Activa las 4 restricciones
    estructurales duras de structures_core (fluencia a sobrevelocidad,
    margen de estallido, límite AN², raíz de álabe) dentro de la
    dominancia restringida. `material=None` las desactiva (solo
    aerodinámica — para estudios/debug)."""

    def __init__(self, PR_target=4.0, massflow=25.0, RPM_max=20_000,
                 U_tip_max=480.0, r_tip_max_mm=400.0, n_stages_max=8,
                 power_max_W=None, T0_in=288.15, P0_in=101_325.0,
                 material="Ti-6Al-4V"):
        self.PR_target = PR_target
        self.massflow = massflow
        self.RPM_max = RPM_max
        self.U_tip_max = U_tip_max
        self.r_tip_max_mm = r_tip_max_mm
        self.n_stages_max = n_stages_max
        self.power_max_W = power_max_W
        self.T0_in = T0_in
        self.P0_in = P0_in
        if material is not None and material not in structures_core.MATERIALS:
            raise ValueError(
                f"material '{material}' no existe — disponibles: "
                f"{', '.join(structures_core.list_materials())} (o None)")
        self.material = material

    def fix_operating_point(self, theta: np.ndarray) -> np.ndarray:
        theta = np.asarray(theta, dtype=float).copy()
        theta[10], theta[11], theta[12] = (self.T0_in, self.P0_in,
                                           self.massflow)
        return theta

    # ---- objetivos (a MINIMIZAR) y restricciones (g<=0) sobre un record --
    def objectives(self, perf: dict) -> np.ndarray:
        f1 = -perf["eta_poly"]                                    # max eta
        f2 = abs(perf["PR"] - self.PR_target) / self.PR_target    # match PR
        return np.array([f1, f2])

    def constraints(self, perf: dict, theta: np.ndarray) -> np.ndarray:
        """[g_phys aero (8), U_tip, RPM, r_tip, n_stages, (potencia),
        g_struct (4)].

        g_struct es EXACTO (discos por etapa sobre el stage_table del
        meanline, ~1 ms), no una predicción del surrogate. Requiere que
        `perf` traiga stage_table/AN2 (todo record físico los trae;
        NSGA2._eval_pop los copia del meanline L0 del individuo)."""
        g_phys = np.asarray(perf.get("g", [0.0] * N_CONSTRAINTS), dtype=float)
        g_spec = [perf["U2"] / self.U_tip_max - 1.0,
                  theta[1] / self.RPM_max - 1.0,
                  perf.get("r_tip_in_mm", 0.0) / self.r_tip_max_mm - 1.0,
                  theta[0] / self.n_stages_max - 1.0]
        if self.power_max_W:
            g_spec.append(perf["power_W"] / self.power_max_W - 1.0)
        parts = [g_phys, np.array(g_spec)]
        if self.material is not None:
            try:
                s = structures_core.evaluate_structural(
                    theta, perf, material=self.material)
                parts.append(np.asarray(s["g_struct"], dtype=float))
            except Exception:
                # theta degenerado: violación grande y continua, no crash
                parts.append(np.full(structures_core.N_STRUCT_CONSTRAINTS,
                                     10.0))
        return np.concatenate(parts)

    @staticmethod
    def violation(g: np.ndarray) -> float:
        return float(np.sum(np.maximum(g, 0.0)))


# ==========================================================================
# 4. NSGA-II con dominancia restringida de Deb — portado verbatim
# ==========================================================================
def _dominates(f1, v1, f2, v2) -> bool:
    """Dominancia restringida: factible >> infactible; entre infactibles
    gana la menor violación; entre factibles, Pareto clásico."""
    feas1, feas2 = v1 <= 1e-9, v2 <= 1e-9
    if feas1 and not feas2:
        return True
    if feas2 and not feas1:
        return False
    if not feas1 and not feas2:
        return v1 < v2
    return bool(np.all(f1 <= f2) and np.any(f1 < f2))


def _fast_nondominated_sort(F, V):
    n = len(F)
    S = [[] for _ in range(n)]
    n_dom = np.zeros(n, dtype=int)
    fronts = [[]]
    for p in range(n):
        for q in range(n):
            if p == q:
                continue
            if _dominates(F[p], V[p], F[q], V[q]):
                S[p].append(q)
            elif _dominates(F[q], V[q], F[p], V[p]):
                n_dom[p] += 1
        if n_dom[p] == 0:
            fronts[0].append(p)
    i = 0
    while fronts[i]:
        nxt = []
        for p in fronts[i]:
            for q in S[p]:
                n_dom[q] -= 1
                if n_dom[q] == 0:
                    nxt.append(q)
        i += 1
        fronts.append(nxt)
    return fronts[:-1]


def _crowding(F, idxs):
    m = F.shape[1]
    d = np.zeros(len(idxs))
    for j in range(m):
        order = np.argsort(F[idxs, j])
        d[order[0]] = d[order[-1]] = np.inf
        span = F[idxs[order[-1]], j] - F[idxs[order[0]], j] + 1e-12
        for k in range(1, len(idxs) - 1):
            d[order[k]] += (F[idxs[order[k + 1]], j]
                            - F[idxs[order[k - 1]], j]) / span
    return d


class NSGA2:
    """NSGA-II sobre el surrogate (objetivos LCB-conservadores)."""

    def __init__(self, spec: DesignSpec, surrogate: EnsembleSurrogate,
                 pop=96, generations=80, kappa=1.0, seed=SEED):
        self.spec, self.sur = spec, surrogate
        self.pop, self.gens, self.kappa = pop, generations, kappa
        self.rng = np.random.default_rng(seed)

    # -- evaluación barata de una población (solo surrogate + L0 features) -
    def _eval_pop(self, U):
        thetas = np.array([self.spec.fix_operating_point(denormalize(u))
                           for u in U])
        recs0 = [evaluate(t, fidelity=Fidelity.L0, use_cache=False,
                          calibrate=False)
                 for t in thetas]
        X = np.array([physics_features(t, rec)
                      for t, rec in zip(thetas, recs0)])
        mu, sd = self.sur.predict(X)
        F = np.empty((len(U), 2))
        V = np.empty(len(U))
        Sig = np.empty(len(U))
        for i, (t, m, s, rec0) in enumerate(zip(thetas, mu, sd, recs0)):
            # LCB conservador: eta penalizada por sigma, |dPR| ensanchado
            perf = dict(PR=m[0], eta_poly=m[1] - self.kappa * s[1],
                        U2=m[2], Mu=m[3] + self.kappa * s[3],
                        power_W=m[4])
            for k in _REC0_KEYS:
                perf[k] = rec0[k]
            F[i] = self.spec.objectives(perf)
            F[i, 1] += self.kappa * s[0] / max(self.spec.PR_target, 1e-9)
            V[i] = self.spec.violation(self.spec.constraints(perf, t))
            Sig[i] = float(np.mean(s / (np.abs(m) + 1e-9)))
        return thetas, F, V, Sig

    def run(self):
        d = NDIM
        P = self.rng.random((self.pop, d))
        thetas, F, V, Sig = self._eval_pop(P)
        for _g in range(self.gens):
            children = np.empty_like(P)
            fronts = _fast_nondominated_sort(F, V)
            rank = np.empty(self.pop, dtype=int)
            for fi, fr in enumerate(fronts):
                rank[fr] = fi
            for i in range(self.pop):
                a, b = self.rng.integers(self.pop, size=2)
                p1 = a if (rank[a], V[a]) < (rank[b], V[b]) else b
                a, b = self.rng.integers(self.pop, size=2)
                p2 = a if (rank[a], V[a]) < (rank[b], V[b]) else b
                w = self.rng.random(d)
                child = w * P[p1] + (1 - w) * P[p2]
                mut = self.rng.random(d) < (1.0 / d) * 2.0
                child[mut] += self.rng.normal(0, 0.08, mut.sum())
                children[i] = np.clip(child, 0, 1)
            cthetas, cF, cV, cSig = self._eval_pop(children)
            allP = np.vstack([P, children])
            allF = np.vstack([F, cF])
            allV = np.concatenate([V, cV])
            allS = np.concatenate([Sig, cSig])
            allT = np.vstack([thetas, cthetas])
            fronts = _fast_nondominated_sort(allF, allV)
            keep = []
            for fr in fronts:
                if len(keep) + len(fr) <= self.pop:
                    keep.extend(fr)
                else:
                    cd = _crowding(allF, np.array(fr))
                    order = np.argsort(-cd)
                    keep.extend(np.array(fr)[order[: self.pop - len(keep)]])
                    break
            keep = np.array(keep)
            P, F, V, Sig, thetas = (allP[keep], allF[keep], allV[keep],
                                    allS[keep], allT[keep])
        fronts = _fast_nondominated_sort(F, V)
        pareto = np.array(fronts[0])
        return dict(U=P, thetas=thetas, F=F, V=V, Sig=Sig, pareto_idx=pareto)


# ==========================================================================
# 5. Adquisición: LCB optimista + diversidad k-means — portado verbatim
# ==========================================================================
def select_acquisition_batch(result: dict, spec: DesignSpec,
                             batch_size=14, explore_frac=0.3,
                             seed=SEED) -> np.ndarray:
    """Elige el batch a evaluar con la física real: mezcla de (a) frente de
    Pareto predicho (explotación) y (b) puntos de alta incertidumbre
    relativa (exploración), des-duplicados por k-means."""
    r = np.random.default_rng(seed)
    U, F, V, Sig = result["U"], result["F"], result["V"], result["Sig"]
    pareto = result["pareto_idx"]

    n_explore = max(1, int(batch_size * explore_frac))

    cand = list(pareto[np.argsort(V[pareto])])
    near = np.where(V <= np.quantile(V, 0.5))[0]
    cand_exp = list(near[np.argsort(-Sig[near])][: 4 * n_explore])

    pool_idx = list(dict.fromkeys(cand + cand_exp))   # únicos, ordenados
    pool = U[pool_idx]
    if len(pool) <= batch_size:
        return pool

    k = batch_size
    centers = pool[r.choice(len(pool), size=k, replace=False)]
    for _ in range(12):
        dist = np.linalg.norm(pool[:, None, :] - centers[None], axis=2)
        lab = dist.argmin(axis=1)
        for j in range(k):
            pts = pool[lab == j]
            if len(pts):
                centers[j] = pts.mean(axis=0)
    chosen_idx: list[int] = []
    for j in range(k):
        members = np.where(lab == j)[0]
        if len(members):
            chosen_idx.append(int(members[0]))
    for i in range(len(pool)):
        if len(chosen_idx) >= batch_size:
            break
        if i not in chosen_idx:
            chosen_idx.append(i)
    return pool[np.array(chosen_idx[:batch_size])]


# ==========================================================================
# 6. Orquestador: diseño inverso autónomo con active learning
# ==========================================================================
class AutonomousAxialDesigner:
    def __init__(self, spec: DesignSpec, fidelity=Fidelity.L1,
                 n_workers=1, log=print, seed=SEED):
        self.spec = spec
        self.fidelity = fidelity
        self.n_workers = n_workers
        self.log = log
        self.rng = np.random.default_rng(seed)
        self.sur = EnsembleSurrogate(K=5)
        self.X: list[np.ndarray] = []     # features
        self.recs: list[dict] = []        # records físicos
        self.history: list[dict] = []

    # ---- evaluación física + registro ------------------------------------
    def _eval_and_store(self, thetas: np.ndarray) -> list[dict]:
        thetas = np.array([self.spec.fix_operating_point(t) for t in thetas])
        recs = evaluate_batch(thetas, fidelity=self.fidelity,
                              n_workers=self.n_workers)
        for t, r in zip(thetas, recs):
            self.X.append(physics_features(t))
            r["_theta"] = t.tolist()
            self.recs.append(r)
        return recs

    @staticmethod
    def _lhs(n, d, r):
        cut = np.linspace(0, 1, n + 1)
        u = r.random((n, d))
        s = cut[:n, None] + u / n
        for j in range(d):
            s[:, j] = r.permutation(s[:, j])
        return s

    def _train(self):
        X = np.array(self.X)
        Y = EnsembleSurrogate._targets(self.recs)
        self.log(f"[surrogate] training ensemble on {len(X)} samples...")
        self.sur.fit(X, Y, log=self.log, rand=self.rng)
        if not self.sur.gate_ok():
            self.log("      ⚠ quality gate FAILED → extra round of "
                     "space-filling sampling before optimizing.")
        return self.sur.gate_ok()

    def _scalar_fitness(self, rec: dict, theta: np.ndarray) -> float:
        """Fitness escalar SOLO para reportar/elegir incumbente."""
        g = self.spec.constraints(rec, theta)
        v = self.spec.violation(g)
        if v > 1e-9:
            return -10.0 * v
        pr_err = abs(rec["PR"] - self.spec.PR_target) / self.spec.PR_target
        return rec["eta_poly"] - 2.0 * pr_err

    def _update_best(self, recs: list[dict], best):
        for r in recs:
            t = np.array(r["_theta"])
            fit = self._scalar_fitness(r, t)
            if best is None or fit > best[2]:
                best = (t, r, fit)
        return best

    def run(self, n_init=320, rounds=5, batch_size=14,
            checkpoint_path: str | None = None) -> dict:
        if n_init < 1:
            raise ValueError("n_init debe ser >= 1")
        t0 = time.time()
        self.log(f"[1/4] Initial LHS sampling: {n_init} physical evaluations "
                 f"({Fidelity(self.fidelity).name})...")
        U0 = self._lhs(n_init, NDIM, self.rng)
        recs0 = self._eval_and_store(np.array([denormalize(u) for u in U0]))

        gate = self._train()
        best = self._update_best(recs0, None)

        for rd in range(rounds):
            self.log(f"\n===== Round {rd + 1}/{rounds} =====")
            if not gate:
                extra = self._lhs(batch_size * 2, NDIM, self.rng)
                recs = self._eval_and_store(
                    np.array([denormalize(u) for u in extra]))
                best = self._update_best(recs, best)
                gate = self._train()
                continue

            self.log("[2/4] Constrained NSGA-II over surrogate...")
            ga = NSGA2(self.spec, self.sur, pop=96, generations=60,
                       seed=SEED + rd)
            result = ga.run()
            n_feas = int(np.sum(result["V"][result["pareto_idx"]] <= 1e-9))
            self.log(f"      Predicted front: {len(result['pareto_idx'])} "
                     f"points ({n_feas} feasible).")

            self.log("[3/4] Acquisition (LCB + diversity): "
                     f"{batch_size} physical evaluations...")
            batchU = select_acquisition_batch(result, self.spec,
                                              batch_size=batch_size,
                                              seed=SEED + 100 * rd)
            recs = self._eval_and_store(
                np.array([denormalize(u) for u in batchU]))

            scored = [(r, np.array(r["_theta"]),
                       self._scalar_fitness(r, np.array(r["_theta"])))
                      for r in recs]
            scored.sort(key=lambda z: z[2], reverse=True)
            cand_rec, cand_theta, cand_fit = scored[0]
            best = self._update_best(recs, best)
            self.history.append(dict(
                round=rd + 1, best_fitness=cand_fit,
                PR=cand_rec["PR"], eta=cand_rec["eta_poly"],
                Mu=cand_rec["Mu"], n_evals=len(self.recs),
                surrogate=self.sur.metrics))
            self.log(f"      Best verified in round: "
                     f"PR={cand_rec['PR']:.3f}  eta={cand_rec['eta_poly']:.3f}"
                     f"  Mu={cand_rec['Mu']:.3f}  fit={cand_fit:.3f}")

            gate = self._train()
            if checkpoint_path:
                self.save(checkpoint_path)

        if checkpoint_path:
            self.save(checkpoint_path)   # siempre: aun sin gate hay dataset
        self.log(f"\n[4/4] Completed in {time.time() - t0:.1f}s · "
                 f"{len(self.recs)} total physical evaluations.")
        theta, rec, fit = best
        return dict(theta=theta, record=rec, fitness=fit,
                    history=self.history,
                    pareto_front=self._final_pareto())

    def _final_pareto(self) -> list[dict]:
        """Frente de Pareto VERIFICADO (solo evaluaciones físicas reales)."""
        F, V, keep = [], [], []
        for r in self.recs:
            t = np.array(r["_theta"])
            F.append(self.spec.objectives(r))
            V.append(self.spec.violation(self.spec.constraints(r, t)))
            keep.append(r)
        F, V = np.array(F), np.array(V)
        fronts = _fast_nondominated_sort(F, V)
        out = []
        for i in fronts[0]:
            r = keep[i]
            out.append(dict(theta=r["_theta"], PR=r["PR"],
                            eta_poly=r["eta_poly"], Mu=r["Mu"],
                            n_stages=r.get("n_stages"),
                            power_W=r["power_W"], feasible=bool(V[i] <= 1e-9),
                            source=r.get("source", "")))
        return out

    def save(self, path: str):
        payload = dict(
            spec=vars(self.spec), seed=SEED,
            n_evals=len(self.recs), history=self.history,
            dataset=[dict(theta=r["_theta"],
                          **{k: r[k] for k in OUT_KEYS},
                          g=r.get("g"), source=r.get("source"))
                     for r in self.recs])
        with open(path, "w") as f:
            json.dump(payload, f, indent=1)


# alias con el nombre del patrón Phy-CC (scripts genéricos)
AutonomousCompressorDesigner = AutonomousAxialDesigner


# ==========================================================================
# 7. API de producto
# ==========================================================================
def design(spec: DesignSpec, rounds=5, n_init=320, batch_size=14,
           fidelity=Fidelity.L1, n_workers=1,
           checkpoint_path="phyac_run.json", log=print, seed=SEED) -> dict:
    """Punto de entrada Phy-AC: espec → diseño verificado + Pareto."""
    d = AutonomousAxialDesigner(spec, fidelity=fidelity,
                                n_workers=n_workers, log=log, seed=seed)
    return d.run(n_init=n_init, rounds=rounds, batch_size=batch_size,
                 checkpoint_path=checkpoint_path)


if __name__ == "__main__":
    import sys
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    # Espec de demo: compresor industrial de ~PR 4
    spec = DesignSpec(PR_target=4.0, massflow=25.0, RPM_max=18_000,
                      U_tip_max=460.0, r_tip_max_mm=400.0)
    out = design(spec, rounds=2, n_init=120, batch_size=10,
                 fidelity=Fidelity.L0)

    print("\n" + "=" * 64)
    print("VERIFIED AUTONOMOUS OPTIMAL DESIGN — Quasar Phy-AC")
    print("=" * 64)
    for name, val in zip(VAR_NAMES, out["theta"]):
        print(f"  {name:10s} = {val:12.3f}")
    print("-" * 64)
    rec = out["record"]
    for k in ["PR", "eta_poly", "eta_isen", "U_tip", "Mu", "M_rel_tip1",
              "min_SM", "n_stages", "length_mm", "power_W", "source"]:
        if k in rec:
            v = rec[k]
            print(f"  {k:14s} = {v:.3f}" if isinstance(v, float)
                  else f"  {k:14s} = {v}")
    print("-" * 64)
    feas = [p for p in out["pareto_front"] if p["feasible"]]
    print(f"Verified Pareto front: {len(out['pareto_front'])} points "
          f"({len(feas)} feasible). Checkpoint: phyac_run.json")
