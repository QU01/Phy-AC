#!/usr/bin/env python3
"""
QUASAR Phy-AC · validation/bench_scm.py — banco de pruebas de la
fidelidad L1 (through-flow por curvatura de líneas de corriente)
=====================================================================
Contesta cinco preguntas sobre `scm_core.py`, todas con número:

  1. COBERTURA   ¿qué fracción del espacio de diseño resuelve, y cuándo
                 no lo hace, POR QUÉ? Un solver que se planta en la mitad
                 de los diseños no es un peldaño de fidelidad usable.
  2. COSTE       segundos por máquina y cómo escala con etapas y con
                 líneas de corriente. L1 dejó de ser gratis cuando dejó
                 de estar muerto; conviene saber cuánto cuesta de verdad.
  3. MALLA       ¿cambia el resultado con el número de líneas? Si cambia,
                 lo que se está midiendo es la malla, no la máquina.
  4. CONVERGENCIA iteraciones y residual del lazo exterior.
  5. RESIDUAL    distribución de L1 − L0. Es la razón de existir de la
                 escalera: si L1 devolviera siempre lo mismo que L0, la
                 capa 2 no tendría residual que aprender y sobraría.

NO es validación: aquí no hay ninguna máquina medida. Esto caracteriza el
solver contra sí mismo y contra L0. La calificación contra medida es
F-02 y sigue abierta.

Uso:
    python validation/bench_scm.py                 # 80 diseños, semilla 71
    python validation/bench_scm.py --n 200 --seed 3
    python validation/bench_scm.py --quick         # 25 diseños, sin malla
"""

from __future__ import annotations

import argparse
import json
import math
import os
import statistics
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import physics_core as pc          # noqa: E402
import scm_core as scm             # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPORT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "BENCH_SCM.md")


def _lhs(n: int, d: int, rng) -> np.ndarray:
    """Latin hypercube en [0,1]^d — el mismo esquema que el optimizador."""
    u = np.empty((n, d))
    for j in range(d):
        u[:, j] = (rng.permutation(n) + rng.random(n)) / n
    return u


def sample_designs(n: int, seed: int) -> list[np.ndarray]:
    """Diseños FACTIBLES en L0. L1 solo corre sobre factibles (el g de L0
    es la puerta), así que el banco tiene que muestrear lo mismo."""
    rng = np.random.default_rng(seed)
    out: list[np.ndarray] = []
    tried = 0
    while len(out) < n and tried < 60 * n:
        u = _lhs(min(4 * n, 400), pc.NDIM, rng)
        for row in u:
            tried += 1
            th = pc.denormalize(row)
            th[0] = float(round(th[0]))          # n_stages entero
            rec = pc.evaluate(th, fidelity=pc.Fidelity.L0, use_cache=False)
            if rec["feasible"]:
                out.append(th)
                if len(out) >= n:
                    break
    return out


def classify(err: str) -> str:
    if "BLOQUEADA" in err:
        return "estacion_bloqueada"
    if "limitador" in err:
        return "limitador_activo"
    if "DERIVA" in err:
        return "deriva_vs_L0"
    if "sin converger" in err:
        return "sin_converger"
    return "otro"


# ---------------------------------------------------------------------------
# 1-2-4-5 · cobertura, coste, convergencia y residual
# ---------------------------------------------------------------------------
def bench_coverage(designs: list[np.ndarray], vortex_n: float) -> dict:
    vn0 = pc.VORTEX_N
    pc.VORTEX_N = vortex_n
    rows = []
    fails: dict[str, int] = {}
    for th in designs:
        rec = pc.evaluate(th, fidelity=pc.Fidelity.L0, use_cache=False)
        if not rec["feasible"]:
            continue
        t0 = time.perf_counter()
        try:
            r = scm.solve(th, rec)
            dt = time.perf_counter() - t0
            rows.append(dict(
                ok=True, dt=dt, n_st=int(rec["n_stages"]),
                it=r["scm"]["iterations"], res=r["scm"]["residual"],
                dPR=r["PR"] / rec["PR"] - 1.0,
                deta=r["eta_poly"] - rec["eta_poly"],
                spread=r["scm"]["work_spread"],
                M_rel=rec["M_rel_tip1"]))
        except scm.SCMDiverged as e:
            dt = time.perf_counter() - t0
            k = classify(str(e))
            fails[k] = fails.get(k, 0) + 1
            rows.append(dict(ok=False, dt=dt, n_st=int(rec["n_stages"]),
                             mode=k))
    pc.VORTEX_N = vn0
    ok = [r for r in rows if r["ok"]]
    return dict(n=len(rows), n_ok=len(ok), fails=fails, rows=rows,
                vortex_n=vortex_n)


# ---------------------------------------------------------------------------
# 3 · independencia de malla
# ---------------------------------------------------------------------------
def bench_grid(designs: list[np.ndarray],
               counts=(5, 7, 9, 11, 13)) -> list[dict]:
    out = []
    for th in designs:
        rec = pc.evaluate(th, fidelity=pc.Fidelity.L0, use_cache=False)
        if not rec["feasible"]:
            continue
        pr, dt = {}, {}
        for c in counts:
            try:
                t0 = time.perf_counter()
                pr[c] = scm.solve(th, rec, n_streamlines=c)["PR"]
                dt[c] = time.perf_counter() - t0
            except scm.SCMDiverged:
                pr[c] = None
                dt[c] = time.perf_counter() - t0
        vals = [v for v in pr.values() if v]
        if len(vals) >= 2:
            ref = pr.get(13) or pr.get(11) or vals[-1]
            out.append(dict(n_st=int(rec["n_stages"]), pr=pr, dt=dt,
                            spread=(max(vals) - min(vals)) / ref,
                            ref=ref))
    return out


# ---------------------------------------------------------------------------
# Informe
# ---------------------------------------------------------------------------
def _stats(xs):
    if not xs:
        return dict(n=0)
    xs = sorted(xs)
    return dict(n=len(xs), med=statistics.median(xs),
                p10=xs[int(0.10 * (len(xs) - 1))],
                p90=xs[int(0.90 * (len(xs) - 1))],
                mn=xs[0], mx=xs[-1],
                mean=sum(xs) / len(xs))


def write_report(cov: list[dict], grid: list[dict], meta: dict) -> str:
    L = ["# Phy-AC — banco de pruebas de la fidelidad L1 (SCM)", "",
         f"Generado por `validation/bench_scm.py` · {meta['n']} diseños "
         f"factibles muestreados por LHS (semilla {meta['seed']}) sobre el "
         "espacio de diseño completo.", "",
         "**No es validación**: no hay ninguna máquina medida aquí. Esto "
         "caracteriza el solver contra sí mismo y contra L0. La "
         "calificación de L1 CONTRA MEDIDA la hace `validate.py` desde la "
         "fase 12.3 (tabla «Mismas máquinas a L1»), y el mapa fuera de "
         "diseño es F-02.", ""]

    # --- 1. cobertura -----------------------------------------------------
    L += ["## 1. Cobertura", "",
          "| Ley de torbellino | diseños | resueltos | estación bloqueada "
          "| limitador activo | deriva vs L0 | sin converger |",
          "|---|---|---|---|---|---|---|"]
    for c in cov:
        f = c["fails"]
        L.append(f"| n = {c['vortex_n']:+.1f} | {c['n']} | "
                 f"**{100 * c['n_ok'] / max(c['n'], 1):.0f}%** "
                 f"({c['n_ok']}) | {f.get('estacion_bloqueada', 0)} | "
                 f"{f.get('limitador_activo', 0)} | "
                 f"{f.get('deriva_vs_L0', 0)} | "
                 f"{f.get('sin_converger', 0)} |")
    L += ["", "Cuando no resuelve lo DICE y `evaluate` degrada a L0 con el "
          "motivo en `source`; nunca devuelve un número que nadie debería "
          "usar.", "",
          "### Cobertura por número de etapas (vórtice libre)", "",
          "| etapas | diseños | resueltos |", "|---|---|---|"]
    by_st_ok: dict[int, list[bool]] = {}
    for r in cov[0]["rows"]:
        by_st_ok.setdefault(r["n_st"], []).append(r["ok"])
    for k in sorted(by_st_ok):
        v = by_st_ok[k]
        L.append(f"| {k} | {len(v)} | {100 * sum(v) / len(v):.0f}% |")
    L += ["", "**Este es el hallazgo central del banco.** La cobertura cae "
          "con el número de etapas, y la razón es estructural, no "
          "numérica: el annulus lo dimensiona L0 con su Cx uniforme y su "
          "densidad media, L1 resuelve un perfil, y el álabe —de ángulo "
          "fijo— convierte esa diferencia en trabajo, que cambia la "
          "densidad, que cambia la siguiente estación. En 1-4 etapas es "
          "ruido; en 7-8 se compone. La cura de fondo es que el annulus "
          "salga del MISMO solver que lo usa; mientras tanto la guarda "
          "`PR_WINDOW` rechaza el punto en vez de devolver el número.", ""]

    # --- 2. coste ---------------------------------------------------------
    base = cov[0]
    ok = [r for r in base["rows"] if r["ok"]]
    s_dt = _stats([r["dt"] for r in ok])
    L += ["## 2. Coste", "",
          f"Sobre los {s_dt['n']} diseños resueltos con vórtice libre:", "",
          f"- mediana **{s_dt['med']:.2f} s** por máquina "
          f"(p10 {s_dt['p10']:.2f}, p90 {s_dt['p90']:.2f}, "
          f"máx {s_dt['mx']:.2f})", ""]
    by_st: dict[int, list[float]] = {}
    for r in ok:
        by_st.setdefault(r["n_st"], []).append(r["dt"])
    L += ["| etapas | n | s/máquina (mediana) | s por etapa |",
          "|---|---|---|---|"]
    for k in sorted(by_st):
        m = statistics.median(by_st[k])
        L.append(f"| {k} | {len(by_st[k])} | {m:.2f} | {m / k:.2f} |")
    L += ["", "El coste crece de forma aproximadamente LINEAL con el número "
          "de etapas: el lazo exterior itera sobre estaciones y cada etapa "
          "añade dos.", ""]

    # --- 3. malla ---------------------------------------------------------
    if grid:
        L += ["## 3. Independencia de malla", "",
              "PR con 5, 7, 9, 11 y 13 líneas de corriente sobre el mismo "
              "diseño. La dispersión es respecto al caso más fino que "
              "resolvió.", "",
              "| etapas | PR 5 | PR 7 | PR 9 | PR 11 | PR 13 | dispersión |",
              "|---|---|---|---|---|---|---|"]
        for g in grid:
            cells = " | ".join(
                f"{g['pr'][c]:.4f}" if g["pr"].get(c) else "—"
                for c in (5, 7, 9, 11, 13))
            L.append(f"| {g['n_st']} | {cells} | "
                     f"**{100 * g['spread']:.2f}%** |")
        sp = _stats([g["spread"] for g in grid])
        L += ["", f"Dispersión mediana **{100 * sp['med']:.2f}%**, "
              f"máxima {100 * sp['mx']:.2f}%. El default son 9 líneas: "
              "por debajo de 7 la derivada radial pierde resolución y por "
              "encima de 11 el coste sube sin que el resultado se mueva.",
              ""]
        # coste vs líneas
        dts: dict[int, list[float]] = {}
        for g in grid:
            for c, v in g["dt"].items():
                dts.setdefault(c, []).append(v)
        L += ["| líneas | s/máquina (mediana) |", "|---|---|"]
        for c in sorted(dts):
            L.append(f"| {c} | {statistics.median(dts[c]):.2f} |")
        L += [""]

    # --- 4. convergencia --------------------------------------------------
    s_it = _stats([r["it"] for r in ok])
    L += ["## 4. Convergencia", "",
          f"Lazo exterior: mediana **{s_it['med']:.0f} iteraciones** "
          f"(p90 {s_it['p90']:.0f}, máx {s_it['mx']:.0f} de "
          f"{scm.SCM_MAX_ITER} permitidas), tolerancia "
          f"{scm.SCM_TOL:g} en movimiento radial relativo.", ""]

    # --- 5. residual ------------------------------------------------------
    L += ["## 5. Residual L1 − L0", "",
          "Es la razón de existir de la escalera de fidelidad: si L1 "
          "devolviera siempre lo mismo que L0, la capa 2 no tendría "
          "residual que aprender.", "",
          "| Ley de torbellino | ΔPR mediana | ΔPR p10–p90 | Δη mediana "
          "| Δη p10–p90 | reparto de trabajo en el span |",
          "|---|---|---|---|---|---|"]
    for c in cov:
        o = [r for r in c["rows"] if r["ok"]]
        if not o:
            continue
        sp = _stats([r["dPR"] for r in o])
        se = _stats([r["deta"] for r in o])
        sw = _stats([r["spread"] for r in o])
        L.append(f"| n = {c['vortex_n']:+.1f} | {100 * sp['med']:+.2f}% | "
                 f"{100 * sp['p10']:+.2f}% … {100 * sp['p90']:+.2f}% | "
                 f"{100 * se['med']:+.2f} pts | "
                 f"{100 * se['p10']:+.2f} … {100 * se['p90']:+.2f} pts | "
                 f"{100 * sw['med']:.1f}% |")
    L += ["", "Con **vórtice libre** el residual pequeño es el resultado "
          "ESPERADO y sirve de verificación: es el caso donde las "
          "hipótesis del meanline (Cx radialmente constante, trabajo "
          "uniforme) son exactas. Con **torbellino controlado** la forma "
          "cerrada del meanline es solo aproximada y ahí aparece el "
          "residual de verdad.", ""]

    # correlación con el Mach relativo de punta
    o = [r for r in base["rows"] if r["ok"]]
    if len(o) > 8:
        x = np.array([r["M_rel"] for r in o])
        y = np.array([r["dPR"] for r in o])
        rr = float(np.corrcoef(x, y)[0, 1])
        L += [f"Correlación de ΔPR con el Mach relativo de punta de la "
              f"etapa 1: **r = {rr:+.2f}** ({len(o)} puntos).", ""]

    L += ["---", "",
          f"Reproducir: `python validation/bench_scm.py --n {meta['n']} "
          f"--seed {meta['seed']}`  ·  "
          f"tiempo total de esta corrida: {meta['wall'] / 60:.1f} min.", ""]
    txt = "\n".join(L)
    with open(REPORT, "w", encoding="utf-8", newline="\n") as f:
        f.write(txt)
    return txt


def main(argv=None) -> int:
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    ap = argparse.ArgumentParser(description="Banco de pruebas de L1 (SCM)")
    ap.add_argument("--n", type=int, default=80)
    ap.add_argument("--seed", type=int, default=71)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-grid", action="store_true")
    a = ap.parse_args(argv)
    n = 25 if a.quick else a.n

    t0 = time.perf_counter()
    print(f"— muestreando {n} diseños factibles (semilla {a.seed})…")
    designs = sample_designs(n, a.seed)
    print(f"  {len(designs)} diseños")

    cov = []
    for vn in (-1.0, -0.5):
        print(f"— cobertura y residual con n = {vn:+.1f}…")
        c = bench_coverage(designs, vn)
        print(f"  {c['n_ok']}/{c['n']} resueltos  fallos={c['fails']}")
        cov.append(c)

    grid = []
    if not (a.quick or a.no_grid):
        print("— independencia de malla (5/7/9/11/13 líneas)…")
        grid = bench_grid(designs[:6])
        print(f"  {len(grid)} diseños")

    meta = dict(n=len(designs), seed=a.seed, wall=time.perf_counter() - t0)
    txt = write_report(cov, grid, meta)
    print()
    print(txt)
    print(f"→ {os.path.relpath(REPORT, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
