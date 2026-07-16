#!/usr/bin/env python3
"""
QUASAR Phy-AC · test_phyac.py — suite de VERIFICACIÓN (VV&UQ).

Distinción del proyecto: este script VERIFICA (¿resolvemos bien las
ecuaciones?); validation/validate.py VALIDA (¿las ecuaciones correctas?
— contra máquinas NASA medidas). CI corre ambos.

    python test_phyac.py            # todos los checks, exit 1 si falla
"""

from __future__ import annotations

import json
import math
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import physics_core as pc
import structures_core as sc
import blade_profiles as bp
import geometry_generator as gg
import data_pipeline as dp
import neural_optimizer as no

THETA_REF = np.array([4.0, 12_500.0, 0.62, 0.55, 0.30, -0.10,
                      0.60, 1.20, 1.10, 2.20, 288.15, 101_325.0, 25.0])

_n_pass = _n_fail = 0


def check(name: str, ok, detail: str = ""):
    global _n_pass, _n_fail
    ok = bool(ok)
    _n_pass += ok
    _n_fail += not ok
    mark = "OK " if ok else "FAIL"
    print(f"  [{mark}] {name}" + (f"  ({detail})" if detail else ""))


# ==========================================================================
print("— T1 · espacio de diseño")
u = np.random.default_rng(0).random(pc.NDIM)
check("normalize∘denormalize = identidad",
      np.allclose(pc.normalize(pc.denormalize(u)), u, atol=1e-12))
check("bounds coherentes (lo < hi)", np.all(pc.BOUNDS_LO < pc.BOUNDS_HI))
check("NDIM = 13", pc.NDIM == 13)

# ==========================================================================
print("— T2 · identidades de triángulos y consistencia termodinámica")
rec = pc._meanline(THETA_REF)
ok_tri = True
for s in rec["stage_table"]:
    lhs = math.tan(math.radians(s["beta1_deg"])) + \
        math.tan(math.radians(s["alpha1_deg"]))
    ok_tri &= abs(lhs - 1.0 / s["phi"]) < 1e-6
check("tanβ1 + tanα1 = 1/φ en todas las etapas", ok_tri)

dh_sum = sum(s["dh0"] for s in rec["stage_table"])
dT = rec["T0_out"] - THETA_REF[10]
check("Euler: ΣΔh0 = cp·ΔT0 (±0.1%)",
      abs(dh_sum - pc.CP * dT) / (pc.CP * dT) < 1e-3,
      f"{dh_sum:.0f} vs {pc.CP * dT:.0f} J/kg")

pr_prod = float(np.prod([s["PR"] for s in rec["stage_table"]]))
check("PR etapas ≈ PR máquina (antes del OGV, ±2%)",
      abs(pr_prod - rec["PR"]) / rec["PR"] < 0.02,
      f"{pr_prod:.3f} vs {rec['PR']:.3f}")

# límite sin pérdidas → η→1
_k = (pc.K_PROFILE, pc.K_ENDWALL, pc.K_SHOCK, pc.K_TIP_CLEARANCE)
pc.K_PROFILE = pc.K_ENDWALL = pc.K_SHOCK = pc.K_TIP_CLEARANCE = 0.0
rec_is = pc._meanline(THETA_REF)
pc.K_PROFILE, pc.K_ENDWALL, pc.K_SHOCK, pc.K_TIP_CLEARANCE = _k
check("límite isentrópico: η ≥ 0.995 sin pérdidas (pre-OGV η_tt)",
      all(s["eta_tt"] >= 0.995 for s in rec_is["stage_table"]),
      f"η_poly máquina {rec_is['eta_poly']:.4f}")

# ==========================================================================
print("— T3 · sanidad física")
th_hi = THETA_REF.copy()
th_hi[4] = 0.65                     # ψ muy por encima del óptimo de Smith
check("ψ≫óptimo ⇒ η↓",
      pc._meanline(th_hi)["eta_poly"] < rec["eta_poly"])
th_rpm = THETA_REF.copy()
th_rpm[1] *= 1.15
check("RPM↑ (mismo ṁ) ⇒ PR↑", pc._meanline(th_rpm)["PR"] > rec["PR"])
check("reacción de hub < reacción media (free vortex)",
      all(s["Rx_hub"] < s["Rx"] for s in rec["stage_table"]))
check("altura de álabe decrece hacia atrás",
      all(rec["stage_table"][i]["h_blade_mm"]
          > rec["stage_table"][i + 1]["h_blade_mm"]
          for i in range(len(rec["stage_table"]) - 1)))
check("Mach axial subsónico en todas las estaciones",
      all(s["Mx"] < 1.0 for s in rec["stage_table"]))

# ==========================================================================
print("— T4 · robustez del vector g (dominancia de Deb)")
rng = np.random.default_rng(42)
n_exc = 0
n_feas = 0
gs = []
for _ in range(500):
    th = pc.denormalize(rng.random(pc.NDIM))
    try:
        r = pc._meanline(th)
        gs.append(r["g"])
        n_feas += r["feasible"]
    except Exception:
        n_exc += 1
gs = np.array(gs)
check("LHS(500): 0 excepciones", n_exc == 0, f"{n_exc}")
check("LHS(500): g siempre finito", np.all(np.isfinite(gs)))
check("LHS(500): ≥ 10% factible", n_feas >= 50, f"{n_feas / 5:.0f}%")

# continuidad de g cerca de choke (barrido de ṁ)
th = THETA_REF.copy()
g4_prev, max_jump = None, 0.0
for mm in np.linspace(20.0, 60.0, 81):
    th[12] = mm
    g4 = pc._meanline(th)["g"][4]
    if g4_prev is not None:
        max_jump = max(max_jump, abs(g4 - g4_prev))
    g4_prev = g4
check("g[4] continuo a través del choke (salto máx < 0.6)",
      max_jump < 0.6, f"salto máx {max_jump:.3f}")

th_a = THETA_REF.copy(); th_a[0] = 3.4
th_b = THETA_REF.copy(); th_b[0] = 3.6
check("n_stages se redondea (3.4→3, 3.6→4)",
      pc._meanline(th_a)["n_stages"] == 3
      and pc._meanline(th_b)["n_stages"] == 4)

# ==========================================================================
print("— T5 · off-design y mapa")
od = pc.offdesign(THETA_REF, rpm=float(THETA_REF[1]),
                  mdot=float(THETA_REF[12]))
check("off-design en el punto de diseño reproduce PR (±1%)",
      abs(od["PR"] - rec["PR"]) / rec["PR"] < 0.01,
      f"{od['PR']:.3f} vs {rec['PR']:.3f}")
check("off-design en diseño: incidencia ≈ 0",
      all(abs(s["incidence_deg"]) < 0.5 for s in od["stage_table"]))
mp = pc.compressor_map(THETA_REF, speed_fracs=(0.9, 1.0), n_points=11)
sl = mp["speedlines"][-1]
check("mapa: speedline de diseño tiene rango estable",
      sl["mdot_surge"] is not None and sl["mdot_choke"] is not None)
if sl["mdot_surge"] is not None:
    check("mapa: ṁ_choke > ṁ_surge", sl["mdot_choke"] >= sl["mdot_surge"])

# ==========================================================================
print("— T6 · perfiles y ángulos metálicos")
for name, fn in (("NACA65", bp.naca65_profile), ("DCA", bp.dca_profile)):
    p = fn(25.0, 0.08)
    x, y = p[:, 0], p[:, 1]
    area = 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)
    check(f"{name}: {bp.N_POINTS} pts, CCW, área positiva",
          len(p) == bp.N_POINTS and area > 0.0, f"área {area:.4f}")
devs = [bp.carter_deviation(c, 30.0, 1.2) for c in (10, 20, 30, 40)]
check("desviación de Carter monótona en comba",
      all(devs[i] < devs[i + 1] for i in range(3)))
ma = bp.metal_angles(55.0, 35.0, 1.2, 0.08)
check("punto fijo de ángulos metálicos: comba > giro de flujo",
      ma["camber_deg"] > (55.0 - 35.0),
      f"θ={ma['camber_deg']:.1f}° vs Δβ=20°")
check("desviación positiva y < 15°", 0 < ma["deviation_deg"] < 15)

# ==========================================================================
print("— T7 · contrato de geometría (capa 5a → 5c)")
contract = gg.build_contract(THETA_REF, rec, None, None)
check("schema phyac-axial-1", contract["schema"] == "phyac-axial-1")
ok_pts = ok_cam = ok_fv = True
for st in contract["stages"]:
    for kind in ("rotor", "stator"):
        secs = st[kind]["sections"]
        ok_pts &= all(len(s["points"]) == bp.N_POINTS for s in secs)
        ok_cam &= len({len(s["camber_points"]) for s in secs}) == 1
        # free vortex: r·Cu const ⇒ ángulo de entrada crece con r (rotor)
    r_secs = st["rotor"]["sections"]
    ok_fv &= all(r_secs[i]["metal_in_deg"] < r_secs[i + 1]["metal_in_deg"]
                 for i in range(len(r_secs) - 1))
check("todas las secciones con conteo de puntos idéntico", ok_pts)
check("líneas de comba con conteo idéntico (loft 5c)", ok_cam)
check("free vortex: β1 metálico crece hub→punta (rotor)", ok_fv)
check("estátor con stagger negativo (convención firmada)",
      all(s["stagger_deg"] < 0
          for s in contract["stages"][0]["stator"]["sections"]))
# sin solape axial entre filas consecutivas (la extensión del contrato ya
# incluye la inflación del shell ±½t — regresión del defecto visto en STL)
row_seq = []
for st in contract["stages"]:
    row_seq += [st["rotor"], st["stator"]]
check("filas sin solape axial (incluida inflación del shell)",
      all(row_seq[i + 1]["z_le_mm"] >= row_seq[i]["z_te_mm"] - 1e-6
          for i in range(len(row_seq) - 1)))
# la extensión axial de cada sección cabe en su hueco
ok_ext = True
for rw in row_seq:
    half = 0.5 * (rw["z_te_mm"] - rw["z_le_mm"])
    for sec in rw["sections"]:
        cz = sec["chord_mm"] * math.cos(math.radians(sec["stagger_deg"]))
        ok_ext &= 0.5 * cz + 0.5 * sec["thickness_mm"] <= half + 1e-6
check("cuerda axial + ½t de toda sección cabe en su hueco", ok_ext)
ann = contract["annulus"]
check("annulus: z monótono y r_tip > r_hub",
      all(ann["tip"][i][0] < ann["tip"][i + 1][0]
          for i in range(len(ann["tip"]) - 1))
      and all(t[1] > h[1] for t, h in zip(ann["tip"], ann["hub"])))

# ==========================================================================
print("— T8 · estructura (solver de disco vs Timoshenko)")
r = np.linspace(0.02, 0.10, 240)
rho_m, E_m, nu_m, om = 4430.0, 113.8e9, 0.342, 1200.0
u_d, s_r, s_t = sc.solve_disk(r, np.full_like(r, 0.01), rho_m, E_m, nu_m, om)
s_r_ex, s_t_ex = sc.analytic_annular_disk(r, r[0], r[-1], rho_m, nu_m, om)
err_t = np.max(np.abs(s_t - s_t_ex)) / np.max(np.abs(s_t_ex))
check("disco anular h const: σ_θ vs Timoshenko < 1%", err_t < 0.01,
      f"{err_t:.2%}")
st1 = sc.evaluate_structural(THETA_REF, rec)
th_fast = THETA_REF.copy(); th_fast[1] *= 1.3
st2 = sc.evaluate_structural(th_fast, pc._meanline(th_fast))
check("esfuerzos crecen con RPM",
      st2["sigma_vm_max_MPa"] > st1["sigma_vm_max_MPa"])
check("g_struct tiene 4 componentes finitas",
      len(st1["g_struct"]) == 4
      and all(np.isfinite(g) for g in st1["g_struct"]))

# ==========================================================================
print("— T9 · pipeline de datos y calibración")
checks_dp = dp.validate()
check("data_pipeline.validate: todo OK",
      all(c.ok for c in checks_dp),
      f"{sum(c.ok for c in checks_dp)}/{len(checks_dp)}")
cal = pc.HiFiCalibration()
for xm, yh in [(2.0, 2.2), (3.0, 3.25), (4.0, 4.3)]:
    cal.register(dict(PR=xm, eta_poly=0.9), dict(PR=yh, eta_poly=0.88))
a, b = cal.coef["PR"]
check("HiFiCalibration recupera la afín (a≈1.05)",
      abs(a - 1.05) < 0.01, f"a={a:.4f} b={b:.4f}")

# ==========================================================================
print("— T10 · optimizador (núcleo portado)")
F = np.array([[1.0, 2.0], [2.0, 1.0], [1.5, 1.5], [3.0, 3.0]])
V = np.array([0.0, 0.0, 0.0, 0.0])
fronts = no._fast_nondominated_sort(F, V)
check("non-dominated sort: frente 0 = {0,1,2}, dominado = {3}",
      sorted(fronts[0]) == [0, 1, 2] and fronts[1] == [3])
check("dominancia de Deb: factible >> infactible",
      no._dominates(np.array([9, 9]), 0.0, np.array([0, 0]), 1.0))

rng = np.random.default_rng(3)
Xs = rng.random((80, 4))
Ys = np.column_stack([
    1.5 + Xs[:, 0],                 # 'PR'
    0.8 + 0.1 * Xs[:, 1],           # 'eta'
    300 + 100 * Xs[:, 2],           # 'U2'
    1.0 + 0.5 * Xs[:, 3],           # 'Mu'
    1e6 * (1 + Xs[:, 0]),           # 'power'
])
sur = no.EnsembleSurrogate(K=3)
sur.fit(Xs, Ys, log=lambda *a: None, rand=rng)
check("EnsembleSurrogate: R²(PR) > 0.9 en función suave",
      sur.metrics["r2"]["PR"] > 0.9, f"{sur.metrics['r2']['PR']:.3f}")

spec = no.DesignSpec(PR_target=4.0, massflow=25.0, material=None)
g_all = spec.constraints(rec, THETA_REF)
check("constraints: aero(8) + espec(4) sin material",
      len(g_all) == pc.N_CONSTRAINTS + 4, f"{len(g_all)}")
spec_m = no.DesignSpec(PR_target=4.0, massflow=25.0)
g_m = spec_m.constraints(rec, THETA_REF)
check("constraints: + g_struct(4) con material",
      len(g_m) == pc.N_CONSTRAINTS + 4 + sc.N_STRUCT_CONSTRAINTS)

with tempfile.TemporaryDirectory() as td:
    d = no.AutonomousAxialDesigner(spec, fidelity=pc.Fidelity.L0,
                                   log=lambda *a: None)
    d._eval_and_store(np.array([THETA_REF]))
    ck = os.path.join(td, "ck.json")
    d.save(ck)
    data = json.load(open(ck))
    check("checkpoint: dataset con theta y OUT_KEYS",
          len(data["dataset"]) == 1
          and all(k in data["dataset"][0] for k in no.OUT_KEYS))

# ==========================================================================
print("— T11 · caché e invariantes de evaluate")
r1 = pc.evaluate(THETA_REF, fidelity=pc.Fidelity.L0, use_cache=True,
                 calibrate=False)
r2 = pc.evaluate(THETA_REF, fidelity=pc.Fidelity.L0, use_cache=True,
                 calibrate=False)
check("caché: misma evaluación devuelve mismo PR", r1["PR"] == r2["PR"])
check("record trae claves del patrón Phy-CC (U2/Mu/T02 alias)",
      all(k in r1 for k in ("U2", "Mu", "T02", "g", "feasible")))

# ==========================================================================
print(f"\n{_n_pass}/{_n_pass + _n_fail} checks OK"
      + (f" — {_n_fail} FALLOS" if _n_fail else ""))
sys.exit(1 if _n_fail else 0)
