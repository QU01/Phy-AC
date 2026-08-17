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
import re
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
import contract_schema as cx

THETA_REF = np.array([4.0, 12_500.0, 0.62, 0.55, 0.30, -0.10,
                      0.60, 1.20, 1.10, 2.20, 288.15, 101_325.0, 25.0])

_n_pass = _n_fail = 0

# Los bloques que dependen de un extra opcional (CadQuery para todo el CAD
# detallado, turbo-design para la fidelidad L1) se saltan solos cuando la
# librería no está. En una máquina de desarrollo es razonable; en CI es
# exactamente lo contrario de lo que queremos, porque un job sin CadQuery
# pasa en verde sin haber comprobado NADA de la geometría. Con
# PHYAC_REQUIRE_STEP=1 la ausencia deja de ser un
# SKIP y pasa a ser un FALLO. El CI los pone a 1 en el job que instala los
# extras (ver .github/workflows/ci.yml).
REQUIRE_STEP = os.environ.get("PHYAC_REQUIRE_STEP", "") not in ("", "0")
# PHYAC_REQUIRE_L1 desapareció en la fase 11: L1 es propio y no puede
# faltar. Se conserva PHYAC_REQUIRE_STEP (CadQuery) y PHYAC_REQUIRE_STL
# (capa 5c compilada), que sí son opcionales de verdad.
# Igual para la capa 5c (C#): la paridad STL↔STEP solo se puede comprobar
# con la solucion compilada. PHYAC_REQUIRE_STL=1 en el job que la compila.
REQUIRE_STL = os.environ.get("PHYAC_REQUIRE_STL", "") not in ("", "0")


def _ok_call(fn) -> bool:
    """True si `fn()` no lanza. Para comprobar que algo FALLA de verdad."""
    try:
        fn()
        return True
    except Exception:
        return False


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
check("NDIM = 15 (13 legacy + phi_slope/Rx_slope, fase 8)",
      pc.NDIM == 15 and pc.NDIM_LEGACY == 13
      and pc.VAR_NAMES[13] == "phi_slope" and pc.VAR_NAMES[14] == "Rx_slope"
      and pc.VAR_NAMES[10] == "T0_in")   # el punto de operación no se movió

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
dh_h = pc.h_air(rec["T0_out"]) - pc.h_air(THETA_REF[10])
# Con gas caloríficamente imperfecto (fase 9) la identidad de Euler es
# ΣΔh0 = h(T0_out) − h(T0_in), NO cp·ΔT0 con un cp de referencia.
check("Euler: ΣΔh0 = h(T0_out) − h(T0_in) (±0.1%)",
      abs(dh_sum - dh_h) / dh_h < 1e-3,
      f"{dh_sum:.0f} vs {dh_h:.0f} J/kg")
_cp_eff = dh_sum / (rec["T0_out"] - THETA_REF[10])
check("cp efectivo por encima del de referencia (gas imperfecto)",
      pc.CP < _cp_eff < 1.10 * pc.CP, f"{_cp_eff:.1f} vs {pc.CP:.1f} J/kgK")

pr_prod = float(np.prod([s["PR"] for s in rec["stage_table"]]))
check("PR etapas ≈ PR máquina (antes del OGV, ±2%)",
      abs(pr_prod - rec["PR"]) / rec["PR"] < 0.02,
      f"{pr_prod:.3f} vs {rec['PR']:.3f}")

# límite sin pérdidas → η→1
# EW_A = 0 apaga también el débito de pared de Koch & Smith (fase 9):
# sin él la etapa nunca llega a η = 1 aunque las ω̄ sean cero.
_k = (pc.K_PROFILE, pc.K_ENDWALL, pc.K_SHOCK, pc.K_TIP_CLEARANCE,
      pc.EW_A, pc.TIP_CLEARANCE_MM)
pc.K_PROFILE = pc.K_ENDWALL = pc.K_SHOCK = pc.K_TIP_CLEARANCE = 0.0
pc.EW_A = pc.TIP_CLEARANCE_MM = 0.0
rec_is = pc._meanline(THETA_REF)
(pc.K_PROFILE, pc.K_ENDWALL, pc.K_SHOCK, pc.K_TIP_CLEARANCE,
 pc.EW_A, pc.TIP_CLEARANCE_MM) = _k
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

# corrección de Reynolds: f_Re continua, monótona y activa en máquinas
# pequeñas (Koch & Smith 1976 / Wassell 1968 / Schäffler 1980)
f_hi = pc._re_correction(pc.RE_REF)
f_mid_a = pc._re_correction(pc.RE_LAM * 1.0001)
f_mid_b = pc._re_correction(pc.RE_LAM * 0.9999)
check("f_Re: continua en RE_LAM/RE_REF y sin crédito sobre RE_REF",
      abs(f_hi - 1.0) < 1e-12 and abs(f_mid_a - f_mid_b) < 1e-3
      and pc._re_correction(2 * pc.RE_REF) == 1.0)
check("f_Re: monótona decreciente con Re",
      pc._re_correction(5e4) > pc._re_correction(3e5)
      > pc._re_correction(9e5) > 1.0 - 1e-12)
th_small = THETA_REF.copy()
th_small[12] = 2.5                  # máquina 10× más pequeña ⇒ Re_c ↓
rec_small = pc._meanline(th_small)
check("máquina pequeña (Re↓) ⇒ η↓ (corrección de Reynolds activa)",
      rec_small["eta_poly"] < rec["eta_poly"] - 0.002
      and any(s["losses"]["f_re_rotor"] > 1.0
              for s in rec_small["stage_table"]),
      f"η {rec_small['eta_poly']:.4f} vs {rec['eta_poly']:.4f}")

# holgura de punta por fila: ε absoluta en mm ⇒ ε/h crece hacia atrás
eh = [s["losses"]["eps_over_h"] for s in rec["stage_table"]]
fr = [s["losses"]["clearance"] / max(s["dh0"], 1e-3)
      for s in rec["stage_table"]]
check("holgura: ε/h crece hacia atrás y la última etapa pierde más",
      all(eh[i] < eh[i + 1] for i in range(len(eh) - 1))
      and fr[-1] > fr[0],
      f"ε/h {eh[0]:.4f}→{eh[-1]:.4f}, frac {fr[0]:.4f}→{fr[-1]:.4f}")
sweep = np.linspace(0.0, 0.06, 121)
vals = [pc._clearance_loss_frac(e) for e in sweep]
check("_clearance_loss_frac: continua y monótona en [0, 0.06]",
      all(vals[i + 1] >= vals[i] for i in range(len(vals) - 1))
      and max(abs(vals[i + 1] - vals[i]) for i in range(len(vals) - 1))
      < 2.5 * pc.K_TIP_CLEARANCE * (sweep[1] - sweep[0]))
_eps = pc.TIP_CLEARANCE_MM
pc.TIP_CLEARANCE_MM = 3.0 * _eps
rec_open = pc._meanline(THETA_REF)
pc.TIP_CLEARANCE_MM = _eps
check("holgura: ε↑ ⇒ η↓ (sensibilidad restaurable por parámetro)",
      rec_open["eta_poly"] < rec["eta_poly"] - 0.002,
      f"η {rec_open['eta_poly']:.4f} vs {rec['eta_poly']:.4f}")

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

# bucket de incidencia W(M): se estrecha con Mach, rama negativa tolerante
check("bucket: se estrecha con M y es 1 en incidencia nula",
      pc._incidence_bucket(5.0, 0.9) > pc._incidence_bucket(5.0, 0.3) > 1.0
      and pc._incidence_bucket(0.0, 0.9) == 1.0
      and pc._incidence_bucket(60.0, 0.9) == 4.0)
check("bucket: rama negativa 1.5× más tolerante",
      pc._incidence_bucket(-5.0, 0.9) < pc._incidence_bucket(5.0, 0.9))
# desviación off-design (Creveling): con incidencia positiva la ψ lograda
# cae respecto al modelo congelado puro (sub-giro creciente hacia stall)
od_ns = pc.offdesign(THETA_REF, rpm=float(THETA_REF[1]), mdot=21.0)
_k_dev = pc.K_DEV_OFFDESIGN
pc.K_DEV_OFFDESIGN = 0.0
od_ns0 = pc.offdesign(THETA_REF, rpm=float(THETA_REF[1]), mdot=21.0)
pc.K_DEV_OFFDESIGN = _k_dev
check("desviación off-design: ψ y PR caen cerca de stall (i>0)",
      od_ns["stage_table"][0]["incidence_deg"] > 2.0
      and od_ns["stage_table"][0]["psi"] < od_ns0["stage_table"][0]["psi"]
      and od_ns["PR"] < od_ns0["PR"],
      f"ψ0 {od_ns['stage_table'][0]['psi']:.3f} vs "
      f"{od_ns0['stage_table'][0]['psi']:.3f}")
# VSV: a 0.7N el schedule auto descarga las etapas frontales; a 1.0N nada
mp_fix = pc.compressor_map(THETA_REF, speed_fracs=(0.7, 1.0), n_points=11)
mp_vsv = pc.compressor_map(THETA_REF, speed_fracs=(0.7, 1.0), n_points=11,
                           vsv="auto")
_stall_fix = sum(p["stall"] for p in mp_fix["speedlines"][0]["points"])
_stall_vsv = sum(p["stall"] for p in mp_vsv["speedlines"][0]["points"])
check("VSV auto: menos puntos en stall a 0.7N que geometría fija",
      mp_vsv["speedlines"][0]["vsv_deg"] > 0.0
      and _stall_vsv < _stall_fix,
      f"{_stall_vsv} vs {_stall_fix} (VSV "
      f"{mp_vsv['speedlines'][0]['vsv_deg']:.0f}°)")
check("VSV auto: speedline de diseño idéntica (VSV=0 a N=Nd)",
      mp_vsv["speedlines"][1]["vsv_deg"] == 0.0
      and all(abs(a["PR"] - b["PR"]) < 1e-12
              for a, b in zip(mp_fix["speedlines"][1]["points"],
                              mp_vsv["speedlines"][1]["points"])))

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
check("schema phyac-axial-2", contract["schema"] == "phyac-axial-2")
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
# IGV/OGV: las filas que la física asume ahora existen en el contrato
igv, ogv = contract["igv"], contract["ogv"]
a1_st0 = rec["stage_table"][0]["alpha1_deg"]
a1_last = rec["stage_table"][-1]["alpha1_deg"]
mid_igv = igv["sections"][len(igv["sections"]) // 2]
mid_ogv = ogv["sections"][len(ogv["sections"]) // 2]
check("IGV: entrada axial y salida sobre-girada más allá de α1 (Carter)",
      abs(mid_igv["metal_in_deg"]) < 3.0
      and abs(mid_igv["metal_out_deg"]) > a1_st0 - 1.0,
      f"χ2={mid_igv['metal_out_deg']:.1f}° vs α1={a1_st0:.1f}°")
check("OGV: entrada ≈ α1 residual y salida ≈ axial (± desviación)",
      abs(abs(mid_ogv["metal_in_deg"]) - a1_last) < 6.0
      and abs(mid_ogv["metal_out_deg"]) < 12.0,
      f"χ1={mid_ogv['metal_in_deg']:.1f}° vs α1={a1_last:.1f}°")
check("IGV/OGV: mismas invariantes de sección que las etapas",
      all(len(s["points"]) == bp.N_POINTS for s in igv["sections"])
      and len({len(s["camber_points"])
               for s in igv["sections"] + ogv["sections"]}) == 1
      and not igv["rotating"] and not ogv["rotating"])

# sin solape axial entre filas consecutivas (la extensión del contrato ya
# incluye la inflación del shell ±½t — regresión del defecto visto en STL)
row_seq = [contract["igv"]]
for st in contract["stages"]:
    row_seq += [st["rotor"], st["stator"]]
row_seq.append(contract["ogv"])
check("filas sin solape axial (IGV→etapas→OGV, incl. inflación del shell)",
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
check("g_struct tiene 5 componentes finitas (+ Campbell, fase 7)",
      len(st1["g_struct"]) == sc.N_STRUCT_CONSTRAINTS == 5
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
check("constraints: aero(9) + espec(5) sin material",
      len(g_all) == pc.N_CONSTRAINTS + 5, f"{len(g_all)}")
spec_m = no.DesignSpec(PR_target=4.0, massflow=25.0)
g_m = spec_m.constraints(rec, THETA_REF)
check("constraints: + g_struct(5) con material",
      len(g_m) == pc.N_CONSTRAINTS + 5 + sc.N_STRUCT_CONSTRAINTS)

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
print("— T12 · controlabilidad (fixed vars, eval directa, warm start)")
spec_fx = no.DesignSpec(PR_target=4.0, massflow=25.0, material=None,
                        fixed_vars={"n_stages": 5, "phi1": 0.60})
th_fx = spec_fx.fix_operating_point(THETA_REF)
check("fixed_vars: fija n_stages y phi1 en TODA evaluación",
      th_fx[0] == 5.0 and th_fx[3] == 0.60
      and th_fx[12] == spec_fx.massflow)
try:
    no.DesignSpec(material=None, fixed_vars={"massflow": 10.0})
    ok_reject = False
except ValueError:
    ok_reject = True
try:
    no.DesignSpec(material=None, fixed_vars={"no_existe": 1.0})
    ok_reject2 = False
except ValueError:
    ok_reject2 = True
check("fixed_vars: rechaza punto de operación y nombres inválidos",
      ok_reject and ok_reject2)

res_ev = no.evaluate_design(spec_fx, THETA_REF[:10],
                            fidelity=pc.Fidelity.L0)
check("evaluate_design: 10 valores + spec ⇒ record completo y front de 1",
      res_ev["record"]["PR"] > 1.0 and len(res_ev["pareto_front"]) == 1
      and res_ev["theta"][0] == 5.0,        # fixed_vars también aplican
      f"PR={res_ev['record']['PR']:.3f}")

with tempfile.TemporaryDirectory() as td:
    d1 = no.AutonomousAxialDesigner(spec, fidelity=pc.Fidelity.L0,
                                    log=lambda *a: None, seed=42)
    thetas3 = np.array([THETA_REF,
                        THETA_REF * [1, 1.05, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        THETA_REF * [1, 0.95, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])
    d1._eval_and_store(thetas3)
    ck2 = os.path.join(td, "warm.json")
    d1.save(ck2)
    d2 = no.AutonomousAxialDesigner.load(ck2, fidelity=pc.Fidelity.L0,
                                         log=lambda *a: None)
    check("warm start: load() restaura n_evals, semilla y salidas",
          len(d2.recs) == 3 and d2.seed == 42
          and abs(d2.recs[0]["PR"] - d1.recs[0]["PR"]) < 1e-9
          and d2.recs[0].get("stage_table") is not None)
    sp2, front2 = no.pareto_from_checkpoint(ck2)
    check("pareto_from_checkpoint: front verificado con θ completo",
          len(front2) >= 1
          and all(len(q["theta"]) == pc.NDIM for q in front2))

import phyac_cli
fx = phyac_cli.parse_fixed(["phi1=0.6", "n_stages=5"])
check("CLI parse_fixed: VAR=VALOR → dict numérico",
      fx == {"phi1": 0.6, "n_stages": 5.0})
a_cli = phyac_cli.parse_args(["--pr", "4", "--fix", "phi1=0.6",
                              "--seed", "7", "--eval-theta", "1,2,3"])
check("CLI: --fix/--seed/--eval-theta parsean",
      a_cli.fix == ["phi1=0.6"] and a_cli.seed == 7
      and a_cli.eval_theta == "1,2,3")

# ==========================================================================
print("— T13 · dinámica de álabes (fase 7: Southwell, Campbell, K_t)")
# f1 estática vs fórmula analítica de viga uniforme (misma sección rect)
_stage0 = rec["stage_table"][0]
_md = sc.blade_modes(_stage0, "Ti-6Al-4V", RPM=0.0, kind="rotor")
_c = _stage0["chord_rotor_mm"] * 1e-3
_t = 0.10 * _c                                # TC_ROOT_R
_h = _stage0["h_blade_mm"] * 1e-3
_E, _rho = sc.MATERIALS["Ti-6Al-4V"]["E"], sc.MATERIALS["Ti-6Al-4V"]["rho"]
_f_exact = (sc.LAMBDA1_SQ / (2 * math.pi)) * math.sqrt(
    _E * (_c * _t ** 3 / 12) / (_rho * _c * _t * _h ** 4))
check("blade_modes: f1 estática = viga empotrada analítica (<2%)",
      abs(_md["f1_static_Hz"] - _f_exact) / _f_exact < 0.02,
      f"{_md['f1_static_Hz']:.1f} vs {_f_exact:.1f} Hz")
_f_lo = sc.blade_modes(_stage0, "Ti-6Al-4V", 8000.0)["f1_Hz"]
_f_hi = sc.blade_modes(_stage0, "Ti-6Al-4V", 16000.0)["f1_Hz"]
check("blade_modes: f1 crece con RPM (Southwell) y > f1 estática",
      _f_hi > _f_lo > _md["f1_static_Hz"])
_s_short = dict(_stage0, h_blade_mm=0.5 * _stage0["h_blade_mm"])
check("blade_modes: f1 cae con la altura de álabe (∝1/h²)",
      sc.blade_modes(_s_short, "Ti-6Al-4V", 12500.0)["f1_static_Hz"]
      > 3.0 * _md["f1_static_Hz"])
_kts = [sc.k_t_root(r, 4.8) for r in (1.0, 2.0, 3.0, 4.0)]
check("k_t_root: monótono decreciente con r_fillet y en banda [1.3, 2.5]",
      all(_kts[i] > _kts[i + 1] for i in range(3))
      and all(1.3 <= k <= 2.5 for k in _kts),
      f"K_t={[round(k, 2) for k in _kts]}")
_marg = [sc.evaluate_structural(THETA_REF, rec)["campbell_margin_min"]]
_st_r = sc.evaluate_structural(THETA_REF, rec)
check("evaluate_structural: Campbell/Goodman/flutter expuestos y finitos",
      np.isfinite(_st_r["campbell_margin_min"])
      and _st_r["sigma_alt_allow_MPa"] > 0
      and len(_st_r["flutter"]) == 2 * len(rec["stage_table"])
      and _st_r["k_t_root"] > 1.0)
_prev_m = None
_ok_cont = True
for _rpm_f in np.linspace(0.9, 1.1, 9):
    _th_c = THETA_REF.copy()
    _th_c[1] *= _rpm_f
    _m = sc.evaluate_structural(_th_c, pc._meanline(_th_c))[
        "campbell_margin_min"]
    if _prev_m is not None:
        _ok_cont &= abs(_m - _prev_m) < 0.35
    _prev_m = _m
check("margen de Campbell continuo con RPM (sin saltos > 0.35)", _ok_cont)

# ==========================================================================
print("— T14 · θ 15-D: padding, pendientes por etapa y per_stage (fase 8)")
_t15 = pc.pad_theta(THETA_REF)
_r13 = pc._meanline(THETA_REF)
_r15 = pc._meanline(_t15)
check("pad_theta: θ 13-D ≡ 15-D con pendientes 0 (bit-exacto)",
      len(_t15) == 15 and _t15[13] == 0.0 and _t15[14] == 0.0
      and abs(_r13["PR"] - _r15["PR"]) < 1e-12
      and abs(_r13["eta_poly"] - _r15["eta_poly"]) < 1e-12
      and max(abs(a - b) for a, b in zip(_r13["g"], _r15["g"])) < 1e-12)

_ts = _t15.copy(); _ts[13] = -0.15; _ts[14] = 0.10
_rs = pc._meanline(_ts)
_phis = [s["phi"] for s in _rs["stage_table"]]
_rxs = [s["Rx"] for s in _rs["stage_table"]]
check("pendientes activas: φ cae y Rx sube hacia atrás (ξ lineal)",
      all(_phis[i] > _phis[i + 1] for i in range(len(_phis) - 1))
      and all(_rxs[i] < _rxs[i + 1] for i in range(len(_rxs) - 1))
      and abs(_rs["PR"] - _r15["PR"]) > 1e-3,
      f"φ {_phis[0]:.3f}→{_phis[-1]:.3f}  Rx {_rxs[0]:.3f}→{_rxs[-1]:.3f}")

_n4 = _rs["n_stages"]
_ps = dict(phi=[float(THETA_REF[3])] * _n4, Rx=[float(THETA_REF[6])] * _n4)
_rp = pc._meanline(_t15, per_stage=_ps)
check("per_stage constante ≡ distribuciones escalares",
      abs(_rp["PR"] - _r15["PR"]) < 1e-12
      and abs(_rp["eta_poly"] - _r15["eta_poly"]) < 1e-12)

spec_s = no.DesignSpec(PR_target=4.0, massflow=25.0, material=None,
                       fixed_vars={"phi_slope": 0.1, "Rx_slope": -0.05})
th_s = spec_s.fix_operating_point(THETA_REF)      # legacy 13 → pad + fix
check("fixed_vars: pendientes fijables; fix_operating_point paddea",
      len(th_s) == 15 and th_s[13] == 0.1 and th_s[14] == -0.05
      and th_s[12] == spec_s.massflow)

res12 = no.evaluate_design(no.DesignSpec(material=None),
                           np.concatenate([THETA_REF[:10], [-0.15, 0.10]]),
                           fidelity=pc.Fidelity.L0)
check("evaluate_design: 12 valores (10 diseño + 2 pendientes)",
      len(res12["theta"]) == 15 and res12["theta"][13] == -0.15
      and res12["record"]["PR"] > 1.0)

# STEP de ensamble (fase 8.2) — condicional a la dependencia opcional
try:
    import cadquery as _cq_mod          # noqa: F401
    _has_cq = True
except Exception:
    _has_cq = False
if _has_cq:
    with tempfile.TemporaryDirectory() as td:
        _paths = gg.export_step(contract, td, mode="parts")
        _names = {os.path.basename(p) for p in _paths}
        check("export_step(parts): Shaft + piezas por etapa + README",
              "Shaft.step" in _names and "README.txt" in _names
              and any(n.startswith("RotorStage") for n in _names)
              and all(os.path.getsize(p) > 0 for p in _paths))
else:
    check("export_step: cadquery no instalado — degradación silenciosa "
          + ("(EXIGIDO por PHYAC_REQUIRE_STEP)" if REQUIRE_STEP
             else "(SKIP)"),
          not REQUIRE_STEP
          and gg.export_step(contract, ".", mode="parts") == [])

# ==========================================================================
print("— T15 · fase 9: gas real, entropía, Koch, endwall, sangrado, rango")

# --- gas caloríficamente imperfecto ---------------------------------------
check("cp(T) creciente y en banda de tablas (250-1000 K)",
      all(pc.cp_air(T) < pc.cp_air(T + 50) for T in (250, 400, 600, 800))
      and abs(pc.cp_air(288.15) - 1004.0) < 4.0
      and abs(pc.cp_air(700.0) - 1075.0) < 12.0,
      f"cp(288)={pc.cp_air(288.15):.1f}  cp(700)={pc.cp_air(700):.1f}")
check("γ(T) cae con T y vale 1.400 en 288 K",
      abs(pc.gamma_air(288.15) - 1.400) < 0.003
      and pc.gamma_air(800.0) < pc.gamma_air(288.15),
      f"γ(288)={pc.gamma_air(288.15):.4f}  γ(800)={pc.gamma_air(800):.4f}")
_Ts = [250.0, 400.0, 650.0, 900.0]
check("h(T) y phi(T) invierten exacto (Newton, <1e-6 K)",
      all(abs(pc.T_from_h(pc.h_air(T), 500.0) - T) < 1e-6 for T in _Ts)
      and all(abs(pc.T_from_phi(pc.phi_air(T), 500.0) - T) < 1e-6
              for T in _Ts))
_gv = pc.GAS_VARIABLE
pc.GAS_VARIABLE = False
check("modo gas perfecto colapsa a (CP, GAMMA) exactos",
      pc.cp_air(900.0) == pc.CP and pc.gamma_air(900.0) == pc.GAMMA
      and abs(pc.h_air(400.0) - pc.CP * 400.0) < 1e-9)
pc.GAS_VARIABLE = _gv

# --- Koch 1981 -------------------------------------------------------------
check("koch_L_over_g2 crece con comba y solidez, cae con β2",
      pc.koch_L_over_g2(30, 1.2, 35) > pc.koch_L_over_g2(10, 1.2, 35)
      and pc.koch_L_over_g2(20, 1.5, 35) > pc.koch_L_over_g2(20, 1.0, 35)
      and pc.koch_L_over_g2(20, 1.2, 55) > pc.koch_L_over_g2(20, 1.2, 25))
_chs = [pc.koch_ch_stall(x) for x in (0.7, 1.0, 2.0, 3.0)]
check("Ch_stall monótona en L/g₂ y en la banda medida de Koch",
      all(a < b for a, b in zip(_chs, _chs[1:]))
      and 0.25 < _chs[0] and _chs[-1] < 0.58,
      f"Ch_stall = {[round(c, 3) for c in _chs]}")
check("Ch_stall: Re bajo y holgura grande la degradan (Figs. 4 y 5)",
      pc.koch_ch_stall(1.5, re_c=2e4) < pc.koch_ch_stall(1.5, re_c=5e5)
      and pc.koch_ch_stall(1.5, eps_over_g=0.15)
      < pc.koch_ch_stall(1.5, eps_over_g=0.01))
check("𝔉_ef en la banda medida (0.78-1.11) para un triángulo típico",
      0.78 <= pc.koch_f_ef(20.0, 55.0, 300.0, 260.0) <= 1.11,
      f"{pc.koch_f_ef(20.0, 55.0, 300.0, 260.0):.3f}")

# --- endwall de Koch & Smith ----------------------------------------------
check("2δ*/g crece con la carga y con la holgura (ec. 3, Fig. 8)",
      pc.endwall_blockage(0.95, 0.02) > pc.endwall_blockage(0.70, 0.02)
      and pc.endwall_blockage(0.85, 0.10) > pc.endwall_blockage(0.85, 0.01)
      and pc.endwall_blockage(0.0, 0.0) == 0.0)
_r9 = pc.evaluate(THETA_REF, fidelity=pc.Fidelity.L0, use_cache=False)
check("bloqueo EMERGE de la carga: KB baja y la etapa cargada bloquea más",
      all(s["blockage_out"] < 1.0 for s in _r9["stage_table"])
      and _r9["stage_table"][0]["blockage_out"]
      < _r9["stage_table"][-1]["blockage_out"],
      f"KB {_r9['stage_table'][0]['blockage_out']:.3f} → "
      f"{_r9['stage_table'][-1]['blockage_out']:.3f}")
check("cada etapa reporta Ch, Ch_stall, L/g₂, 𝔉_ef y ds finitos",
      all(np.isfinite([s["Ch"], s["Ch_stall"], s["L_over_g2"], s["f_ef"],
                       s["ds_stage_J_kgK"]]).all()
          for s in _r9["stage_table"]))

# --- entropía: coherencia marcha ↔ rendimiento -----------------------------
_ds_sum = sum(s["ds_stage_J_kgK"] for s in _r9["stage_table"])
check("Σds de etapas ≤ ds de máquina (IGV/OGV añaden, nunca restan)",
      _ds_sum <= _r9["ds_machine_J_kgK"] + 1e-9,
      f"{_ds_sum:.2f} vs {_r9['ds_machine_J_kgK']:.2f} J/kgK")
check("η_poly coincide con R·lnPR/Δphi (definición de gas imperfecto)",
      abs(_r9["eta_poly"] - pc.RGAS * math.log(_r9["PR"])
          / (pc.phi_air(_r9["T0_out"]) - pc.phi_air(THETA_REF[10]))) < 1e-9)

# --- IGV y sangrado --------------------------------------------------------
check("el IGV existe en la física y paga pérdida y longitud",
      _r9["igv"]["present"] and _r9["igv"]["dP0_Pa"] > 0.0
      and _r9["igv"]["length_mm"] > 0.0,
      f"ΔP0_IGV = {_r9['igv']['dP0_Pa']:.0f} Pa")
_r_bl = pc._meanline(THETA_REF, bleed={1: 0.05})
check("sangrado: sale menos gasto y menos potencia que sin sangrar",
      _r_bl["bleed"]["mdot_exit_kgs"] < 0.96 * THETA_REF[12]
      and _r_bl["power_W"] < _r9["power_W"],
      f"ṁ_salida {_r_bl['bleed']['mdot_exit_kgs']:.2f} kg/s")
check("sin sangrado el gasto se conserva exactamente",
      abs(_r9["bleed"]["total_frac"]) < 1e-12)

# --- solidez real vs optimizada -------------------------------------------
check("σ real (nº entero de álabes) reportada y cercana a la de θ",
      all(abs(s["sigma_rotor_actual"] - THETA_REF[7]) < 0.25
          for s in _r9["stage_table"]),
      f"σ_θ={THETA_REF[7]:.2f} σ_real="
      f"{_r9['stage_table'][0]['sigma_rotor_actual']:.3f}")

# --- rango operativo: el hallazgo que motivó la fase 9 ---------------------
check("g tiene 9 componentes y la 9ª es la reacción de raíz",
      len(_r9["g"]) == pc.N_CONSTRAINTS == 9
      and abs(_r9["g"][8] - (pc.RX_HUB_MIN - _r9["Rx_hub_min"])
              / pc.RX_HUB_MIN) < 1e-9)
_sm = pc.surge_margin(THETA_REF, _r9)
check("surge_margin devuelve un punto de trabajo con margen positivo",
      _sm["ok"] and 0.0 < _sm["sm_flow"] < 1.0
      and _sm["mdot_surge"] < _sm["mdot_wl"] <= _sm["mdot_choke"],
      f"SM_gasto = {100 * _sm['sm_flow']:.1f}%  "
      f"(surge {_sm['mdot_surge']:.1f} < wl {_sm['mdot_wl']:.1f} "
      f"≤ choke {_sm['mdot_choke']:.1f} kg/s)")
check("línea de trabajo: ṁ crece con PR y pasa por el diseño",
      abs(pc.working_line_mdot(_r9["PR"], _r9["PR"], 25.0,
                               _r9["eta_poly"]) - 25.0) < 1e-9
      and pc.working_line_mdot(1.2 * _r9["PR"], _r9["PR"], 25.0, 0.9) > 25.0)
_map9 = pc.compressor_map(THETA_REF, speed_fracs=(0.8, 1.0))
check("el mapa tiene ANCHO real en todas las speedlines (regresión H1)",
      all(sl["mdot_choke"] > 1.05 * sl["mdot_surge"] for sl in
          _map9["speedlines"]),
      "  ".join(f"N={sl['speed_frac']:.2f}:"
                f"{100 * (sl['mdot_choke'] / sl['mdot_surge'] - 1):.0f}%"
                for sl in _map9["speedlines"]))
check("el choke LIMITA el gasto: ningún punto del mapa tiene PR < 1",
      all(p["PR"] >= 1.0 for sl in _map9["speedlines"]
          for p in sl["points"]))
check("cada speedline trae punto de trabajo y etapa crítica",
      all(sl["working_point"] is not None
          and 0 <= sl["working_point"]["stage_stall_first"]
          < len(_r9["stage_table"]) for sl in _map9["speedlines"]))
# Los VSV solo tienen sentido en la máquina que los NECESITA: una de PR
# bajo ya tiene margen de sobra a velocidad parcial. Se prueba en una de
# PR ≈ 4.5 y 6 etapas, que es donde vive el problema histórico (J79).
_TH_HI = np.array([6.0, 14_000.0, 0.55, 0.55, 0.34, -0.10, 0.62,
                   1.25, 1.15, 2.00, 288.15, 101_325.0, 25.0])
_hi_no = pc.compressor_map(_TH_HI, speed_fracs=(0.70, 0.85), vsv="none")
_hi_yes = pc.compressor_map(_TH_HI, speed_fracs=(0.70, 0.85), vsv="auto")
check("VSV: cerrar los estátores frontales corre el bombeo a menos gasto",
      _hi_yes["speedlines"][0]["mdot_surge"]
      < 0.95 * _hi_no["speedlines"][0]["mdot_surge"],
      f"ṁ_surge 70%: {_hi_no['speedlines'][0]['mdot_surge']:.2f} → "
      f"{_hi_yes['speedlines'][0]['mdot_surge']:.2f} kg/s")
check("VSV: y con ello mejoran el margen en la línea de trabajo (85%)",
      _hi_yes["speedlines"][1]["working_point"]["sm_flow"]
      > _hi_no["speedlines"][1]["working_point"]["sm_flow"],
      f"SM 85%: "
      f"{100 * _hi_no['speedlines'][1]['working_point']['sm_flow']:.1f}% → "
      f"{100 * _hi_yes['speedlines'][1]['working_point']['sm_flow']:.1f}%")
check("sin VSV la máquina de PR alto se queda SIN margen al 70% "
      "(la razón histórica de los estátores variables)",
      _hi_no["speedlines"][0]["working_point"]["sm_flow"] < 0.05,
      f"SM 70% sin VSV = "
      f"{100 * _hi_no['speedlines'][0]['working_point']['sm_flow']:.1f}%")

# --- retención de abeto (fase 9, paridad con turbodesigner) ----------------
_ft = gg.firtree_profile(6.0, 18.0, n_teeth=3, taper_deg=8.0)
_fu = [u for u, v in _ft]
_fv = [v for u, v in _ft]
check("firtree: contorno cerrado, simétrico y dentro de la envolvente",
      len(_ft) >= 12 and abs(max(_fu) + min(_fu)) < 1e-9
      and abs(max(_fu) - 3.0) < 1e-9 and abs(max(_fv) - 18.0) < 1e-9,
      f"{len(_ft)} puntos, u ±{max(_fu):.2f} mm, v hasta {max(_fv):.1f} mm")
_right = [(u, v) for u, v in _ft if u > 0]        # solo el flanco +u
_ru = [u for u, v in _right]
_w_top = _ru[0]
_w_bot = _ru[-1]
_necks = sum(1 for a, b, c in zip(_ru, _ru[1:], _ru[2:]) if b < a and b < c)
check("firtree: el árbol se estrecha hacia el fondo y tiene 3 cuellos",
      _w_bot < _w_top and _necks >= 2,
      f"semiancho {_w_top:.2f} → {_w_bot:.2f} mm, {_necks} cuellos")
_ft_slot = gg.firtree_profile(6.0, 18.0, clearance_mm=0.08)
check("firtree: la ranura del disco es MAYOR que la raíz (juego de montaje)",
      max(u for u, v in _ft_slot) > max(_fu),
      f"ranura ±{max(u for u, v in _ft_slot):.3f} vs raíz ±{max(_fu):.3f} mm")
check("contrato: bloque firtree presente y coherente con el álabe",
      "firtree" in contract["assembly"]
      and 4.0 <= contract["assembly"]["firtree"]["depth_mm"]
      < contract["stages"][0]["rotor"]["sections"][0]["r_mm"]
      and contract["assembly"]["firtree"]["n_teeth"] >= 2,
      f"depth = {contract['assembly']['firtree']['depth_mm']:.1f} mm")
check("contrato: el puerto de sangrado apunta detrás de la etapa crítica",
      contract["assembly"]["bleed_stage"]
      == min(max(rec["stage_stall_first"] + 2, 1), rec["n_stages"]),
      f"bleed_stage = {contract['assembly']['bleed_stage']} "
      f"(1ª en stall: etapa {rec['stage_stall_first'] + 1})")

# ==========================================================================
print("— T16 · fase 9.1: torbellino controlado y ensamble detallado")

# --- exponente de torbellino ----------------------------------------------
check("vórtice libre (n=−1) es exactamente r·Cu = cte y Cx constante",
      all(abs(pc.vortex_cu(100.0, r, 1.0, -1.0) * r - 100.0) < 1e-12
          for r in (0.6, 1.0, 1.4))
      and all(pc.vortex_cx(150.0, 90.0, r, 1.0, -1.0) == 150.0
              for r in (0.6, 1.0, 1.4)))
check("vórtice forzado (n=+1) da reacción CONSTANTE en el span",
      all(abs(pc.vortex_reaction(0.6, r, 1.0, 1.0) - 0.6) < 1e-12
          for r in (0.5, 1.0, 1.5)))
check("n=−1 reproduce la fórmula clásica Rx_h = 1−(1−Rx_m)(r_m/r)²",
      abs(pc.vortex_reaction(0.6, 0.8, 1.0, -1.0)
          - (1.0 - 0.4 / 0.8 ** 2)) < 1e-12)
check("Cx(r) continua al cruzar n = 0 (rama logarítmica del SRE)",
      abs(pc.vortex_cx(150.0, 90.0, 1.3, 1.0, -1e-5)
          - pc.vortex_cx(150.0, 90.0, 1.3, 1.0, +1e-5)) < 1e-3)
_rx = [pc.vortex_reaction(0.6, 0.7, 1.0, n) for n in (-1.0, -0.5, 0.0)]
check("subir n rescata la reacción de RAÍZ (lo que un fan necesita)",
      _rx[0] < _rx[1] < _rx[2],
      f"Rx_hub(n=−1,−0.5,0) = {[round(x, 3) for x in _rx]}")

_vn = pc.VORTEX_N
_r_free = pc.evaluate(THETA_REF, fidelity=pc.Fidelity.L0)
pc.VORTEX_N = -0.4
_r_ctrl = pc.evaluate(THETA_REF, fidelity=pc.Fidelity.L0)
pc.VORTEX_N = _vn
check("el caché distingue configuraciones de módulo (no devuelve la otra)",
      _r_ctrl["Rx_hub_min"] > _r_free["Rx_hub_min"] + 0.01,
      f"Rx_hub {_r_free['Rx_hub_min']:.3f} → {_r_ctrl['Rx_hub_min']:.3f}")
check("con n=−1 el record es BIT-EXACTO al de antes de la fase 9.1",
      pc.evaluate(THETA_REF, fidelity=pc.Fidelity.L0)["PR"]
      == _r_free["PR"])

# --- IGV solo si hay pre-swirl --------------------------------------------
_th_fan = np.array([1.0, 6200.0, 0.36, 0.55, 0.375, 0.0, 1.0 - 0.375 / 2,
                    1.35, 1.25, 2.7, 288.15, 101_325.0, 120.0])
_rec_fan = pc.evaluate(_th_fan, fidelity=pc.Fidelity.L0, use_cache=False)
_c_fan = gg.build_contract(_th_fan, _rec_fan, None, None)
check("Rx = 1−ψ/2 ⇒ α₁ = 0 exacto ⇒ la máquina NO lleva IGV",
      abs(_rec_fan["stage_table"][0]["alpha1_deg"]) < 1e-9
      and not _rec_fan["igv"]["present"] and _c_fan.get("igv") is None
      and _c_fan.get("ogv") is not None,
      f"α₁ = {_rec_fan['stage_table'][0]['alpha1_deg']:.2e}°")
check("con pre-swirl el IGV sí aparece en física y en contrato",
      rec["igv"]["present"] and contract.get("igv") is not None)

# --- ensamble detallado ----------------------------------------------------
if gg._cq() is None:
    check("ensamble detallado: cadquery no instalado "
          + ("(EXIGIDO por PHYAC_REQUIRE_STEP)" if REQUIRE_STEP
             else "(SKIP)"), not REQUIRE_STEP)
else:
    _cq = gg._cq()
    _asm = gg._cq_detailed_machine(_cq, contract)
    _subs = {ch.name: ch for ch in _asm.children}
    check("ensamble: tres subensambles ROTOR / STATOR / FASTENERS",
          set(_subs) == {"ROTOR", "STATOR", "FASTENERS"})
    _rnames = [ch.name for ch in _subs["ROTOR"].children]
    check("ROTOR: eje, un disco POR ETAPA y los álabes separados",
          "Shaft" in _rnames
          and sum(n.startswith("Disc") for n in _rnames) == rec["n_stages"]
          and sum(n.startswith("Rotor") for n in _rnames)
          == sum(s["rotor"]["n_blades"] for s in contract["stages"]),
          f"{len(_rnames)} componentes de rotor")
    _snames = [ch.name for ch in _subs["STATOR"].children]
    check("STATOR: un anillo de carcasa POR ETAPA de estátor",
          sum(n.startswith("CasingRing") for n in _snames)
          == rec["n_stages"])
    _fnames = [ch.name for ch in _subs["FASTENERS"].children]
    check("FASTENERS: tirantes del rotor y pernos en cada junta",
          sum(n.startswith("TieBolt") for n in _fnames)
          == contract["assembly"]["tie_bolt_count"]
          and sum(n.startswith("CasingBolt") for n in _fnames)
          == (rec["n_stages"] - 1)
          * contract["assembly"]["flange_bolt_count"],
          f"{len(_fnames)} elementos de unión")
    # el disco DEBE tener menos material que el mismo disco sin ranurar
    def _find(asm, name):
        for ch in asm.children:
            if ch.name == name:
                return ch
        return None

    _d1 = _find(_subs["ROTOR"], "Disc1")
    _c2 = json.loads(json.dumps(contract))
    _c2["assembly"]["firtree"]["enabled"] = False
    _asm0 = gg._cq_detailed_machine(_cq, _c2)
    _d0 = _find(_find(_asm0, "ROTOR"), "Disc1")
    check("el disco lleva las ranuras REALMENTE brochadas (menos volumen)",
          _d1.obj.val().Volume() < 0.995 * _d0.obj.val().Volume(),
          f"{_d1.obj.val().Volume() / 1000:.0f} vs "
          f"{_d0.obj.val().Volume() / 1000:.0f} cm³")
    # --- REGRESIÓN: ningún álabe puede atravesar carcasa ni tambor -------
    # Bug reportado sobre el ensamble: las palas se salían del casing por
    # fuera y del cubo por dentro. Tres causas encadenadas — secciones
    # planas sin compensar el radio, filas talladas con los radios de la
    # ENTRADA de su etapa (el annulus se contrae) y ninguna holgura de
    # marcha aplicada en el CAD. Máximos medidos entonces: +8.5 mm fuera
    # de la carcasa y −8.8 mm dentro del tambor.
    _hub_l, _tip_l = contract["annulus"]["hub"], contract["annulus"]["tip"]
    _eps = contract["annulus"]["tip_clearance_mm"]
    _zh = [p[0] for p in _hub_l]
    _rh = [p[1] for p in _hub_l]
    _zt = [p[0] for p in _tip_l]
    _rt = [p[1] for p in _tip_l]
    _worst_out = -9e9
    _worst_in = 9e9
    _rows = [(s["rotor"], "rotor") for s in contract["stages"]] \
        + [(s["stator"], "stator") for s in contract["stages"]] \
        + [(contract[k], "stator") for k in ("igv", "ogv") if contract.get(k)]
    for _row, _kind in _rows:
        _sh = gg._cq_blade_trimmed(_cq, contract, _row, _kind, 2, 2).val()
        for _v in _sh.Vertices():
            _r = math.hypot(_v.X, _v.Y)
            _worst_out = max(_worst_out,
                             _r - (float(np.interp(_v.Z, _zt, _rt)) + _eps))
            _worst_in = min(_worst_in,
                            _r - float(np.interp(_v.Z, _zh, _rh)))
    check("ningún álabe atraviesa la CARCASA (recorte contra la vena)",
          _worst_out <= 1e-6,
          f"máx. invasión = {_worst_out:+.4f} mm (ε = {_eps:.2f} mm)")
    check("ningún álabe atraviesa el TAMBOR por dentro",
          _worst_in >= -1e-6,
          f"mín. holgura al cubo = {_worst_in:+.4f} mm")
    _rot0 = contract["stages"][0]["rotor"]
    _sh_r = gg._cq_blade_trimmed(_cq, contract, _rot0, "rotor", 2, 2).val()
    _r_tip_real = max(math.hypot(v.X, v.Y) for v in _sh_r.Vertices())
    # El álabe LLENA la vena: su punta alcanza la línea de punta en su
    # estación (el límite superior ya lo garantiza el check de invasión;
    # esto verifica que el recorte no se pasa y deja el álabe corto).
    _r_case_c = float(np.interp(_rot0["z_center_mm"], _zt, _rt))
    check("la punta del rotor llega a la vena (no queda corta)",
          _r_tip_real >= _r_case_c - 0.25,
          f"r_punta = {_r_tip_real:.2f} mm vs vena {_r_case_c:.2f} mm "
          f"en la misma estación (ε = {_eps:.2f} mm a la carcasa)")
    check("todos los sólidos de álabe son VÁLIDOS para booleanas",
          all(gg._cq_blade_trimmed(_cq, contract, r, k, 2, 2).val().isValid()
              for r, k in _rows[:4]))

    # --- G-04: INTERFERENCIAS DEL ENSAMBLE MONTADO ---------------------
    # Los checks de arriba miden filas SUELTAS contra las polilíneas de
    # annulus. Ninguno de los tres fallos que aparecieron en el STEP era
    # de esa forma: eran pieza contra pieza ya montada (rim de disco
    # dentro del álabe, casco de tambor a través del disco, tirantes
    # atravesando el alma de los discos delanteros). Esto los busca.
    _parts_asm = gg._asm_parts(_asm)
    _smp = gg._interference_sample(_parts_asm)
    _hits = gg.assembly_interferences(contract, parts=_parts_asm)
    check("ninguna pareja de piezas del ensamble comparte material",
          not _hits,
          gg.format_interferences(_hits, len(_smp)).replace("\n", " · "))

    # control NEGATIVO: un chequeo que nunca encuentra nada pasa igual de
    # verde que uno correcto. Se le mete un solape fabricado a mano.
    _d1_shape = [ch for ch in _subs["ROTOR"].children
                 if ch.name == "Disc1"][0].obj.val()
    _fake = [("A/Disc1", _d1_shape),
             ("A/Disc1_desplazado", _d1_shape.translate((0, 0, 1.0)))]
    check("el chequeo DETECTA un solape fabricado a propósito",
          len(gg.assembly_interferences(contract, parts=_fake)) == 1,
          f"{len(gg.assembly_interferences(contract, parts=_fake))} par(es)")
    check("la muestra colapsa las coronas patronadas pero no los discos",
          sum(n.endswith(("blade001", "blade002")) for n, _ in _smp)
          == rec["n_stages"] * 2
          and sum("Disc" in n for n, _ in _smp) == rec["n_stages"]
          and not any(n.endswith("blade003") for n, _ in _smp),
          f"{len(_smp)} de {len(_parts_asm)} piezas")

    # el álabe DEBE bajar por debajo del hub (lleva raíz)
    _b1 = [ch for ch in _subs["ROTOR"].children
           if ch.name == "Rotor1_blade001"][0]
    _r_lo = min(math.hypot(v.X, v.Y) for v in _b1.obj.val().Vertices())
    check("el álabe baja POR DEBAJO de la línea de cubo (tiene raíz)",
          _r_lo < contract["stages"][0]["rotor"]["sections"][0]["r_mm"] - 1.0,
          f"r_mín = {_r_lo:.1f} mm vs hub "
          f"{contract['stages'][0]['rotor']['sections'][0]['r_mm']:.1f} mm")

# ==========================================================================
print("— T17 · contrato versionado phyac-axial-2 (fase 10 · S-01)")

_schema = cx.load_schema()
check("el JSON Schema publicado existe y se declara a sí mismo",
      _schema["properties"]["schema"]["const"] == cx.SCHEMA_ID
      and os.path.getsize(cx.schema_path()) > 4000,
      os.path.relpath(cx.schema_path()))

# --- el contrato real valida ----------------------------------------------
# Se valida el JSON REDONDEADO (dump → load), que es lo que ve un
# consumidor: de paso prueba que todo el contrato es serializable.
_c_ref = json.loads(json.dumps(contract))
_errs = cx.validate(_c_ref)
check("el contrato de referencia valida contra phyac-axial-2",
      not _errs, "; ".join(_errs[:3]) or "0 incumplimientos")

_st_ref = sc.evaluate_structural(THETA_REF, rec)
_c_full = json.loads(json.dumps(
    gg.build_contract(THETA_REF, rec, _st_ref, {"run": "t17"})))
_c_fan2 = json.loads(json.dumps(_c_fan))
_e_full, _e_fan = cx.validate(_c_full), cx.validate(_c_fan2)
check("validan también las variantes: con bloque structural y sin IGV",
      not _e_full and not _e_fan and _c_full["structural"]
      and _c_fan2.get("igv") is None,
      "; ".join((_e_full + _e_fan)[:3]) or "0 incumplimientos")

# --- lo que motivó subir de versión ---------------------------------------
# El contrato ganó estos campos DENTRO de phyac-axial-1, sin avisar. Que
# ahora sean OBLIGATORIOS es justo lo que hace que un consumidor v1 no
# pueda leer un fichero v2 creyendo que lo entiende.
_req_d = set(_schema["properties"]["derived"]["required"])
_req_a = set(_schema["properties"]["assembly"]["required"])
check("los campos de la fase 9 son OBLIGATORIOS en el esquema, no opcionales",
      {"vortex_n", "gas_model", "bleed", "stage_stall_first",
       "ds_machine_J_kgK", "blockage_exit"} <= _req_d
      and {"firtree", "tie_bolt_count", "tie_bolt_d_mm"} <= _req_a)

# --- el validador detecta DE VERDAD ---------------------------------------
# Un validador que acepta todo pasa igual de verde que uno correcto. Se le
# dan seis contratos rotos de formas distintas y tiene que cazar los seis.
def _broken(mut):
    c = json.loads(json.dumps(contract))
    mut(c)
    return len(cx.validate(c))


def _del_firtree(c):
    del c["assembly"]["firtree"]


def _bad_type(c):
    c["derived"]["n_stages"] = "cuatro"


def _out_of_range(c):
    c["derived"]["eta_poly"] = 1.4


def _bad_enum(c):
    c["stages"][0]["rotor"]["sections"][0]["profile"] = "NACA-4412"


def _short_polyline(c):
    c["annulus"]["hub"] = [[0.0, 100.0]]


def _bad_point(c):
    c["annulus"]["tip"][1] = [0.0, 100.0, 7.0]


_caught = [_broken(f) > 0 for f in (_del_firtree, _bad_type, _out_of_range,
                                    _bad_enum, _short_polyline, _bad_point)]
check("el validador caza los 6 modos de rotura (falta, tipo, rango, "
      "enum, longitud, forma)", all(_caught),
      f"{sum(_caught)}/6")
check("el booleano no es booleano: n_stages=True NO es un entero válido",
      _broken(lambda c: c["derived"].update(n_stages=True)) > 0)

# --- versión mayor y frontera con la capa 5c ------------------------------
check("la versión mayor se lee del identificador y es la 2",
      cx.schema_major(contract) == cx.SCHEMA_MAJOR == 2
      and cx.schema_major({"schema": "phyac-axial-7"}) == 7
      and cx.schema_major({"schema": "otra-cosa"}) is None)

_cs_src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "AxialCompressorDesigner", "src",
                            "AxialCompressorDesigner", "PhyACImport.cs"),
               encoding="utf-8").read()
_m_cs = re.search(r"nSCHEMA_MAJOR\s*=\s*(\d+)", _cs_src)
check("la capa 5c (C#) declara la MISMA versión mayor que Python",
      _m_cs is not None and int(_m_cs.group(1)) == cx.SCHEMA_MAJOR,
      f"C# = {_m_cs.group(1) if _m_cs else '?'}, "
      f"Python = {cx.SCHEMA_MAJOR}")
check("el lector C# RECHAZA una versión mayor desconocida en vez de "
      "aplicar defaults",
      "NotSupportedException" in _cs_src
      and "nMajor != nSCHEMA_MAJOR" in _cs_src)

# un contrato de la versión anterior tiene que fallar por el identificador,
# no colarse porque «casi» encaja
_c_v1 = json.loads(json.dumps(contract))
_c_v1["schema"] = "phyac-axial-1"
del _c_v1["assembly"]["firtree"]
check("un contrato v1 se rechaza (identificador y campos que le faltan)",
      cx.schema_major(_c_v1) == 1 and len(cx.validate(_c_v1)) >= 2)

# --- el CLI del validador ---------------------------------------------------
with tempfile.TemporaryDirectory() as _td17:
    _p_ok = os.path.join(_td17, "axial_compressor.json")
    with open(_p_ok, "w", encoding="utf-8") as _f:
        json.dump(_c_ref, _f)
    _p_bad = os.path.join(_td17, "v1.json")
    with open(_p_bad, "w", encoding="utf-8") as _f:
        json.dump(_c_v1, _f)
    check("CLI del validador: 0 en el contrato bueno, 1 en el v1",
          cx._main(["contract_schema.py", _p_ok]) == 0
          and cx._main(["contract_schema.py", _p_bad]) == 1)

# ==========================================================================
print("— T18 · fidelidad L1: SCM propio (fase 11 · F-01/H2)")
import scm_core as sc1     # noqa: E402

check("L1 ya NO depende de nada externo: está siempre disponible",
      pc.l1_available() and pc.l1_unavailable_reason() is None)

_r_l1 = pc.evaluate(THETA_REF, fidelity=pc.Fidelity.L1, use_cache=False)
_r_l0 = pc.evaluate(THETA_REF, fidelity=pc.Fidelity.L0, use_cache=False)
check("L1 resuelve de verdad: el record viene del SCM, no del meanline",
      _r_l1["source"] == "scm_L1", _r_l1["source"])
check("el record L1 conserva el contrato de campos del L0",
      set(_r_l0) <= set(_r_l1)
      and _r_l1["n_stages"] == _r_l0["n_stages"])
check("un diseño INVIABLE no paga el SCM (se queda en L0)",
      pc.evaluate(np.array([8.0, 30_000.0, 0.90, 1.20, 0.90, 0.0, 0.30,
                            2.5, 2.5, 1.0, 288.15, 101_325.0, 60.0]),
                  fidelity=pc.Fidelity.L1,
                  use_cache=False)["source"].startswith("meanline_L0"))

# ==========================================================================
print("— T21 · equilibrio radial y líneas de corriente (fase 11)")

# --- VERIFICACIÓN: el ODE reproduce la forma cerrada de la fase 9.1 -------
# Sin pérdidas ni curvatura, la ecuación de equilibrio radial integrada
# numéricamente TIENE que dar el mismo perfil que
# physics_core.vortex_cx, que es su solución analítica para Cu ∝ r^n. Si
# esto falla, el solver no está resolviendo la ecuación que dice resolver.
_r21 = np.linspace(0.18, 0.28, 9)
_rm21, _cum21, _cxm21 = 0.23, 90.0, 160.0
_A21 = math.pi * (0.28 ** 2 - 0.18 ** 2)
_md21 = (101_325.0 / (287.0 * 288.15) * 0.85) * _cxm21 * _A21
_worst21 = {}
for _n21 in (-1.0, -0.5, 0.0, 0.5):
    _cu21 = np.array([pc.vortex_cu(_cum21, x, _rm21, _n21) for x in _r21])
    _cm21, _, _lim21 = sc1._solve_station_cm(
        _r21, _cu21, np.full(9, pc.h_air(288.15)), np.zeros(9),
        np.full(9, 288.15), np.full(9, _cxm21), np.zeros(9), np.ones(9),
        _md21, 1.0, 288.15, 101_325.0, _cxm21)
    _ref21 = np.array([pc.vortex_cx(_cxm21, _cum21, x, _rm21, _n21)
                       for x in _r21])
    _worst21[_n21] = float(np.max(np.abs(
        _cm21 / np.interp(_rm21, _r21, _cm21)
        - _ref21 / np.interp(_rm21, _r21, _ref21))))
check("el ODE de equilibrio radial reproduce la forma cerrada de vortex_cx",
      max(_worst21.values()) < 5e-3 and _worst21[-1.0] < 1e-12,
      "máx. desvío de forma: "
      + ", ".join(f"n={k:+.1f}: {v:.1e}" for k, v in _worst21.items()))

_scm = _r_l1["scm"]
check("el SCM resuelve sobre N líneas de corriente y 2n+1 estaciones",
      _scm["n_streamlines"] == sc1.N_STREAMLINES >= 5
      and _scm["n_stations"] == 2 * _r_l1["n_stages"] + 1,
      f"{_scm['n_streamlines']} líneas × {_scm['n_stations']} estaciones "
      f"en {_scm['iterations']} iteraciones")

# --- continuidad POR ESTACIÓN --------------------------------------------
_sl = np.array(_scm["streamlines_mm"]) / 1000.0
check("las líneas de corriente son monótonas y viven dentro de la vena",
      all(np.all(np.diff(_sl[k]) > 0) for k in range(len(_sl)))
      and all(abs(_sl[k][0] - _scm["r_hub_mm"][k] / 1000.0) < 1e-9
              and abs(_sl[k][-1] - _scm["r_tip_mm"][k] / 1000.0) < 1e-9
              for k in range(len(_sl))))

# --- el vórtice libre DEBE reproducir el meanline -------------------------
# Es la condición en la que las hipótesis del meanline son exactas (Cx
# radialmente constante, trabajo uniforme): si ahí las dos capas no
# coinciden, la diferencia no es fidelidad, es un error.
check("con vórtice libre L1 reproduce el meanline (±1% PR, ±1 pt η)",
      abs(_r_l1["PR"] / _r_l0["PR"] - 1.0) < 0.01
      and abs(_r_l1["eta_poly"] - _r_l0["eta_poly"]) < 0.01,
      f"PR {_r_l0['PR']:.4f} → {_r_l1['PR']:.4f} "
      f"({100 * (_r_l1['PR'] / _r_l0['PR'] - 1):+.2f}%), "
      f"η {_r_l0['eta_poly']:.4f} → {_r_l1['eta_poly']:.4f}")

# --- y con torbellino controlado NO debe reproducirlo ---------------------
# Ahí la forma cerrada del meanline es una aproximación y el trabajo deja
# de ser uniforme en el span: es donde L1 aporta el residual que la capa 2
# necesita aprender. Un L1 que devuelve siempre lo mismo que L0 no es un
# peldaño de fidelidad.
_vn21 = pc.VORTEX_N
pc.VORTEX_N = -0.5
_r0_cv = pc.evaluate(THETA_REF, fidelity=pc.Fidelity.L0, use_cache=False)
_r1_cv = pc.evaluate(THETA_REF, fidelity=pc.Fidelity.L1, use_cache=False)
pc.VORTEX_N = _vn21
check("con torbellino controlado (n=−0.5) L1 SÍ se separa de L0",
      _r1_cv["source"] == "scm_L1"
      and abs(_r1_cv["PR"] / _r0_cv["PR"] - 1.0) > 0.01,
      f"PR {_r0_cv['PR']:.4f} → {_r1_cv['PR']:.4f} "
      f"({100 * (_r1_cv['PR'] / _r0_cv['PR'] - 1):+.2f}%)")

# --- independencia del número de líneas ----------------------------------
_pr_sl = {}
for _n_sl in (7, 11):
    _pr_sl[_n_sl] = sc1.solve(THETA_REF, _r_l0, n_streamlines=_n_sl)["PR"]
check("el resultado no depende del nº de líneas (7 vs 11, <0.5%)",
      abs(_pr_sl[7] / _pr_sl[11] - 1.0) < 0.005,
      f"PR 7 líneas {_pr_sl[7]:.4f} vs 11 líneas {_pr_sl[11]:.4f}")
check("menos de 5 líneas se RECHAZA (la derivada radial no tendría sentido)",
      not _ok_call(lambda: sc1.solve(THETA_REF, _r_l0, n_streamlines=3)))

# --- el trabajo deja de ser uniforme en el span ---------------------------
_dh21 = np.array(_scm["dh0_span"])
check("L1 entrega el reparto RADIAL de trabajo que el meanline no tiene",
      len(_dh21) == _scm["n_streamlines"] and _scm["work_spread"] > 0.005
      and np.all(_dh21 > 0),
      f"Δh₀ de {_dh21.min() / 1000:.1f} a {_dh21.max() / 1000:.1f} kJ/kg "
      f"({100 * _scm['work_spread']:.1f}% del trabajo medio)")
check("y el perfil de Cm tiene déficit en LAS DOS paredes",
      (np.array(_scm["cm_exit"])[0] < max(_scm["cm_exit"])
       and np.array(_scm["cm_exit"])[-1] < max(_scm["cm_exit"])),
      "Cm salida: " + " ".join(f"{v:.0f}" for v in _scm["cm_exit"]))

# --- la divergencia se REPORTA, no se traga -------------------------------
check("SCMDiverged es un error propio, no un None silencioso",
      issubclass(sc1.SCMDiverged, Exception)
      and not _ok_call(lambda: sc1.solve(
          THETA_REF, _r_l0, n_streamlines=4)))

# --- la guarda de DERIVA vs L0 --------------------------------------------
# El annulus lo dimensiona L0 con su Cx uniforme; sobre muchas etapas la
# diferencia con el perfil de L1 se compone a través del álabe de ángulo
# fijo (medido hasta ±27% de PR en máquinas de 7-8 etapas, ver
# validation/BENCH_SCM.md). Fuera de la ventana el punto se rechaza.
_r_fake = json.loads(json.dumps({k: v for k, v in _r_l0.items()
                                 if k not in ("frozen", "stage_table")}))
_r_fake["frozen"] = _r_l0["frozen"]
_r_fake["stage_table"] = _r_l0["stage_table"]
_r_fake["PR"] = _r_l0["PR"] * 2.0        # L0 "dice" el doble
check("una DERIVA grande respecto a L0 se rechaza, no se devuelve",
      not _ok_call(lambda: sc1.solve(THETA_REF, _r_fake)),
      f"ventana ±{100 * sc1.PR_WINDOW:.0f}%")
check("y dentro de la ventana el mismo diseño SÍ resuelve",
      _ok_call(lambda: sc1.solve(THETA_REF, _r_l0)))

# ==========================================================================
print("— T20 · paridad capa 5c (C#/PicoGK) ↔ via CadQuery (fase 10 · G-01)")
# El STL es la ruta de fabricacion y el STEP la de re-CAD. Entre la fase
# 8.2 y la 10 se separaron sin que nada avisara. Esto compara la parte
# BARATA (matematica pura: perfil de abeto, longitud de raiz, pendiente
# de plataforma, tornilleria) en cada corrida de la suite; la paridad
# cara por volumenes vive en validation/parity_stl_step.py.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "validation"))
import parity_stl_step as ps        # noqa: E402

_cs_ok, _cs_why = ps.csharp_available()
if not _cs_ok:
    check("paridad 5c: capa C# no compilada "
          + ("(EXIGIDO por PHYAC_REQUIRE_STL)" if REQUIRE_STL else "(SKIP)"),
          not REQUIRE_STL, _cs_why)
else:
    with tempfile.TemporaryDirectory() as _td20:
        _p20 = os.path.join(_td20, "axial_compressor.json")
        with open(_p20, "w", encoding="utf-8") as _f:
            json.dump(_c_ref, _f)
        _cs20 = ps.csharp_profile(_p20)
        _py20 = ps.python_profile(_c_ref)
        _e20 = ps.compare_profiles(_cs20, _py20)
        check("las dos vias construyen la MISMA raiz de abeto "
              "(perfil punto a punto)", not _e20,
              "; ".join(_e20[:3]) or f"{len(_cs20['root_profile'])} puntos, "
              f"{len(_cs20['rows'])} filas de rotor")
        check("la capa 5c lee tirantes y rebaje de rim del contrato",
              _cs20["tie_bolt_count"]
              == _c_ref["assembly"]["tie_bolt_count"]
              and abs(_cs20["tie_bolt_d_mm"]
                      - _c_ref["assembly"]["tie_bolt_d_mm"]) < 1e-6
              and abs(_cs20["rim_relief_mm"]
                      - _c_ref["assembly"]["rim_relief_mm"]) < 1e-6,
              f"{_cs20['tie_bolt_count']} tirantes de "
              f"{_cs20['tie_bolt_d_mm']:g} mm")
        check("el numero de secciones del loft de abeto coincide "
              "(la ranura no puede encajar si no)",
              _cs20["n_sections"] == gg.N_FIRTREE_SECTIONS
              and abs(_cs20["slot_axial_factor"] - 1.05) < 1e-9,
              f"{_cs20['n_sections']} secciones")
        # el comparador tiene que DETECTAR una separacion real
        _py_bad = json.loads(json.dumps(_py20))
        _py_bad["root_profile"][3][0] += 0.05
        _py_bad["rows"][0]["root_axial_mm"] += 1.0
        check("el comparador de paridad caza una separacion fabricada",
              len(ps.compare_profiles(_cs20, _py_bad)) >= 2)

# ==========================================================================
print("— T22 · fuera de diseño calificado contra medida (fase 12 · F-02)")
# Hasta la fase 12 la campaña de validación solo calificaba el PUNTO DE
# DISEÑO. El margen de bombeo, que desde la fase 9 es la restricción dura
# que más recorta el espacio, no estaba contrastado con NINGÚN dato
# medido. Esto comprueba que la maquinaria de calificación del mapa
# existe, corre y no se degrada en silencio.
from validation.machines import OFFDESIGN as _OD      # noqa: E402
from validation.validate import run_offdesign as _ro  # noqa: E402

check("hay al menos un caso de FUERA DE DISEÑO con fuente citada",
      len(_OD) >= 1
      and all(o.get("source") and o.get("measured") for o in _OD),
      f"{len(_OD)} caso(s)")
_r22 = [_ro(o) for o in _OD]
_keys22 = {it["key"] for r in _r22 for it in r["items"]}
check("se califican los DOS extremos del rango de gasto, no solo uno",
      {"mdot_choke", "stall_over_choke"} <= _keys22,
      ", ".join(sorted(_keys22)))
check("cada cantidad declara tolerancia OBJETIVO y guarda interina",
      all(0 < it["tol"] <= it["guard"] for r in _r22 for it in r["items"]))
# Regresión: los deltas medidos hoy. No se comprueba PASS (el choke FALLA
# el objetivo del 5% y eso está documentado), sino que no se muevan a
# peor en silencio.
_d22 = {it["key"]: it["d_rel"] for r in _r22 for it in r["items"]}
check("el criterio de BOMBEO sitúa el stall donde lo midió NASA "
      "(ṁ_stall/ṁ_choke)",
      abs(_d22["stall_over_choke"]) < 0.03,
      f"modelo/medido = {1 + _d22['stall_over_choke']:.3f} "
      f"({100 * _d22['stall_over_choke']:+.2f}% vs 0.925 medido)")
check("el criterio de CHOKE sigue siendo optimista y no ha empeorado",
      0.0 < _d22["mdot_choke"] < 0.12,
      f"{100 * _d22['mdot_choke']:+.2f}% (objetivo ±5%, guarda ±12%) — "
      "brecha ABIERTA, ver docs/VALIDATION.md")
check("todas las cantidades dentro de la guarda con la que corre --strict",
      all(it["ok_guard"] for r in _r22 for it in r["items"]))

# --- la SPEEDLINE completa, punto a punto ---------------------------------
from validation.machines import SPEEDLINES as _SL     # noqa: E402
from validation.validate import (run_speedline as _rs,  # noqa: E402
                                 _speedline_items as _si)

check("hay una speedline medida con sus puntos y su fuente",
      len(_SL) >= 1
      and all(len(s["points"]) >= 5 and s.get("source") for s in
              _SL.values()),
      ", ".join(f"{k}: {len(v['points'])} pts" for k, v in _SL.items()))
_sp22 = [_rs(k) for k in _SL]
_it22 = {it["key"]: it for sp in _sp22 for it in _si(sp)}
check("se califican las cinco métricas de la característica",
      {"pr_max_abs", "eta_max_abs", "pr_peak", "eta_peak", "slope_rel"}
      <= set(_it22))
# Regresión sobre lo medido hoy: no se exige PASS (η punto a punto FALLA
# el objetivo y está documentado), se exige que no empeore en silencio.
check("el modelo sigue la característica en PR a lo largo de la línea",
      _it22["pr_max_abs"]["val"] < 0.05,
      f"máx |ΔPR| = {100 * _it22['pr_max_abs']['val']:.2f}% sobre "
      f"{len(_sp22[0]['rows'])} puntos")
check("y acierta la PENDIENTE, que es la que mueve el punto de operación",
      _it22["slope_rel"]["val"] < 0.20,
      f"dPR/dṁ modelo {_sp22[0]['slope_model']:+.4f} vs medida "
      f"{_sp22[0]['slope_meas']:+.4f} /(kg/s) "
      f"({100 * _sp22[0]['slope_rel']:+.1f}%)")
check("el η punto a punto sigue fuera de objetivo y no ha empeorado",
      0.02 < _it22["eta_max_abs"]["val"] < 0.06,
      f"máx |Δη| = {100 * _it22['eta_max_abs']['val']:.2f} pts en el "
      "extremo de choke (objetivo 3 pts) — brecha ABIERTA, "
      "ver docs/VALIDATION.md")
check("toda la speedline dentro de la guarda de --strict",
      all(it["ok_guard"] for sp in _sp22 for it in _si(sp)))

# ==========================================================================
print(f"\n{_n_pass}/{_n_pass + _n_fail} checks OK"
      + (f" — {_n_fail} FALLOS" if _n_fail else ""))
sys.exit(1 if _n_fail else 0)
