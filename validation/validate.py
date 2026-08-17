#!/usr/bin/env python3
"""
QUASAR Phy-AC · validation/validate.py
======================================
Corre la campaña de VALIDACIÓN (¿resolvemos las ecuaciones CORRECTAS?)
contra las máquinas de `machines.py` y regenera `validation/RESULTS.md`.

    python validation/validate.py                  # tabla + RESULTS.md
    python validation/validate.py --strict         # exit 1 si PR falla, un
                                                   # ancla se mueve o η
                                                   # excede la guarda
    python validation/validate.py --freeze-anchors # congela las anclas de
                                                   # regresión al valor
                                                   # actual del meanline

Distinción VV&UQ del proyecto: `test_phyac.py` VERIFICA (ecuaciones bien
resueltas); este script VALIDA (predicciones vs máquinas medidas).

Mapeo espec → θ (ver machines.py): r_tip desde U_tip/ω; φ1 invertido por
bisección para que la continuidad reproduzca ese r_tip; ψ_mid desde el
trabajo MEDIDO (ΔT0 de PR y η publicados). Así el modelo recibe el trabajo
real y se califica su predicción de PÉRDIDAS → (η, PR).
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
from datetime import date

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import numpy as np

import physics_core as pc
from physics_core import CP, GAMMA, Fidelity, _meanline, evaluate
from validation.machines import (MACHINES, OFFDESIGN, SPEEDLINES,
                                 REGRESSION_ANCHORS, PR_TOL_REL,
                                 ETA_TOL_PTS, ETA_GUARD_PTS)

KIND_LABEL = {
    "rotor": "rotor aislado (t-t)",
    "stage": "etapa completa (t-t)",
    "machine": "máquina multietapa (t-t)",
}


# ---------------------------------------------------------------------------
# Espec publicada → θ 13-D
# ---------------------------------------------------------------------------
def measured_dh0(spec: dict, measured: dict) -> float:
    """Δh0 total [J/kg] desde (PR, η) medidos — el trabajo real.

    Desde la fase 9 el meanline usa gas caloríficamente imperfecto, así
    que el trabajo que se le INYECTA tiene que derivarse con el mismo gas
    (si no, se le pide a la máquina una ψ inconsistente con su propio cp
    y el PR sale sesgado). Con la función de entropía phi(T):

        η_isen dada:  phi(T2s) = phi(T1) + R·ln(PR)
                      Δh0 = [h(T2s) − h(T1)] / η_isen
        η_poly dada:  phi(T2)  = phi(T1) + R·ln(PR)/η_poly
                      Δh0 = h(T2) − h(T1)
    """
    PR, T0 = measured["PR"], spec["T0"]
    h1, p1 = pc.h_air(T0), pc.phi_air(T0)
    if "eta_isen" in measured:
        T2s = pc.T_from_phi(p1 + pc.RGAS * math.log(PR), T0 * 1.2)
        return (pc.h_air(T2s) - h1) / measured["eta_isen"]
    T2 = pc.T_from_phi(p1 + pc.RGAS * math.log(PR) / measured["eta_poly"],
                       T0 * 1.5)
    return pc.h_air(T2) - h1


def measured_dT0(spec: dict, measured: dict) -> float:
    """ΔT0 equivalente al trabajo medido (solo para reporte)."""
    T0 = spec["T0"]
    return pc.T_from_h(pc.h_air(T0) + measured_dh0(spec, measured), T0) - T0


def build_theta(spec: dict, measured: dict,
                slopes: dict | None = None) -> np.ndarray:
    """Construye θ que reproduce annulus (r_tip vía φ1) y trabajo (ψ_mid).

    `slopes` (fase 8, multietapa): {"phi_slope", "Rx_slope"} describe la
    distribución POR ETAPA de la máquina (θ 15-D); la bisección de φ1
    sigue anclando la etapa frontal al annulus publicado. Sin slopes el θ
    es el legacy de 13 (pad_theta pone pendientes 0)."""
    omega = spec["RPM"] * 2 * math.pi / 60.0
    r_tip_target = spec["U_tip"] / omega
    r_mean = 0.5 * (1.0 + spec["HTR"]) * r_tip_target
    Um = omega * r_mean
    psi = measured_dh0(spec, measured) / (spec["n_stages"] * Um ** 2)
    tail = ([slopes.get("phi_slope", 0.0), slopes.get("Rx_slope", 0.0)]
            if slopes else [])

    def theta_of(phi1: float) -> np.ndarray:
        return np.array([spec["n_stages"], spec["RPM"], spec["HTR"], phi1,
                         psi, 0.0, spec["Rx_est"], spec["sigma_r"],
                         spec["sigma_s"], spec["AR"], spec["T0"],
                         spec["P0"], spec["mdot"]] + tail)

    lo, hi = 0.25, 1.10
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        r_tip = _meanline(theta_of(mid))["r_tip_in_mm"] / 1000.0
        # φ1 ↑ ⇒ Cx ↑ ⇒ annulus más pequeño ⇒ r_tip ↓
        if r_tip > r_tip_target:
            lo = mid
        else:
            hi = mid
    return theta_of(0.5 * (lo + hi))


# ---------------------------------------------------------------------------
# Métricas del modelo según el plano/kind de calificación
# ---------------------------------------------------------------------------
def model_metrics(rec: dict, spec: dict, kind: str) -> dict:
    """PR y η (isen/poly) del modelo en el plano medido."""
    T0 = spec["T0"]
    if kind == "rotor":
        # PR/η del ROTOR de la etapa 1, derivados del desglose de pérdidas.
        # El endwall de la etapa se reparte por cabeza dinámica (el rotor
        # se lleva su parte); la fuga de punta es toda suya.
        s = rec["stage_table"][0]
        L = s["losses"]
        dh0 = s["dh0"]
        ds_ew = max(s["ds_stage_J_kgK"] - L["ds_rotor"] - L["ds_stator"]
                    - L["ds_clearance"], 0.0)
        f_rot = L["ds_rotor"] / max(L["ds_rotor"] + L["ds_stator"], 1e-9)
        dh_rot = s["T0_out_K"] * (L["ds_rotor"] + L["ds_clearance"]
                                  + f_rot * ds_ew)
        eta_rot = float(np.clip((dh0 - dh_rot) / max(dh0, 1e-3), 0.05, 0.99))
        T0_out = pc.T_from_h(pc.h_air(T0) + dh0, T0 + dh0 / 1005.0)
        PR = pc.pr_from_work(T0, eta_rot * dh0)
    else:
        PR = rec["PR"]
        T0_out = rec["T0_out"]
    # η exactas para gas imperfecto (mismas definiciones que physics_core)
    if PR > 1.0 + 1e-9 and T0_out > T0 + 1e-9:
        dphi = pc.phi_air(T0_out) - pc.phi_air(T0)
        eta_poly = pc.RGAS * math.log(PR) / max(dphi, 1e-9)
        T_isen = pc.T_from_phi(pc.phi_air(T0) + pc.RGAS * math.log(PR), T0_out)
        eta_isen = (pc.h_air(T_isen) - pc.h_air(T0)) / \
            max(pc.h_air(T0_out) - pc.h_air(T0), 1e-6)
    else:
        eta_poly = eta_isen = 0.0
    return dict(PR=PR, eta_isen=eta_isen, eta_poly=eta_poly)


def run_machine(m: dict) -> dict:
    # Inyecta la holgura de punta PUBLICADA de la máquina (ε absoluta en
    # mm — parámetro de módulo del meanline) y restaura siempre.
    import physics_core as _pc
    eps_saved = _pc.TIP_CLEARANCE_MM
    if "eps_tip_mm" in m:
        _pc.TIP_CLEARANCE_MM = float(m["eps_tip_mm"])
    try:
        theta = build_theta(m["spec"], m["measured"], m.get("slopes"))
        rec = evaluate(theta, fidelity=Fidelity.L0, use_cache=False,
                       calibrate=False)
    finally:
        _pc.TIP_CLEARANCE_MM = eps_saved
    model = model_metrics(rec, m["spec"], m["kind"])

    meas = m["measured"]
    eta_key = "eta_isen" if "eta_isen" in meas else "eta_poly"
    dPR_rel = (model["PR"] - meas["PR"]) / meas["PR"]
    deta = model[eta_key] - meas[eta_key]
    tol = m.get("tol", {})
    pr_tol = tol.get("PR_rel", PR_TOL_REL)
    eta_tol = tol.get("eta_pts", ETA_TOL_PTS)

    return dict(
        name=m["name"], kind=m["kind"], eta_key=eta_key, theta=theta,
        PR_model=model["PR"], PR_meas=meas["PR"], dPR_rel=dPR_rel,
        eta_model=model[eta_key], eta_meas=meas[eta_key], deta=deta,
        PR_pass=abs(dPR_rel) < pr_tol,
        eta_pass=abs(deta) < eta_tol,
        eta_guard=abs(deta) < ETA_GUARD_PTS,
        record=rec, machine=m,
    )


def run_offdesign(o: dict) -> dict:
    """Califica una cantidad de FUERA DE DISEÑO contra medida (F-02).

    Reusa la entrada de MACHINES a la que se engancha: mismo θ invertido,
    misma holgura publicada. Lo que cambia es QUÉ se compara — aquí el
    límite de gasto del mapa, no el punto de diseño.
    """
    m = next((x for x in MACHINES if x["name"] == o["machine"]), None)
    if m is None:
        raise KeyError(f"OFFDESIGN apunta a una máquina que no existe: "
                       f"{o['machine']}")
    eps_saved = pc.TIP_CLEARANCE_MM
    if "eps_tip_mm" in m:
        pc.TIP_CLEARANCE_MM = float(m["eps_tip_mm"])
    try:
        theta = build_theta(m["spec"], m["measured"], m.get("slopes"))
        rec = evaluate(theta, fidelity=Fidelity.L0, use_cache=False,
                       calibrate=False)
        rpm = float(theta[1]) * float(o.get("speed_frac", 1.0))
        mdot_d = float(theta[12])
        mdot_choke = pc._choke_mdot(theta, rpm, rec["frozen"], mdot_d)
        # Bombeo por el MISMO criterio dual que usa el mapa (capacidad de
        # Koch agotada O pendiente nula de la speedline), sobre la misma
        # speedline. Comparar el ancho del mapa exige que los dos extremos
        # salgan del mismo barrido.
        mdot_surge, _srg = pc._surge_mdot(theta, rpm, rec["frozen"], mdot_d,
                                          mdot_choke)
        model = dict(mdot_choke=mdot_choke,
                     mdot_stall=mdot_surge,
                     stall_over_choke=mdot_surge / max(mdot_choke, 1e-9))
    finally:
        pc.TIP_CLEARANCE_MM = eps_saved

    items = []
    for key, meas in o["measured"].items():
        mod = model[key]
        d = (mod - meas) / meas
        tol = o.get("tol", {}).get(f"{key}_rel", 0.05)
        guard = o.get("guard", {}).get(f"{key}_rel", tol)
        items.append(dict(key=key, model=mod, meas=meas, d_rel=d,
                          tol=tol, guard=guard,
                          ok=abs(d) < tol, ok_guard=abs(d) < guard))
    return dict(name=f"{o['machine']} · {100 * o.get('speed_frac', 1.0):.0f}% N",
                machine=o["machine"], items=items, offdesign=o,
                mdot_design=mdot_d, record=rec)


def run_speedline(key: str) -> dict:
    """Califica la SPEEDLINE completa punto a punto contra medida (F-02).

    Los extremos del rango de gasto los cubre `run_offdesign`. Esto
    responde las otras dos preguntas de F-02: ¿acierta el modelo el PR y
    el η a lo largo de la característica, y acierta su PENDIENTE? La
    pendiente es la que gobierna cómo se mueve el punto de operación
    cuando cambia la contrapresión, así que es la que decide si el margen
    de bombeo significa algo.

    Se evalúa con la GEOMETRÍA CONGELADA del diseño (igual que el mapa) y
    se extrae la misma métrica que el punto de diseño de esa máquina
    (para un caso `rotor`, el PR/η del rotor desde el desglose de
    pérdidas), porque los valores medidos son del rotor aislado.
    """
    sl = SPEEDLINES[key]
    m = next((x for x in MACHINES if x["name"] == sl["machine"]), None)
    if m is None:
        raise KeyError(f"SPEEDLINES['{key}'] apunta a una máquina que no "
                       f"existe: {sl['machine']}")
    eps_saved = pc.TIP_CLEARANCE_MM
    if "eps_tip_mm" in m:
        pc.TIP_CLEARANCE_MM = float(m["eps_tip_mm"])
    try:
        theta = build_theta(m["spec"], m["measured"], m.get("slopes"))
        rec_d = evaluate(theta, fidelity=Fidelity.L0, use_cache=False,
                         calibrate=False)
        rpm = float(theta[1]) * float(sl.get("speed_frac", 1.0))
        rows = []
        for lbm, kg, pr_m, tr_m, eta_m, tag in sl["points"]:
            rec = pc.offdesign(theta, rpm, kg, frozen=rec_d["frozen"])
            mod = model_metrics(rec, m["spec"], m["kind"])
            rows.append(dict(mdot=kg, tag=tag,
                             PR_meas=pr_m, PR_model=mod["PR"],
                             dPR_rel=mod["PR"] / pr_m - 1.0,
                             eta_meas=eta_m, eta_model=mod["eta_isen"],
                             deta=mod["eta_isen"] - eta_m,
                             min_SM=rec["min_SM"]))
    finally:
        pc.TIP_CLEARANCE_MM = eps_saved

    mdots = np.array([r["mdot"] for r in rows])
    # Pendiente dPR/dṁ por mínimos cuadrados sobre TODO el rango medido:
    # es la forma de la característica, no la diferencia entre extremos.
    sl_meas = float(np.polyfit(mdots, [r["PR_meas"] for r in rows], 1)[0])
    sl_mod = float(np.polyfit(mdots, [r["PR_model"] for r in rows], 1)[0])
    return dict(
        name=f"{sl['machine']} · {100 * sl.get('speed_frac', 1.0):.0f}% N "
             f"· {len(rows)} puntos",
        key=key, speedline=sl, rows=rows,
        pr_max_abs=max(abs(r["dPR_rel"]) for r in rows),
        pr_mean=sum(r["dPR_rel"] for r in rows) / len(rows),
        eta_max_abs=max(abs(r["deta"]) for r in rows),
        pr_peak_meas=max(r["PR_meas"] for r in rows),
        pr_peak_model=max(r["PR_model"] for r in rows),
        eta_peak_meas=max(r["eta_meas"] for r in rows),
        eta_peak_model=max(r["eta_model"] for r in rows),
        slope_meas=sl_meas, slope_model=sl_mod,
        slope_rel=sl_mod / sl_meas - 1.0,
        tol=sl.get("tol", {}), guard=sl.get("guard", {}))


def run_anchor(a: dict) -> dict:
    theta = np.asarray(a["theta"], dtype=float)
    rec = evaluate(theta, fidelity=Fidelity.L0, use_cache=False,
                   calibrate=False)
    if not a["expect"]:
        return dict(name=a["name"], ok=True, diffs={}, frozen=False,
                    feasible=rec["feasible"], machine=a)
    diffs = {k: (rec[k], v, abs(rec[k] - v) / max(abs(v), 1e-9))
             for k, v in a["expect"].items()}
    ok = (all(d[2] < a["rtol"] for d in diffs.values())
          and rec["feasible"] == a["feasible"])
    return dict(name=a["name"], ok=ok, diffs=diffs, frozen=True,
                feasible=rec["feasible"], machine=a)


ANCHOR_KEYS = ("PR", "eta_poly", "eta_isen", "T0_out", "U_tip")


def freeze_anchors() -> None:
    """Reescribe machines.py con los `expect` actuales del meanline."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "machines.py")
    src = open(path, encoding="utf-8").read()
    for a in REGRESSION_ANCHORS:
        rec = evaluate(np.asarray(a["theta"], dtype=float),
                       fidelity=Fidelity.L0, use_cache=False,
                       calibrate=False)
        expect = ", ".join(f"{k}={rec[k]:.6f}" for k in ANCHOR_KEYS)
        src = re.sub(r"expect=dict\([^)]*\)",
                     f"expect=dict({expect})", src, count=1)
        src = re.sub(r"feasible=(True|False)",
                     f"feasible={rec['feasible']}", src, count=1)
    open(path, "w", encoding="utf-8").write(src)
    print(f"→ anclas congeladas en {path}")


def _fmt_status(ok: bool) -> str:
    return "PASS" if ok else "FAIL"


def _speedline_items(sp: dict) -> list[dict]:
    """Las cinco métricas de una speedline, con su tolerancia y su guarda."""
    tol, guard = sp["tol"], sp["guard"]
    out = []
    for key, val, txt, d_txt in (
        ("pr_max_abs", sp["pr_max_abs"],
         f"máx |ΔPR| en la línea", f"{100 * sp['pr_max_abs']:.2f}%"),
        ("eta_max_abs", sp["eta_max_abs"],
         f"máx |Δη| en la línea", f"{100 * sp['eta_max_abs']:.2f}p"),
        ("pr_peak", abs(sp["pr_peak_model"] / sp["pr_peak_meas"] - 1.0),
         f"{sp['pr_peak_model']:.3f}/{sp['pr_peak_meas']:.3f}",
         f"{100 * (sp['pr_peak_model'] / sp['pr_peak_meas'] - 1):+.2f}%"),
        ("eta_peak", abs(sp["eta_peak_model"] - sp["eta_peak_meas"]),
         f"{sp['eta_peak_model']:.3f}/{sp['eta_peak_meas']:.3f}",
         f"{100 * (sp['eta_peak_model'] - sp['eta_peak_meas']):+.2f}p"),
        ("slope_rel", abs(sp["slope_rel"]),
         f"{sp['slope_model']:+.4f}/{sp['slope_meas']:+.4f}",
         f"{100 * sp['slope_rel']:+.1f}%"),
    ):
        t = tol.get(key, 0.05)
        g = guard.get(key, t)
        out.append(dict(key=key, val=val, txt=txt, d_txt=d_txt, tol=t,
                        guard=g, ok=val < t, ok_guard=val < g))
    # etiquetas legibles para las dos primeras
    out[0]["txt"] = f"{100 * sp['pr_mean']:+.2f}% media"
    out[1]["txt"] = "13 puntos"
    return out


def write_results_md(results: list[dict], anchors: list[dict],
                     path: str, offd: list[dict] | None = None,
                     speed: list[dict] | None = None) -> None:
    lines = [
        "# Resultados de validación — Quasar Phy-AC",
        "",
        f"Generado por `validation/validate.py` el {date.today()} "
        "(meanline L0, sin calibración afín). **Regenerar tras cualquier "
        "cambio en `physics_core.py`** — este archivo se versiona como "
        "evidencia.",
        "",
        "Metodología: el θ de cada máquina reproduce su annulus (r_tip vía "
        "φ1) y su TRABAJO medido (ψ desde ΔT0 publicado); se califica la "
        "predicción de pérdidas → (η, PR). Tolerancias por máquina en "
        "`machines.py` (monoetapa 5%/2 pts; transónicos y multietapa "
        "relajados — ver notas).",
        "",
        "## Máquinas medidas",
        "",
        "| Máquina | Plano | PR modelo | PR medido | ΔPR | η modelo | "
        "η medida | Δη [pts] | PR | η |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        lines.append(
            f"| {r['name']} | {KIND_LABEL[r['kind']]} | "
            f"{r['PR_model']:.3f} | {r['PR_meas']:.3f} | "
            f"{r['dPR_rel']:+.2%} | {r['eta_model']:.3f} "
            f"({r['eta_key']}) | {r['eta_meas']:.3f} | "
            f"{100 * r['deta']:+.1f} | {_fmt_status(r['PR_pass'])} | "
            f"{_fmt_status(r['eta_pass'])} |")
    # --- fuera de diseño (F-02) -------------------------------------------
    if offd:
        lines += [
            "", "## Fuera de diseño (F-02)", "",
            "Primera calificación del MAPA contra medida. Hasta la fase 12 "
            "la campaña solo calificaba el punto de diseño, y sin embargo "
            "el margen de bombeo es desde la fase 9 la restricción DURA "
            "que más recorta el espacio de diseño.", "",
            "| Caso | Cantidad | Modelo | Medido | Δ | Objetivo | Guarda | |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for o in offd:
            for it in o["items"]:
                lines.append(
                    f"| {o['name']} | `{it['key']}` | {it['model']:.3f} | "
                    f"{it['meas']:.3f} | {it['d_rel']:+.2%} | "
                    f"±{it['tol']:.0%} | ±{it['guard']:.0%} | "
                    f"{_fmt_status(it['ok'])} |")
        for o in offd:
            od = o["offdesign"]
            lines += ["", f"**{o['name']}** — fuente: {od['source']}", ""]
            if od.get("notes"):
                lines += [od["notes"], ""]

    # --- speedline completa (F-02) ----------------------------------------
    if speed:
        for sp in speed:
            lines += [
                "", f"### Speedline medida — {sp['name']}", "",
                f"Fuente: {sp['speedline']['source']}", "",
                "| Métrica | Modelo | Medido | Δ | Objetivo | |",
                "|---|---|---|---|---|---|",
            ]
            for it in _speedline_items(sp):
                lines.append(
                    f"| `{it['key']}` | {it['txt']} | | {it['d_txt']} | "
                    f"{it['tol']:.2f} | {_fmt_status(it['ok'])} |")
            lines += ["", "Punto a punto:", "",
                      "| ṁ [kg/s] | PR med | PR mod | ΔPR | η med | η mod "
                      "| Δη [pts] | SM modelo | |",
                      "|---|---|---|---|---|---|---|---|---|"]
            for r in sp["rows"]:
                lines.append(
                    f"| {r['mdot']:.3f} | {r['PR_meas']:.3f} | "
                    f"{r['PR_model']:.3f} | {r['dPR_rel']:+.2%} | "
                    f"{r['eta_meas']:.3f} | {r['eta_model']:.3f} | "
                    f"{100 * r['deta']:+.2f} | {r['min_SM']:+.3f} | "
                    f"{r['tag']} |")

    lines += ["", "### Detalle por máquina", ""]
    for r in results:
        m = r["machine"]
        t = r["theta"]
        lines += [
            f"**{r['name']}** — fuente: {m['source']}", "",
            f"θ construido: n={t[0]:.0f}, RPM={t[1]:.0f}, HTR={t[2]:.3f}, "
            f"φ1={t[3]:.3f}, ψ={t[4]:.3f}, Rx={t[6]:.2f}, "
            f"σr={t[7]:.2f}, σs={t[8]:.2f}, AR={t[9]:.2f}", "",
            f"Notas: {m['notes']}", ""]
    lines += [
        "## Anclas de regresión internas (no son mediciones)",
        "",
        "| Ancla | Estado | Detalle |",
        "|---|---|---|",
    ]
    for a in anchors:
        if not a["frozen"]:
            det = "sin congelar (correr --freeze-anchors)"
        else:
            det = "; ".join(f"{k}={got:.4f} (esp. {exp:.4f})"
                            for k, (got, exp, _) in a["diffs"].items())
        lines.append(f"| {a['name']} | {_fmt_status(a['ok'])} | {det} |")
    lines += [
        "",
        "## Cómo añadir una máquina",
        "",
        "1. Reunir espec publicada (U_tip, HTR, mdot, RPM, PR, η con tipo "
        "declarado) — sin datos reconstruidos de memoria.",
        "2. Añadir la entrada a `validation/machines.py`.",
        "3. `python validation/validate.py` y versionar este archivo.",
        "",
    ]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def main(argv=None) -> int:
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    ap = argparse.ArgumentParser(description="Campaña de validación Phy-AC")
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    ap.add_argument("--freeze-anchors", action="store_true")
    args = ap.parse_args(argv)

    if args.freeze_anchors:
        freeze_anchors()
        return 0

    results = [run_machine(m) for m in MACHINES]
    offd = [run_offdesign(o) for o in OFFDESIGN]
    speed = [run_speedline(k) for k in SPEEDLINES]
    anchors = [run_anchor(a) for a in REGRESSION_ANCHORS]

    print(f"{'Máquina':26s} {'plano':8s} {'PR mod/med':>14s} {'ΔPR':>8s} "
          f"{'η mod/med':>14s} {'Δη':>7s}  PR   η")
    for r in results:
        print(f"{r['name']:26s} {r['kind']:8s} "
              f"{r['PR_model']:6.3f}/{r['PR_meas']:.3f} "
              f"{r['dPR_rel']:+8.2%} "
              f"{r['eta_model']:6.3f}/{r['eta_meas']:.3f} "
              f"{100 * r['deta']:+6.1f}p  "
              f"{_fmt_status(r['PR_pass'])} {_fmt_status(r['eta_pass'])}")
    for a in anchors:
        print(f"{a['name']:26s} {'ancla':8s} {'':>14s} {'':>8s} {'':>14s} "
              f"{'':>7s}  {_fmt_status(a['ok'])}")

    if offd:
        print()
        print(f"{'Fuera de diseño (F-02)':34s} {'cantidad':13s} "
              f"{'modelo/medido':>16s} {'Δ':>8s}  obj")
        for o in offd:
            for it in o["items"]:
                print(f"{o['name']:34s} {it['key']:13s} "
                      f"{it['model']:7.2f}/{it['meas']:.2f} "
                      f"{it['d_rel']:+8.2%}  "
                      f"{_fmt_status(it['ok'])}"
                      + ("" if it["ok"] else
                         f"  (guarda ±{100 * it['guard']:.0f}%: "
                         f"{_fmt_status(it['ok_guard'])})"))

    if speed:
        print()
        print(f"{'Speedline medida (F-02)':30s} {'métrica':14s} "
              f"{'modelo/medido':>20s} {'Δ':>9s}  obj")
        for sp in speed:
            for it in _speedline_items(sp):
                print(f"{sp['name'][:30]:30s} {it['key']:14s} "
                      f"{it['txt']:>20s} {it['d_txt']:>9s}  "
                      f"{_fmt_status(it['ok'])}"
                      + ("" if it["ok"] else
                         f"  (guarda {it['guard']:.2f}: "
                         f"{_fmt_status(it['ok_guard'])})"))

    if not args.no_write:
        out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "RESULTS.md")
        write_results_md(results, anchors, out, offd, speed)
        print(f"\n→ {out}")

    # Fuera de diseño y speedline: --strict usa la GUARDA interina, igual
    # que η, para no dejar el CI rojo de forma permanente mientras la
    # brecha esté documentada y abierta. El objetivo se reporta siempre.
    hard_fail = (any(not r["PR_pass"] or not r["eta_guard"] for r in results)
                 or any(not a["ok"] for a in anchors)
                 or any(not it["ok_guard"] for o in offd
                        for it in o["items"])
                 or any(not it["ok_guard"] for sp in speed
                        for it in _speedline_items(sp)))
    if hard_fail:
        print("\nSTRICT: fallo duro (PR, guarda de η / fuera de diseño / "
              "speedline, o ancla de regresión).")
    return 1 if (args.strict and hard_fail) else 0


if __name__ == "__main__":
    sys.exit(main())
