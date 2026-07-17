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

from physics_core import CP, GAMMA, Fidelity, _meanline, evaluate
from validation.machines import (MACHINES, REGRESSION_ANCHORS, PR_TOL_REL,
                                 ETA_TOL_PTS, ETA_GUARD_PTS)

KIND_LABEL = {
    "rotor": "rotor aislado (t-t)",
    "stage": "etapa completa (t-t)",
    "machine": "máquina multietapa (t-t)",
}


# ---------------------------------------------------------------------------
# Espec publicada → θ 13-D
# ---------------------------------------------------------------------------
def measured_dT0(spec: dict, measured: dict) -> float:
    """ΔT0 total desde (PR, η) medidos — el trabajo real de la máquina."""
    PR, T0 = measured["PR"], spec["T0"]
    if "eta_isen" in measured:
        return T0 * (PR ** ((GAMMA - 1) / GAMMA) - 1.0) / measured["eta_isen"]
    # η politrópica: τ = PR^((γ−1)/(γ·η_p))
    return T0 * (PR ** ((GAMMA - 1) / (GAMMA * measured["eta_poly"])) - 1.0)


def build_theta(spec: dict, measured: dict) -> np.ndarray:
    """Construye θ que reproduce annulus (r_tip vía φ1) y trabajo (ψ_mid)."""
    omega = spec["RPM"] * 2 * math.pi / 60.0
    r_tip_target = spec["U_tip"] / omega
    r_mean = 0.5 * (1.0 + spec["HTR"]) * r_tip_target
    Um = omega * r_mean
    psi = CP * measured_dT0(spec, measured) / (spec["n_stages"] * Um ** 2)

    def theta_of(phi1: float) -> np.ndarray:
        return np.array([spec["n_stages"], spec["RPM"], spec["HTR"], phi1,
                         psi, 0.0, spec["Rx_est"], spec["sigma_r"],
                         spec["sigma_s"], spec["AR"], spec["T0"],
                         spec["P0"], spec["mdot"]])

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
        # PR/η del ROTOR de la etapa 1, derivados del desglose de pérdidas
        s = rec["stage_table"][0]
        L = s["losses"]
        dh0 = s["dh0"]
        eta_rot = float(np.clip(
            (dh0 - L["rotor"] - L["clearance"]) / max(dh0, 1e-3), 0.05, 0.99))
        tau = 1.0 + dh0 / (CP * T0)
        PR = (1.0 + eta_rot * dh0 / (CP * T0)) ** (GAMMA / (GAMMA - 1))
    else:
        PR = rec["PR"]
        tau = rec["T0_out"] / T0
    eta_isen = (PR ** ((GAMMA - 1) / GAMMA) - 1.0) / max(tau - 1.0, 1e-9)
    eta_poly = ((GAMMA - 1) / GAMMA) * math.log(PR) / \
        max(math.log(tau), 1e-9) if PR > 1 else 0.0
    return dict(PR=PR, eta_isen=eta_isen, eta_poly=eta_poly)


def run_machine(m: dict) -> dict:
    # Inyecta la holgura de punta PUBLICADA de la máquina (ε absoluta en
    # mm — parámetro de módulo del meanline) y restaura siempre.
    import physics_core as _pc
    eps_saved = _pc.TIP_CLEARANCE_MM
    if "eps_tip_mm" in m:
        _pc.TIP_CLEARANCE_MM = float(m["eps_tip_mm"])
    try:
        theta = build_theta(m["spec"], m["measured"])
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


def write_results_md(results: list[dict], anchors: list[dict],
                     path: str) -> None:
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

    if not args.no_write:
        out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "RESULTS.md")
        write_results_md(results, anchors, out)
        print(f"\n→ {out}")

    hard_fail = (any(not r["PR_pass"] or not r["eta_guard"] for r in results)
                 or any(not a["ok"] for a in anchors))
    if hard_fail:
        print("\nSTRICT: fallo duro (PR, guarda de η o ancla de regresión).")
    return 1 if (args.strict and hard_fail) else 0


if __name__ == "__main__":
    sys.exit(main())
