"""
QUASAR Phy-AC · structures_core.py
==================================
Capa 1s del CEM: NÚCLEO ESTRUCTURAL del rotor axial (nivel L0s — analítico,
NumPy puro, <1 ms por evaluación), espejo de la escalera de fidelidad
aerodinámica (L0s analítico → L1s FEA axisimétrico → L2s FEA 3D externo).

Contrato (idéntico a Phy-CC):

    evaluate_structural(theta, record, material="Ti-6Al-4V") -> dict

    theta   vector de diseño 13-D — el mismo del núcleo aerodinámico.
    record  record meanline L0 del MISMO theta (usa stage_table para radios,
            cuerdas, nº de álabes y T0 por etapa; el acoplamiento térmico
            entra por aquí, no se recalcula).

Contenido (L0s):

  * Biblioteca de MATERIALES con derating térmico piecewise-lineal +
    límite AN² por material (práctica de diseño de discos de compresor).
  * Solver de DISCO ROTATORIO 1D (formulación de desplazamiento, espesor
    variable, Thomas O(n)) heredado VERBATIM de Phy-CC — validado contra
    la solución exacta de Timoshenko (anular, espesor constante).
  * Modelo de rotor axial: cada etapa se evalúa como disco anular (alma de
    espesor = fracción de la cuerda, del barreno al hub) con la tracción
    centrífuga de sus álabes en el rim; se reporta la PEOR etapa (las
    frontales tienen más radio/masa de álabe, las traseras más T).
  * Márgenes de producto: fluencia a sobrevelocidad, margen de estallido
    (esfuerzo tangencial promedio), raíz de álabe (columna centrífuga con
    factor de estrechamiento), límite AN².

El vector g_struct <= 0 es RESTRICCIÓN DURA del optimizador — DesignSpec
(neural_optimizer) lo concatena a las restricciones aero/espec.
"""

from __future__ import annotations

import math

import numpy as np

# ==========================================================================
# 1. Biblioteca de materiales
# ==========================================================================
# Propiedades a temperatura ambiente + curva de derating k(T) aplicada a
# sigma_y y sigma_uts. an2_max en in²·rpm² (unidad de la práctica
# industrial); IN2_TO_M2 lo convierte al SI del record. Valores típicos de
# handbook (MMPDS/ASM): suficientes para el proxy L0s.
IN2_TO_M2 = 6.4516e-4

MATERIALS = {
    "Al-2618-T61": dict(
        rho=2760.0, E=74.5e9, nu=0.33, alpha=22.0e-6,
        sigma_y=372e6, sigma_uts=441e6, an2_max_in2rpm2=2.0e10,
        derate=[(293, 1.00), (423, 0.90), (473, 0.65), (523, 0.40),
                (573, 0.25)],
        note="aleación ligera; solo etapas frontales frías"),
    "Al-7075-T6": dict(
        rho=2810.0, E=71.7e9, nu=0.33, alpha=23.4e-6,
        sigma_y=503e6, sigma_uts=572e6, an2_max_in2rpm2=2.5e10,
        derate=[(293, 1.00), (373, 0.95), (423, 0.75), (473, 0.45),
                (523, 0.25)],
        note="alta resistencia en frío; derating térmico severo"),
    "Ti-6Al-4V": dict(
        rho=4430.0, E=113.8e9, nu=0.342, alpha=8.6e-6,
        sigma_y=880e6, sigma_uts=950e6, an2_max_in2rpm2=4.5e10,
        derate=[(293, 1.00), (473, 0.82), (573, 0.72), (673, 0.62)],
        note="default del proyecto (discos y álabes de compresor)"),
    "Inconel-718": dict(
        rho=8190.0, E=200.0e9, nu=0.29, alpha=13.0e-6,
        sigma_y=1034e6, sigma_uts=1240e6, an2_max_in2rpm2=4.0e10,
        derate=[(293, 1.00), (673, 0.94), (923, 0.88), (1023, 0.70)],
        note="etapas traseras calientes; denso (castiga centrífugo)"),
    "17-4PH-H900": dict(
        rho=7750.0, E=196.0e9, nu=0.27, alpha=10.8e-6,
        sigma_y=1170e6, sigma_uts=1310e6, an2_max_in2rpm2=3.5e10,
        derate=[(293, 1.00), (573, 0.85), (673, 0.75)],
        note="acero PH; común en compresores industriales"),
}

# --- Parámetros del modelo (calibrables; documentar cambios) --------------
BORE_RATIO = 0.30         # radio del barreno / radio de hub de la etapa
WEB_T_RATIO = 0.50        # espesor del alma del disco / cuerda del rotor
RIM_T_RATIO = 1.10        # espesor del rim / espesor del alma
BLADE_TAPER_K = 0.55      # factor de estrechamiento de la raíz (álabe cónico)
BLADE_T_OVER_C = 0.08     # espesor máx/cuerda para la masa del álabe
OVERSPEED = 1.05          # sobrevelocidad de verificación de fluencia
BURST_MARGIN_MIN = 1.22   # N_burst/N_max mínimo (práctica API 617 / aero)
BURST_UTIL = 0.95         # utilización de sigma_uts en el criterio de burst
T_METAL_FRAC = 1.00       # T_metal = T0 de salida de la etapa (conservador)

N_STRUCT_CONSTRAINTS = 4


def list_materials() -> list[str]:
    return list(MATERIALS)


def material_at_T(name: str, T_K: float) -> dict:
    """Propiedades del material con sigma_y/sigma_uts derateadas a T_K."""
    if name not in MATERIALS:
        raise ValueError(f"material desconocido '{name}' — disponibles: "
                         f"{', '.join(MATERIALS)}")
    m = dict(MATERIALS[name])
    pts = m["derate"]
    Ts = np.array([p[0] for p in pts])
    ks = np.array([p[1] for p in pts])
    k = float(np.interp(T_K, Ts, ks))       # extrapolación plana en extremos
    m["derate_factor"] = k
    m["sigma_y_T"] = m["sigma_y"] * k
    m["sigma_uts_T"] = m["sigma_uts"] * k
    return m


# ==========================================================================
# 2. Solver de disco rotatorio 1D (heredado verbatim de Phy-CC)
# ==========================================================================
def _thomas(lo: np.ndarray, di: np.ndarray, up: np.ndarray,
            rhs: np.ndarray) -> np.ndarray:
    """Algoritmo de Thomas para sistema tridiagonal (lo[0] y up[-1] se
    ignoran). O(n) — g_struct corre dentro del lazo del optimizador."""
    n = len(di)
    cp = np.empty(n)
    dp = np.empty(n)
    cp[0] = up[0] / di[0]
    dp[0] = rhs[0] / di[0]
    for i in range(1, n):
        m = di[i] - lo[i] * cp[i - 1]
        cp[i] = up[i] / m if i < n - 1 else 0.0
        dp[i] = (rhs[i] - lo[i] * dp[i - 1]) / m
    u = np.empty(n)
    u[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        u[i] = dp[i] - cp[i] * u[i + 1]
    return u


def solve_disk(r: np.ndarray, h: np.ndarray, rho: float, E: float,
               nu: float, omega: float, sigma_rim: float = 0.0):
    """Disco axisimétrico rotatorio en formulación de desplazamiento.

    Equilibrio (plane stress, espesor h(r) variable):
        d/dr [ h·r·σ_r ] − h·σ_θ + ρ·ω²·h·r² = 0
        σ_r = C·(u' + ν·u/r),  σ_θ = C·(u/r + ν·u'),  C = E/(1−ν²)

    Contornos: barreno libre (σ_r=0) y rim con tracción σ_r = sigma_rim
    (tirón centrífugo de los álabes). Devuelve (u, σ_r, σ_θ) en SI.
    Validación: reproduce Timoshenko (disco anular, h const) a <1%.
    """
    r = np.asarray(r, dtype=float)
    h = np.asarray(h, dtype=float)
    n = len(r)
    d = r[1] - r[0]
    C = E / (1.0 - nu ** 2)
    hr = h * r

    lo = np.zeros(n)
    di = np.zeros(n)
    up = np.zeros(n)
    rhs = np.zeros(n)
    Ap = 0.5 * (hr[1:-1] + hr[2:])
    Am = 0.5 * (hr[1:-1] + hr[:-2])
    lo[1:-1] = C * Am / d ** 2 - C * nu * h[:-2] / (2 * d) \
        + C * nu * h[1:-1] / (2 * d)
    di[1:-1] = -C * (Ap + Am) / d ** 2 - C * h[1:-1] / r[1:-1]
    up[1:-1] = C * Ap / d ** 2 + C * nu * h[2:] / (2 * d) \
        - C * nu * h[1:-1] / (2 * d)
    rhs[1:-1] = -rho * omega ** 2 * h[1:-1] * r[1:-1] ** 2

    b0 = -3.0 / (2 * d) + nu / r[0]
    c0 = 4.0 / (2 * d)
    e0 = -1.0 / (2 * d)
    di[0] = b0 - e0 * lo[1] / up[1]
    up[0] = c0 - e0 * di[1] / up[1]
    rhs[0] = -e0 * rhs[1] / up[1]

    b_n = C * (3.0 / (2 * d) + nu / r[-1])
    a_n = C * (-4.0 / (2 * d))
    e_n = C * (1.0 / (2 * d))
    lo[-1] = a_n - e_n * di[-2] / lo[-2]
    di[-1] = b_n - e_n * up[-2] / lo[-2]
    rhs[-1] = sigma_rim - e_n * rhs[-2] / lo[-2]

    u = _thomas(lo, di, up, rhs)
    up_grad = np.gradient(u, r)
    sigma_r = C * (up_grad + nu * u / r)
    sigma_t = C * (u / r + nu * up_grad)
    return u, sigma_r, sigma_t


# ==========================================================================
# 3. Evaluación estructural completa (L0s) — peor etapa del tambor
# ==========================================================================
def _stage_disk(stage: dict, omega: float, rho: float, E: float,
                nu: float) -> dict:
    """Disco anular de una etapa: alma (h const) del barreno al hub, con
    la tracción de los álabes en el rim. Devuelve esfuerzos clave en SI."""
    r_hub = stage["r_hub_mm"] / 1000.0
    r_tip = stage["r_tip_mm"] / 1000.0
    c = stage["chord_rotor_mm"] / 1000.0
    n_bl = int(stage["n_blades_rotor"])
    h_blade = r_tip - r_hub

    # Tirón centrífugo de los álabes → tracción en el rim del disco
    A_cross = BLADE_T_OVER_C * c * c            # sección media del álabe
    vol_blade = BLADE_TAPER_K * A_cross * h_blade
    r_cg = 0.5 * (r_hub + r_tip)
    F_blades = n_bl * rho * vol_blade * omega ** 2 * r_cg
    h_rim = RIM_T_RATIO * WEB_T_RATIO * c
    sigma_pull = F_blades / max(2 * math.pi * r_hub * h_rim, 1e-9)

    r_bore = max(BORE_RATIO * r_hub, 0.004)
    n = 120
    r = np.linspace(r_bore, r_hub, n)
    h = np.full(n, WEB_T_RATIO * c)
    u, s_r, s_t = solve_disk(r, h, rho, E, nu, omega, sigma_rim=sigma_pull)
    s_vm = np.sqrt(s_r ** 2 - s_r * s_t + s_t ** 2)
    s_t_avg = float(np.trapezoid(s_t * h, r) / np.trapezoid(h, r))

    # Raíz de álabe: columna centrífuga con estrechamiento
    sigma_root = BLADE_TAPER_K * rho * omega ** 2 * \
        0.5 * (r_tip ** 2 - r_hub ** 2)

    return dict(sigma_vm_max=float(np.max(s_vm)),
                sigma_bore=float(s_t[0]), sigma_t_avg=s_t_avg,
                sigma_pull=sigma_pull, sigma_root=sigma_root,
                u_rim=float(u[-1]),
                disk_mass=rho * float(np.trapezoid(2 * math.pi * r * h, r)))


def evaluate_structural(theta: np.ndarray, record: dict,
                        material: str = "Ti-6Al-4V",
                        overspeed: float = OVERSPEED) -> dict:
    """Record estructural del rotor axial: evalúa cada etapa (disco + raíz
    de álabe a su temperatura local) y reporta la peor. Los esfuerzos
    escalan con ω² (elasticidad lineal): la verificación a sobrevelocidad
    multiplica por overspeed².
    """
    theta = np.asarray(theta, dtype=float)
    RPM = float(theta[1])
    T0_in = float(theta[10])
    omega = RPM * 2 * math.pi / 60.0
    st = record["stage_table"]

    worst = dict(yield_ratio=-1.0, burst_margin=99.0, blade_ratio=-1.0)
    per_stage = []
    total_mass = 0.0
    for s in st:
        T_metal = T0_in + T_METAL_FRAC * (s["T0_out_K"] - T0_in)
        mat = material_at_T(material, T_metal)
        d = _stage_disk(s, omega, mat["rho"], mat["E"], mat["nu"])
        yr = d["sigma_vm_max"] * overspeed ** 2 / mat["sigma_y_T"]
        bm = math.sqrt(BURST_UTIL * mat["sigma_uts_T"]
                       / max(d["sigma_t_avg"], 1.0))
        br = d["sigma_root"] * overspeed ** 2 / mat["sigma_y_T"]
        total_mass += d["disk_mass"]
        per_stage.append(dict(stage=s["stage"], T_metal_K=T_metal,
                              sigma_vm_MPa=d["sigma_vm_max"] / 1e6,
                              sigma_root_MPa=d["sigma_root"] / 1e6,
                              yield_ratio=yr, burst_margin=bm,
                              blade_ratio=br))
        if yr > worst["yield_ratio"]:
            worst.update(yield_ratio=yr, stage_yield=s["stage"],
                         sigma_vm_MPa=d["sigma_vm_max"] / 1e6,
                         sigma_allow_MPa=mat["sigma_y_T"] / 1e6,
                         derate_factor=mat["derate_factor"])
        worst["burst_margin"] = min(worst["burst_margin"], bm)
        worst["blade_ratio"] = max(worst["blade_ratio"], br)

    # Límite AN² (annulus de salida, la métrica de la práctica industrial)
    AN2 = float(record.get("AN2_m2_rpm2", 0.0))
    mat_cold = MATERIALS[material]
    AN2_max = mat_cold["an2_max_in2rpm2"] * IN2_TO_M2

    warnings: list[str] = []
    if worst.get("derate_factor", 1.0) < 0.6:
        warnings.append(
            f"{material} derateado a {worst['derate_factor']:.0%} en la "
            f"etapa {worst.get('stage_yield', '?')} — considerar otro "
            "material para las etapas traseras")

    # --- Vector g_struct <= 0 (continuo, dominancia de Deb) ---------------
    g_struct = [
        worst["yield_ratio"] - 1.0,             # fluencia a sobrevelocidad
        BURST_MARGIN_MIN - worst["burst_margin"],  # margen de estallido
        AN2 / AN2_max - 1.0,                    # límite AN²
        worst["blade_ratio"] - 1.0,             # raíz de álabe
    ]

    return dict(
        material=material,
        sigma_vm_max_MPa=worst.get("sigma_vm_MPa", 0.0),
        sigma_allow_MPa=worst.get("sigma_allow_MPa", 0.0),
        yield_ratio_overspeed=worst["yield_ratio"],
        stage_worst_yield=worst.get("stage_yield", 0),
        burst_margin=worst["burst_margin"],
        burst_margin_min=BURST_MARGIN_MIN,
        blade_ratio=worst["blade_ratio"],
        AN2_m2_rpm2=AN2, AN2_max_m2_rpm2=AN2_max,
        AN2_in2rpm2=AN2 / IN2_TO_M2,
        rotor_mass_kg=total_mass,
        per_stage=per_stage,
        g_struct=[float(g) for g in g_struct],
        feasible_struct=bool(all(g <= 0.0 for g in g_struct)),
        structural_warnings=warnings,
        params=dict(overspeed=overspeed, bore_ratio=BORE_RATIO,
                    web_t_ratio=WEB_T_RATIO, blade_taper_k=BLADE_TAPER_K,
                    burst_util=BURST_UTIL),
        source="disk1d_L0s",
    )


# ==========================================================================
# 4. Solución analítica exacta (validación del solver — tests)
# ==========================================================================
def analytic_annular_disk(r: np.ndarray, Rb: float, Ro: float, rho: float,
                          nu: float, omega: float):
    """Timoshenko: disco anular de espesor constante, bordes libres.
        σ_r = (3+ν)/8·ρω²·(Rb² + Ro² − Rb²Ro²/r² − r²)
        σ_θ = (3+ν)/8·ρω²·(Rb² + Ro² + Rb²Ro²/r² − (1+3ν)/(3+ν)·r²)
    """
    k = (3.0 + nu) / 8.0 * rho * omega ** 2
    s_r = k * (Rb ** 2 + Ro ** 2 - (Rb * Ro / r) ** 2 - r ** 2)
    s_t = k * (Rb ** 2 + Ro ** 2 + (Rb * Ro / r) ** 2
               - (1.0 + 3.0 * nu) / (3.0 + nu) * r ** 2)
    return s_r, s_t


if __name__ == "__main__":
    import sys
    import time
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    from physics_core import evaluate, Fidelity
    theta_ref = np.array([4.0, 12_500.0, 0.62, 0.55, 0.30, -0.10,
                          0.60, 1.20, 1.10, 2.20, 288.15, 101_325.0, 25.0])
    rec = evaluate(theta_ref, fidelity=Fidelity.L0)
    t0 = time.perf_counter()
    for mat in list_materials():
        s = evaluate_structural(theta_ref, rec, material=mat)
        print(f"{mat:14s} σvm={s['sigma_vm_max_MPa']:7.1f} MPa "
              f"(σy_T={s['sigma_allow_MPa']:6.1f})  "
              f"burst={s['burst_margin']:.2f}  "
              f"AN²={s['AN2_in2rpm2']:.2e} in²rpm²  "
              f"raíz={s['blade_ratio']:.2f}  "
              f"{'OK' if s['feasible_struct'] else 'VIOLA'} "
              f"{s['structural_warnings']}")
    dt = (time.perf_counter() - t0) / len(MATERIALS)
    print(f"[{1000*dt:.2f} ms/eval]")
