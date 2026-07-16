"""
QUASAR Phy-AC · physics_core.py
===============================
Capa 1 del Computational Engineering Model: NÚCLEO FÍSICO MULTI-FIDELIDAD
para compresores AXIALES multietapa (hermano de Phy-CC, centrífugos).

Contrato único (idéntico a Phy-CC):

    evaluate(theta, fidelity=Fidelity.L0|L1) -> dict (PerformanceRecord)

Niveles de fidelidad
--------------------
  L0  Stage-stacking meanline 1-D (~1 ms/punto). Triángulos de velocidad
      medios por etapa desde (φ, ψ, Rx), pérdidas de perfil por Lieblein
      (espesor de momento de estela vs difusión equivalente), secundarias/
      endwall por Howell, holgura de punta, choque normal si M_rel > 1.
      Margen de bombeo por coeficiente de subida de presión estática
      (proxy de Koch 1981). Devuelve el VECTOR DE RESTRICCIONES g <= 0
      (8 aero) y banderas de régimen — la interfaz de la dominancia
      restringida de Deb en la capa de búsqueda.
  L1  Streamline Curvature Method de turbo-design (NASA TD3) con spool
      AXIAL multi-fila y los patches heredados de Phy-CC (frustum de área
      + loss 0-D). Import perezoso: sin turbo-design, L1 degrada a L0.
  L2  CFD 3D (externa). Hook de calibración register_hifi_pair(theta,
      y_hifi): corrección afín por salida sobre L0/L1.

Otras capacidades (heredadas del patrón Phy-CC)
-----------------------------------------------
  * Caché persistente por hash de theta+fidelity (JSONL, solo L1+).
  * evaluate_batch(..., n_workers) — paralelismo por procesos para L1.
  * physics_features(theta) — embedding físico cerrado para el surrogate:
    el ensemble aprende el RESIDUAL sobre la física, no la física desde 0.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import math
import os
import warnings
from concurrent.futures import ProcessPoolExecutor
from enum import IntEnum

import numpy as np

warnings.filterwarnings("ignore")

# ==========================================================================
# 0. Espacio de diseño (14-D) — distribuciones escalares, sin parámetros
#    libres por etapa (ver plan §L0 y docs/Science.md)
# ==========================================================================
DESIGN_VARS = [
    # nombre          min        max          unidad
    ("n_stages",      1.0,       8.0),        # -    (se redondea a entero)
    ("RPM",       5_000.0,  25_000.0),        # rev/min
    ("HTR_in",        0.40,      0.80),       # -    hub-to-tip ratio entrada
    ("phi1",          0.35,      0.80),       # -    coef. de flujo Cx/Um
    ("psi_mid",       0.22,      0.45),       # -    carga media Δh0/Um²
    ("psi_slope",    -0.30,      0.30),       # -    gradiente de ψ frontal→trasero
    ("Rx_mean",       0.50,      0.85),       # -    grado de reacción medio
    ("sigma_r",       0.90,      1.60),       # -    solidez rotor (media)
    ("sigma_s",       0.80,      1.50),       # -    solidez estátor (media)
    ("AR",            1.20,      3.50),       # -    aspect ratio rotor h/c
    ("T0_in",       270.0,     320.0),        # K
    ("P0_in",    85_000.0, 110_000.0),        # Pa
    ("massflow",      5.0,      60.0),        # kg/s
]
# NOTA: r_tip de entrada NO es variable de diseño — queda DERIVADO de la
# continuidad (dados φ1, HTR, RPM, ṁ el annulus de entrada es único). Con
# r_tip libre el espacio quedaba sobredeterminado: la continuidad imponía
# un φ distinto del pedido y toda la parametrización de Smith se rompía.
VAR_NAMES = [v[0] for v in DESIGN_VARS]
BOUNDS_LO = np.array([v[1] for v in DESIGN_VARS])
BOUNDS_HI = np.array([v[2] for v in DESIGN_VARS])
NDIM = len(DESIGN_VARS)

GAMMA, CP, RGAS = 1.4, 1005.0, 287.0
MU_AIR = 1.81e-5          # Pa·s

# --- Parámetros del meanline (calibrables; ver docs/VALIDATION.md) --------
TIP_CLEARANCE_RATIO = 0.015   # ε/h holgura de punta relativa (rotor)
K_TIP_CLEARANCE = 2.0         # Δη ≈ K·(ε/h) — sensibilidad típica 2-3 %/1 %
BLOCKAGE_INIT = 0.98          # KB de la primera etapa
BLOCKAGE_PER_STAGE = 0.005    # crecimiento de bloqueo endwall por etapa
ANNULUS_MODE = "const_mean"   # const_hub | const_mean | const_tip
STATOR_AR_FACTOR = 1.10       # AR_estator = 1.10 · AR_rotor
ROW_GAP_FRACTION = 0.25       # gap axial inter-fila / cuerda
CHORD_TAPER = 0.85            # c_tip/c_hub (usado por la capa 5a)
K_PROFILE = 1.0               # multiplicador de pérdida de perfil (calibración)
K_ENDWALL = 1.0               # multiplicador de pérdida secundaria/endwall
K_SHOCK = 0.70                # multiplicador de pérdida de choque. El
#                               choque real del pasaje es OBLICUO; el modelo
#                               de choque normal al Mach de entrada
#                               sobreestima. Calibrado contra NASA Rotor
#                               37/67 (validation/RESULTS.md) manteniendo
#                               Stage 35 dentro de tolerancia.
MX_CHOKE = 0.78               # Mach axial máx. antes de declarar choke local

# Work-done factor de Howell (mediado por nº de etapa; tabla clásica)
_HOWELL_WDF = [0.982, 0.952, 0.929, 0.910, 0.895, 0.883, 0.875, 0.868]

# Límites de las restricciones aero
DF_MAX_R = 0.55           # factor de difusión de Lieblein, rotor
DF_MAX_S = 0.60           # factor de difusión, estátor
DEHALLER_MIN = 0.72       # W2/W1 (rotor) y C3/C2 (estátor) mínimos
M_REL_TIP_MAX = 1.35      # Mach relativo en punta del primer rotor
CH_STALL_MAX = 0.48       # coef. de subida de presión estática máx. (Koch)
KOCH_SM_MIN = 0.10        # margen mínimo respecto a CH_STALL_MAX
M_EXIT_MAX = 0.55         # Mach absoluto de salida (amigable al difusor)
H_BLADE_MIN_MM = 8.0      # altura mínima del álabe de la última etapa

N_CONSTRAINTS = 8         # dim del vector g(theta) — solo aerodinámica


class Fidelity(IntEnum):
    L0 = 0   # stage-stacking meanline analítico
    L1 = 1   # SCM turbo-design (si está disponible)


def denormalize(u: np.ndarray) -> np.ndarray:
    """[0,1]^14 -> unidades físicas."""
    return BOUNDS_LO + np.asarray(u, dtype=float) * (BOUNDS_HI - BOUNDS_LO)


def normalize(theta: np.ndarray) -> np.ndarray:
    return (np.asarray(theta, dtype=float) - BOUNDS_LO) / (BOUNDS_HI - BOUNDS_LO)


def _theta_key(theta: np.ndarray, fidelity: int) -> str:
    arr = np.round(np.asarray(theta, dtype=float), 9)
    return hashlib.sha1(arr.tobytes() + bytes([fidelity])).hexdigest()


# ==========================================================================
# 1. PATCHES de turbo-design (heredados de Phy-CC, aplicados perezosamente)
# ==========================================================================
_TD_AVAILABLE: bool | None = None
_TD_REASON: str | None = None
_TD_MODS: dict = {}


def _patched_compute_streamline_areas(row):
    """Frustum cónico correcto: A = pi*(r_j + r_{j-1})*sqrt(dx^2+dr^2).

    El código original eleva dx al cuadrado y permite áreas negativas. En
    pasajes axiales (dr≈0) esta forma degenera correctamente al cilindro
    A = 2π·r·dx, así que el patch es seguro para ambos tipos de pasaje.
    """
    total_area = 0.0
    streamline_area = np.zeros(len(row.percent_hub_shroud))
    if len(row.percent_hub_shroud) <= 1:
        if getattr(row, "total_area", None):
            total_area = float(row.total_area)
            streamline_area = np.array([total_area])
        return total_area, streamline_area
    for j in range(1, len(row.percent_hub_shroud)):
        dx = row.x[j] - row.x[j - 1]
        dr = row.r[j] - row.r[j - 1]
        if abs(dx) < 1e-5 and abs(dr) > 1e-9:
            delta = abs(np.pi * (row.r[j] ** 2 - row.r[j - 1] ** 2))
        else:
            dl = math.sqrt(dx ** 2 + dr ** 2)
            delta = abs(np.pi * (row.r[j] + row.r[j - 1]) * dl)
        streamline_area[j] = delta
        total_area += delta
    return total_area, streamline_area


def _try_load_turbodesign() -> bool:
    """Import perezoso de turbo-design + patches. Nunca lanza: la ausencia
    de la librería degrada limpiamente el sistema a L0."""
    global _TD_AVAILABLE, _TD_REASON
    if _TD_AVAILABLE is not None:
        return _TD_AVAILABLE
    try:
        with contextlib.redirect_stdout(io.StringIO()), \
             contextlib.redirect_stderr(io.StringIO()):
            import turbodesign.flow_math as _fm
            import turbodesign.compressor_math as _cm
            import turbodesign.compressor_spool as _cs
            from turbodesign import Passage, PassageType, Inlet, Outlet
            from turbodesign.row_factory import make_rotor_row
            from turbodesign.compressor_spool import CompressorSpool
            from turbodesign.loss.losstype import LossBaseClass
            from turbodesign.loss.fixedpressureloss import FixedPressureLoss
            from turbodesign.enums import LossType
            from cantera import Solution
            # make_stator_row no existe en todas las versiones de TD3 1.4.x;
            # se introspecta y, si falta, el spool L1 corre solo con las
            # filas de rotor (ver _scm_solve).
            try:
                from turbodesign.row_factory import make_stator_row
            except Exception:
                make_stator_row = None
        _fm.compute_streamline_areas = _patched_compute_streamline_areas
        _cm.compute_streamline_areas = _patched_compute_streamline_areas

        # PATCH 3 (verificado en TD3 1.4.2, 2026-07-11): stator_calc hace
        # `row.Yp[:] = 0` pero con pérdida politrópica en filas de estátor
        # Yp puede llegar como np.float64 0-D → TypeError. Se coerce a
        # array 1-D antes de delegar. compressor_spool importa stator_calc
        # POR NOMBRE, así que hay que parchear ambas referencias.
        _orig_stator_calc = _cm.stator_calc

        def _patched_stator_calc(row, upstream, calculate_vm=True):
            if np.ndim(getattr(row, "Yp", None)) == 0:
                row.Yp = np.atleast_1d(np.asarray(row.Yp, dtype=float))
            return _orig_stator_calc(row, upstream, calculate_vm)

        _cm.stator_calc = _patched_stator_calc
        _cs.stator_calc = _patched_stator_calc

        class ScalarPolytropicLoss(LossBaseClass):
            """Eficiencia politrópica fija compatible con TD3 v1.4.2
            (rotor_calc hace float(loss_fn(...)) y solo acepta 0-D)."""

            def __init__(self, eta_poly: float):
                super().__init__(LossType.Polytropic)
                self._eta = float(eta_poly)

            def __call__(self, row, upstream):  # noqa: ARG002
                return np.float64(self._eta)

        _TD_MODS.update(dict(
            Passage=Passage, PassageType=PassageType, Inlet=Inlet,
            Outlet=Outlet, make_rotor_row=make_rotor_row,
            make_stator_row=make_stator_row,
            CompressorSpool=CompressorSpool, Solution=Solution,
            ScalarPolytropicLoss=ScalarPolytropicLoss,
            FixedPressureLoss=FixedPressureLoss,
        ))
        _TD_AVAILABLE = True
    except Exception as e:
        _TD_AVAILABLE = False
        msg = str(e).strip().splitlines()
        _TD_REASON = f"{type(e).__name__}: {msg[0] if msg else ''}"
    return _TD_AVAILABLE


def l1_available() -> bool:
    """True si turbo-design + cantera son importables y parchables."""
    return _try_load_turbodesign()


def l1_unavailable_reason() -> str | None:
    """Primera línea del error de import si L1 no está disponible."""
    _try_load_turbodesign()
    return _TD_REASON


# ==========================================================================
# 2. L0 — Stage-stacking meanline con correlaciones y restricciones
# ==========================================================================
def _normal_shock_p0_ratio(M: float) -> float:
    """P02/P01 a través de un choque normal (M>1)."""
    if M <= 1.0:
        return 1.0
    g = GAMMA
    t1 = ((g + 1) * M * M / ((g - 1) * M * M + 2)) ** (g / (g - 1))
    t2 = ((g + 1) / (2 * g * M * M - (g - 1))) ** (1 / (g - 1))
    return t1 * t2


def _lieblein_theta_c(deq: float) -> float:
    """Espesor de momento de estela θ/c vs difusión equivalente (Lieblein
    1959, fit clásico). Diverge al acercarse Deq→e^(1/0.95): se acota."""
    deq = min(max(deq, 1.0), 2.6)
    denom = 1.0 - 0.95 * math.log(deq)
    return 0.0045 / max(denom, 0.08)


def _shock_omega(M: float) -> float:
    """ω̄ de choque normal referido a la presión dinámica de entrada."""
    if M <= 1.0:
        return 0.0
    dyn = 1.0 - (1.0 + 0.5 * (GAMMA - 1) * M ** 2) ** (-GAMMA / (GAMMA - 1))
    return (1.0 - _normal_shock_p0_ratio(M)) / max(dyn, 1e-3)


def _row_losses(beta1: float, beta2: float, W1: float, W2: float,
                sigma: float, h_over_c: float, M_tip: float,
                M_mean: float) -> tuple[float, float, dict]:
    """Pérdidas de una fila en el marco relativo a la fila.

    beta1/beta2 en rad (ángulos de flujo entrada/salida respecto al eje,
    positivos), W1/W2 velocidades entrada/salida del marco de la fila,
    M_tip/M_mean Mach de entrada en punta y en la línea media (el choque
    se promedia entre ambos, práctica de Miller — el choque solo cubre el
    span exterior). Devuelve (omega_bar_total, dh_loss [J/kg], desglose).
    ω̄ está referida a la presión dinámica de ENTRADA de la fila.
    """
    cb1, cb2 = math.cos(beta1), math.cos(beta2)
    tb1, tb2 = math.tan(beta1), math.tan(beta2)
    # Difusión equivalente de Lieblein (forma de circulación, diseño)
    deq = (cb2 / max(cb1, 1e-3)) * (
        1.12 + 0.61 * (cb1 ** 2 / max(sigma, 1e-3)) * abs(tb1 - tb2))
    theta_c = _lieblein_theta_c(deq)
    # ω̄ de perfil (Lieblein): 2·(θ/c)·(σ/cosβ2)·(cosβ1/cosβ2)²
    om_profile = K_PROFILE * 2.0 * theta_c * (sigma / max(cb2, 1e-3)) * \
        (cb1 / max(cb2, 1e-3)) ** 2
    # Secundarias + annulus (Howell): CDs=0.018·CL², CDa=0.020·(s/h);
    # conversión arrastre→pérdida referida a la entrada:
    # ω̄ = CD·σ·cos²β1/cos³βm  (Dixon & Hall §3; Wm/W1 = cosβ1/cosβm)
    beta_m = math.atan(0.5 * (tb1 + tb2))
    cbm = math.cos(beta_m)
    CL = 2.0 * (1.0 / max(sigma, 1e-3)) * cbm * abs(tb1 - tb2)
    cd_sec = 0.018 * CL ** 2
    cd_ann = 0.020 * (1.0 / max(sigma, 1e-3)) / max(h_over_c, 0.2)
    om_endwall = K_ENDWALL * (cd_sec + cd_ann) * sigma * cb1 ** 2 / \
        max(cbm, 1e-3) ** 3
    # Choque: promedio punta/media (dos puntos del span, Miller-style)
    om_shock = K_SHOCK * 0.5 * (_shock_omega(M_tip) + _shock_omega(M_mean))
    om_total = om_profile + om_endwall + om_shock
    # ΔP0rel = ω̄·½ρW1² → pérdida de entalpía ≈ ΔP0rel/ρ = ω̄·½W1²
    dh = om_total * 0.5 * W1 ** 2
    detail = dict(profile=om_profile, endwall=om_endwall, shock=om_shock,
                  deq=deq, CL=CL)
    return om_total, dh, detail


def _solve_axial_mach(mdot: float, A: float, T0: float, P0: float,
                      kb: float, swirl_tan: float = 0.0) -> tuple[float, bool]:
    """Mach axial subsónico desde continuidad ṁ = ρ(M)·A·KB·Cx(M).

    swirl_tan = tanα del remolino en la estación (C = Cx·√(1+tan²α)).
    Devuelve (Mx, choked). Si no hay raíz subsónica devuelve el Mx del
    máximo de gasto (garganta) y choked=True — g queda continuo.
    """
    a0 = math.sqrt(GAMMA * RGAS * T0)
    rho0 = P0 / (RGAS * T0)
    k = 1.0 + swirl_tan ** 2

    def flow(mx: float) -> float:
        c = mx * a0  # aproximación: Mx referido a a0; corregido abajo
        # iteración corta densidad-estática
        T = T0
        for _ in range(4):
            C2 = k * c * c
            T = max(T0 - C2 / (2 * CP), 0.5 * T0)
            a = math.sqrt(GAMMA * RGAS * T)
            c = mx * a
        rho = rho0 * (T / T0) ** (1.0 / (GAMMA - 1))
        return rho * A * kb * c

    # el gasto crece hasta Mx≈1/√k; bisección sobre la rama subsónica
    mx_hi = min(0.999, 1.0 / math.sqrt(k))
    if flow(mx_hi) < mdot:
        return mx_hi, True
    lo, hi = 1e-4, mx_hi
    for _ in range(48):
        mid = 0.5 * (lo + hi)
        if flow(mid) < mdot:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), False


def _static_from_mach(mx: float, T0: float, P0: float,
                      swirl_tan: float = 0.0) -> tuple[float, float, float, float]:
    """(T, P, rho, Cx) estáticos desde Mach axial + remolino."""
    k = 1.0 + swirl_tan ** 2
    T = T0 / (1.0 + 0.5 * (GAMMA - 1) * k * mx * mx)
    a = math.sqrt(GAMMA * RGAS * T)
    Cx = mx * a
    P = P0 * (T / T0) ** (GAMMA / (GAMMA - 1))
    rho = P / (RGAS * T)
    return T, P, rho, Cx


def _meanline(theta: np.ndarray, rpm: float | None = None,
              mdot: float | None = None,
              frozen: dict | None = None) -> dict:
    """Stage-stacking 1-D en la línea media, termodinámicamente consistente.

    Punto de diseño: `_meanline(theta)` — los ángulos de flujo salen de
    (φ, ψᵢ, Rx) y el álabe se talla a incidencia cero (convención del
    generador de geometría). Fuera de diseño: pasar `rpm`/`mdot` y
    `frozen` (dict devuelto en record["frozen"] del diseño: ángulos de
    flujo y áreas congelados); la incidencia paga una pérdida parabólica.

    Convención de triángulos (Dixon, ángulos desde el eje, Cx constante
    en el diseño): Cu1/U = 1−Rx−ψ/2, Cu2/U = 1−Rx+ψ/2,
    tanβ = (U−Cu)/Cx, tanα = Cu/Cx.
    """
    (n_st_f, RPM_d, HTR, phi1, psi_mid, psi_slope, Rx,
     sigma_r, sigma_s, AR, T0_in, P0_in, mdot_d) = [float(x) for x in theta]
    # clamp de seguridad holgado (los bounds del optimizador son 1-8; la
    # validación evalúa máquinas de hasta 10+ etapas fuera de bounds)
    n_st = int(round(min(max(n_st_f, 1.0), 16.0)))
    RPM = float(rpm) if rpm is not None else RPM_d
    mdot = float(mdot) if mdot is not None else mdot_d
    off_design = frozen is not None
    clipped: list[str] = []

    omega = RPM * 2 * math.pi / 60.0
    a0_in = math.sqrt(GAMMA * RGAS * T0_in)
    rho0_in = P0_in / (RGAS * T0_in)

    # --- Annulus de entrada: r_tip DERIVADO de la continuidad --------------
    # Dados (φ1, HTR, RPM, ṁ): Cx = φ1·ω·r_mean y A = π·r_tip²(1−HTR²)
    # deben transportar ṁ — punto fijo en r_tip (único y contractivo).
    psi_first = psi_mid * (1.0 + psi_slope * (-1.0)) if n_st > 1 else psi_mid
    tan_a1_first = (1.0 - Rx - max(psi_first, 0.05) / 2.0) / phi1
    r_tip = 0.20
    for _ in range(80):
        r_mean_i = 0.5 * (1.0 + HTR) * r_tip
        Cx_i = phi1 * omega * r_mean_i
        C1_i = Cx_i * math.sqrt(1.0 + tan_a1_first ** 2)
        if C1_i > 0.95 * a0_in:
            C1_i = 0.95 * a0_in
        T1_i = max(T0_in - C1_i ** 2 / (2 * CP), 0.5 * T0_in)
        rho1_i = rho0_in * (T1_i / T0_in) ** (1.0 / (GAMMA - 1))
        A_req = mdot / max(rho1_i * Cx_i * BLOCKAGE_INIT, 1e-6)
        r_tip_new = math.sqrt(A_req / (math.pi * max(1.0 - HTR ** 2, 1e-3)))
        if abs(r_tip_new - r_tip) < 1e-9:
            r_tip = r_tip_new
            break
        r_tip = 0.5 * (r_tip + r_tip_new)
    if phi1 * omega * 0.5 * (1.0 + HTR) * r_tip * \
            math.sqrt(1.0 + tan_a1_first ** 2) > 0.95 * a0_in:
        clipped.append("inlet_C1_sonic")
    r_hub = HTR * r_tip
    r_mean = 0.5 * (r_hub + r_tip)          # línea media (const en const_mean)
    A1 = math.pi * (r_tip ** 2 - r_hub ** 2)
    if off_design:
        A1 = frozen["areas"][0]

    T0, P0 = T0_in, P0_in
    kb = BLOCKAGE_INIT
    stage_table: list[dict] = []
    frozen_out: dict = {"areas": [A1], "alpha1": [], "beta2": [],
                        "beta1_deg": [], "r_m": []}
    choke_any = False
    g_soft_pen = 0.0                        # penalización continua de choke
    length_m = 0.0
    U_tip1 = omega * r_tip
    M_rel_tip1 = 0.0
    dh_shaft_total = 0.0
    dh_loss_total = 0.0

    # Ángulo de entrada de la primera etapa (IGV implícito = pre-swirl del
    # triángulo repetitivo de la etapa 1; su pérdida se carga a la etapa 1).
    for i in range(n_st):
        U = omega * (frozen["r_m"][i] if off_design and i < len(frozen["r_m"])
                     else r_mean)
        wdf = _HOWELL_WDF[min(i, len(_HOWELL_WDF) - 1)]
        if n_st > 1:
            psi_i = psi_mid * (1.0 + psi_slope * (2.0 * i / (n_st - 1) - 1.0))
        else:
            psi_i = psi_mid
        psi_i = max(psi_i, 0.05)

        if off_design:
            tan_a1 = frozen["alpha1"][i]
            tan_b2_fr = frozen["beta2"][i]
            A_in = frozen["areas"][i]
        else:
            tan_a1 = (1.0 - Rx - psi_i / 2.0) / phi1
            tan_b2_fr = None
            A_in = frozen_out["areas"][i]

        # Continuidad en la entrada de la etapa (con remolino α1)
        mx, chk = _solve_axial_mach(mdot, A_in, T0, P0, kb, swirl_tan=tan_a1)
        if chk or mx > MX_CHOKE:
            choke_any = True
            g_soft_pen += max(mx / MX_CHOKE - 1.0, 0.0) + (0.5 if chk else 0.0)
            clipped.append(f"choke_st{i}")
        T1, P1, rho1, Cx = _static_from_mach(mx, T0, P0, swirl_tan=tan_a1)
        phi = Cx / max(U, 1e-3)

        if off_design:
            # geometría congelada: α1 y β2 metálicos fijos → ψ emergente
            tan_b2 = tan_b2_fr
            psi_eff = 1.0 - phi * (tan_a1 + tan_b2)
            psi_eff = float(np.clip(psi_eff, -0.2, 0.9))
            Cu1 = Cx * tan_a1
            Cu2 = U - Cx * tan_b2
        else:
            Cu1 = (1.0 - Rx - psi_i / 2.0) * U
            Cu2 = (1.0 - Rx + psi_i / 2.0) * U
            psi_eff = psi_i
            tan_b2 = (U - Cu2) / max(Cx, 1e-3)

        tan_b1 = (U - Cu1) / max(Cx, 1e-3)
        tan_a2 = Cu2 / max(Cx, 1e-3)
        b1, b2r = math.atan(tan_b1), math.atan(tan_b2)
        # α1 con el Cx REAL de continuidad (no el φ1 de diseño): mantiene
        # exacta la identidad tanβ1 + tanα1 = 1/φ del triángulo
        a1r, a2r = math.atan2(Cu1, max(Cx, 1e-3)), math.atan(tan_a2)
        W1 = Cx / max(math.cos(b1), 1e-3)
        W2 = Cx / max(math.cos(b2r), 1e-3)
        C1 = Cx / max(math.cos(a1r), 1e-3)
        C2 = Cx / max(math.cos(a2r), 1e-3)

        a1_sound = math.sqrt(GAMMA * RGAS * T1)
        # Punta local de esta etapa; Mach relativo con free-vortex Cu·r=const
        h_blade = A_in / (2 * math.pi * r_mean)
        r_tip_i = r_mean + 0.5 * h_blade
        r_hub_i = r_mean - 0.5 * h_blade
        U_t = omega * r_tip_i
        Cu1_t = Cu1 * r_mean / max(r_tip_i, 1e-6)
        W1_tip = math.sqrt(Cx ** 2 + (U_t - Cu1_t) ** 2)
        M_rel_tip = W1_tip / a1_sound
        if i == 0:
            M_rel_tip1 = M_rel_tip
            U_tip1 = U_t
        M1_rel_mean = W1 / a1_sound

        # Reacción de hub (free-vortex) — diagnóstico
        Cu1_h = Cu1 * r_mean / max(r_hub_i, 1e-6)
        Cu2_h = Cu2 * r_mean / max(r_hub_i, 1e-6)
        U_h = omega * r_hub_i
        Rx_hub = 1.0 - (Cu1_h + Cu2_h) / max(2.0 * U_h, 1e-6)

        # Trabajo de Euler (ψ ya es la carga LOGRADA; el wdf entra al tallar
        # el álabe con más giro en la capa 5a, no al trabajo aquí)
        dh0 = psi_eff * U * U
        dwu = abs(Cu2 - Cu1)

        # Factores de difusión de Lieblein
        DF_r = 1.0 - W2 / max(W1, 1e-3) + dwu / max(2.0 * sigma_r * W1, 1e-3)
        # estátor: entra C2 (α2), sale C1 de la etapa siguiente (α1)
        DF_s = 1.0 - C1 / max(C2, 1e-3) + dwu / max(2.0 * sigma_s * C2, 1e-3)
        dehaller_r = W2 / max(W1, 1e-3)
        dehaller_s = C1 / max(C2, 1e-3)

        # Geometría de fila (cuerdas, nº de álabes, longitud)
        c_rotor = h_blade / max(AR, 0.5)
        c_stator = h_blade / max(STATOR_AR_FACTOR * AR, 0.5)
        n_bl_r = max(int(round(2 * math.pi * r_mean * sigma_r / max(c_rotor, 1e-4))), 7)
        n_bl_s = max(int(round(2 * math.pi * r_mean * sigma_s / max(c_stator, 1e-4))), 7)
        length_m += (1.0 + ROW_GAP_FRACTION) * (c_rotor + c_stator)

        h_c_r = h_blade / max(c_rotor, 1e-4)
        h_c_s = h_blade / max(c_stator, 1e-4)

        # Pérdidas por fila
        om_r, dh_r, det_r = _row_losses(b1, b2r, W1, W2, sigma_r, h_c_r,
                                        M_rel_tip, M1_rel_mean)
        M2_abs = C2 / a1_sound
        om_s, dh_s, det_s = _row_losses(a2r, a1r, C2, C1, sigma_s, h_c_s,
                                        M2_abs, M2_abs)
        # Fuera de diseño: bucket parabólico de incidencia sobre el rotor
        inc_deg = 0.0
        if off_design:
            b1_design = frozen["beta1_deg"][i]
            inc_deg = math.degrees(b1) - b1_design
            bucket = 1.0 + (inc_deg / 10.0) ** 2
            dh_r *= min(bucket, 4.0)
        # Holgura de punta (solo rotor)
        dh_cl = K_TIP_CLEARANCE * TIP_CLEARANCE_RATIO * dh0 if dh0 > 0 else 0.0

        dh_loss_i = dh_r + dh_s + dh_cl
        eta_tt = float(np.clip((dh0 - dh_loss_i) / max(dh0, 1e-3), 0.05, 1.0)) \
            if dh0 > 1.0 else 0.05
        if dh0 <= 1.0:
            clipped.append(f"psi_collapse_st{i}")

        # Subida de presión estática de la etapa → proxy de Koch
        T0_out_i = T0 + dh0 / CP
        PR_i = (1.0 + eta_tt * dh0 / (CP * T0)) ** (GAMMA / (GAMMA - 1)) \
            if dh0 > 0 else 1.0
        P0_out_i = P0 * PR_i
        # Coef. de subida de presión estática de la etapa (proxy de Koch
        # 1981): Ch = Δh_estática / ½(W1²+C2²). En etapa repetitiva C3=C1
        # ⇒ Δh_estática = Δh0 (la energía cinética no cambia). El límite
        # CH_STALL_MAX≈0.48 reproduce el rango de stall de Koch (0.45-0.55).
        dh_static = dh0
        Ch = dh_static / max(0.5 * (W1 ** 2 + C2 ** 2), 1e-3)
        SM_i = 1.0 - Ch / CH_STALL_MAX

        # Estado estático y área en la SALIDA DEL ROTOR (estación 2) — lo
        # necesita el pasaje L1 para contraer el annulus fila a fila (sin
        # esto TD3 ve el área de entrada en el rotor y sobregira el flujo).
        T02_i = T0 + dh0 / CP
        T2_st = max(T02_i - C2 ** 2 / (2 * CP), 0.5 * T0)
        eta_r_approx = 0.9
        P02_i = P0 * (1.0 + eta_r_approx * dh0 / (CP * T0)) ** (GAMMA / (GAMMA - 1)) \
            if dh0 > 0 else P0
        P2_st = P02_i * (T2_st / T02_i) ** (GAMMA / (GAMMA - 1))
        rho2_st = max(P2_st / (RGAS * T2_st), 1e-3)
        A_rotor_exit = mdot / max(rho2_st * Cx * kb, 1e-6)

        stage_table.append(dict(
            stage=i, phi=phi, psi=psi_eff, Rx=Rx, Rx_hub=Rx_hub, U_m=U,
            A_rotor_exit_m2=A_rotor_exit, A_in_m2=A_in,
            r_mean_mm=r_mean * 1000, r_hub_mm=r_hub_i * 1000,
            r_tip_mm=r_tip_i * 1000, h_blade_mm=h_blade * 1000,
            alpha1_deg=math.degrees(a1r), alpha2_deg=math.degrees(a2r),
            beta1_deg=math.degrees(b1), beta2_deg=math.degrees(b2r),
            W1=W1, W2=W2, C1=C1, C2=C2, Cx=Cx,
            DF_rotor=DF_r, DF_stator=DF_s,
            dehaller_rotor=dehaller_r, dehaller_stator=dehaller_s,
            M_rel_tip=M_rel_tip, M1_rel_mean=M1_rel_mean, Mx=mx,
            chord_rotor_mm=c_rotor * 1000, chord_stator_mm=c_stator * 1000,
            n_blades_rotor=n_bl_r, n_blades_stator=n_bl_s,
            eta_tt=eta_tt, PR=PR_i, T0_in_K=T0, T0_out_K=T0_out_i,
            P0_in_Pa=P0, P0_out_Pa=P0_out_i, dh0=dh0, Ch=Ch, SM=SM_i,
            incidence_deg=inc_deg, wdf=wdf,
            losses=dict(rotor=dh_r, stator=dh_s, clearance=dh_cl,
                        omega_rotor=om_r, omega_stator=om_s,
                        deq_rotor=det_r["deq"], deq_stator=det_s["deq"]),
        ))

        dh_shaft_total += dh0
        dh_loss_total += dh_loss_i
        T0, P0 = T0_out_i, P0_out_i
        kb = max(kb - BLOCKAGE_PER_STAGE, 0.90)

        # Annulus de la siguiente estación por continuidad (Cx de diseño
        # constante; en const_mean, r_mean fijo y h se ajusta)
        if not off_design:
            # Área de la siguiente estación que MANTIENE el Cx de diseño
            # (φ1·Um), evaluada con el estado (T0,P0) tras la etapa. Para la
            # última estación el remolino ya fue quitado por el OGV.
            tan_a_next = tan_a1 if i < n_st - 1 else 0.0
            Cx_dsg = phi1 * omega * r_mean
            C_n = Cx_dsg * math.sqrt(1.0 + tan_a_next ** 2)
            T_n = max(T0 - C_n ** 2 / (2 * CP), 0.5 * T0)
            rho_n = (P0 / (RGAS * T0)) * (T_n / T0) ** (1.0 / (GAMMA - 1))
            A_next = mdot / max(rho_n * Cx_dsg * kb, 1e-6)
            frozen_out["areas"].append(A_next)
            frozen_out["alpha1"].append(tan_a1)
            frozen_out["beta2"].append(tan_b2)
            frozen_out["beta1_deg"].append(math.degrees(b1))
            frozen_out["r_m"].append(r_mean)
            if ANNULUS_MODE == "const_hub":
                h_n = A_next / (2 * math.pi)  # resolver r_tip
                # r_tip² = r_hub² + A/π
                r_tip_n = math.sqrt(r_hub ** 2 + A_next / math.pi)
                r_mean = 0.5 * (r_hub + r_tip_n)
            elif ANNULUS_MODE == "const_tip":
                r_hub_n = math.sqrt(max(r_tip ** 2 - A_next / math.pi, 1e-8))
                r_mean = 0.5 * (r_hub_n + r_tip)
            # const_mean: r_mean no cambia
        else:
            if i + 1 < len(frozen["areas"]):
                pass  # áreas congeladas ya en frozen["areas"]

    # --- OGV: quita el remolino residual α1 → axial (fila de estátor) -----
    s_last = stage_table[-1]
    a_ogv = math.radians(s_last["alpha1_deg"])
    C_ogv_in = s_last["C1"]
    Cx_last = s_last["Cx"]
    om_ogv, _, _ = _row_losses(a_ogv, 0.0, C_ogv_in, Cx_last, sigma_s,
                               s_last["h_blade_mm"] /
                               max(s_last["chord_stator_mm"], 1e-3),
                               0.5, 0.5)
    T_st, _, rho_st, _ = _static_from_mach(s_last["Mx"], T0, P0,
                                           swirl_tan=math.tan(a_ogv))
    P0 = max(P0 - om_ogv * 0.5 * rho_st * C_ogv_in ** 2, 0.5 * P0)
    length_m += (1.0 + ROW_GAP_FRACTION) * s_last["chord_stator_mm"] / 1000.0
    # estado de salida tras OGV (flujo axial)
    mx_e, chk_e = _solve_axial_mach(mdot, frozen_out["areas"][-1] if not off_design
                                    else frozen["areas"][-1], T0, P0, kb, 0.0)
    if chk_e:
        choke_any = True
        clipped.append("choke_exit")
    T_e, P_e, rho_e, Cx_e = _static_from_mach(mx_e, T0, P0, 0.0)
    M_exit = Cx_e / math.sqrt(GAMMA * RGAS * T_e)

    # --- Totales de máquina ------------------------------------------------
    PR = P0 / P0_in
    tau = T0 / T0_in
    if PR > 1.0 + 1e-9 and tau > 1.0 + 1e-9:
        eta_poly = ((GAMMA - 1) / GAMMA) * math.log(PR) / math.log(tau)
        eta_isen = (PR ** ((GAMMA - 1) / GAMMA) - 1.0) / (tau - 1.0)
    else:
        eta_poly, eta_isen = 0.05, 0.05
    eta_poly = float(np.clip(eta_poly, 0.05, 1.0))
    eta_isen = float(np.clip(eta_isen, 0.05, 1.0))
    power = mdot * CP * (T0 - T0_in)
    Mu = U_tip1 / a0_in

    # AN² (annulus de salida): límite estructural clásico [m²·rpm²]
    A_exit = frozen_out["areas"][-1] if not off_design else frozen["areas"][-1]
    AN2 = A_exit * RPM ** 2

    st = stage_table
    max_DF_r = max(s["DF_rotor"] for s in st)
    max_DF_s = max(s["DF_stator"] for s in st)
    min_dh_r = min(s["dehaller_rotor"] for s in st)
    min_dh_s = min(s["dehaller_stator"] for s in st)
    min_SM = min(s["SM"] for s in st)
    h_last = st[-1]["h_blade_mm"]

    # ===================== Vector de restricciones g <= 0 =================
    # Magnitud = grado de violación (continua — dominancia de Deb)
    g = np.array([
        max_DF_r / DF_MAX_R - 1.0,                       # difusión rotor
        max_DF_s / DF_MAX_S - 1.0,                       # difusión estátor
        (DEHALLER_MIN - min_dh_r) / DEHALLER_MIN,        # de Haller rotor
        (DEHALLER_MIN - min_dh_s) / DEHALLER_MIN,        # de Haller estátor
        M_rel_tip1 / M_REL_TIP_MAX - 1.0 + g_soft_pen,   # punta transónica+choke
        KOCH_SM_MIN - min_SM,                            # margen de bombeo
        M_exit / M_EXIT_MAX - 1.0,                       # Mach de salida
        (H_BLADE_MIN_MM - h_last) / H_BLADE_MIN_MM,      # manufacturabilidad
    ])

    def _safe(x, default):
        return float(x) if np.isfinite(x) else default

    out = dict(
        PR=_safe(PR, 1.0), eta_poly=_safe(eta_poly, 0.05),
        eta_isen=_safe(eta_isen, 0.05),
        T0_out=_safe(T0, T0_in), T02=_safe(T0, T0_in),
        power_W=_safe(power, 0.0),
        U_tip=_safe(U_tip1, 0.0), U2=_safe(U_tip1, 0.0),  # alias patrón Phy-CC
        Mu=_safe(Mu, 0.0), r_tip_in_mm=_safe(r_tip * 1000, 0.0),
        n_stages=n_st, length_mm=_safe(length_m * 1000, 0.0),
        AN2_m2_rpm2=_safe(AN2, 0.0),
        M_rel_tip1=_safe(M_rel_tip1, 9.9), M_exit=_safe(M_exit, 9.9),
        max_DF_rotor=_safe(max_DF_r, 9.9), max_DF_stator=_safe(max_DF_s, 9.9),
        min_dehaller_rotor=_safe(min_dh_r, 0.0),
        min_dehaller_stator=_safe(min_dh_s, 0.0),
        min_SM=_safe(min_SM, -9.9), h_blade_last_mm=_safe(h_last, 0.0),
        Rx_hub_min=_safe(min(s["Rx_hub"] for s in st), -9.9),
        stage_table=st,
        frozen=frozen_out if not off_design else None,
        g=[_safe(x, 99.0) for x in g],
        feasible=bool(np.all(np.asarray(g) <= 0.0) and np.all(np.isfinite(g))
                      and not choke_any),
        choke_flag=bool(choke_any),
        stall_flag=bool(min_SM < 0.0),
        clipped=clipped,
        source="meanline_L0" + ("_offdesign" if off_design else ""),
    )
    return out


# ==========================================================================
# 2b. Fuera de diseño y mapa de operación (L0)
# ==========================================================================
def offdesign(theta: np.ndarray, rpm: float, mdot: float,
              frozen: dict | None = None) -> dict:
    """Evalúa el diseño theta en otro punto (RPM, ṁ) con la geometría
    congelada del diseño (ángulos de flujo y áreas)."""
    theta = np.asarray(theta, dtype=float)
    if frozen is None:
        frozen = _meanline(theta)["frozen"]
    frozen = dict(frozen)
    frozen.setdefault("phi_d", float(theta[3]))
    return _meanline(theta, rpm=rpm, mdot=mdot, frozen=frozen)


def compressor_map(theta: np.ndarray,
                   speed_fracs=(0.7, 0.8, 0.9, 1.0, 1.05),
                   mdot_range=(0.55, 1.25), n_points=17) -> dict:
    """Mapa de operación L0: speedlines PR(ṁ) con marcadores de límite.

    Surge: SM de Koch <= 0 o rama de pendiente positiva a la izquierda del
    pico de PR. Choke: continuidad sin raíz subsónica en alguna estación.
    Da la TENDENCIA; recalibrar contra CFD/banco (docs/VALIDATION.md).
    """
    theta = np.asarray(theta, dtype=float)
    base = _meanline(theta)
    frozen = dict(base["frozen"])
    frozen.setdefault("phi_d", float(theta[3]))
    RPM_d, mdot_d = float(theta[1]), float(theta[12])
    out = dict(design=dict(RPM=RPM_d, mdot=mdot_d, PR=base["PR"],
                           eta_poly=base["eta_poly"]), speedlines=[])
    for fr in speed_fracs:
        rpm = fr * RPM_d
        pts = []
        for mm in np.linspace(mdot_range[0], mdot_range[1], n_points) * mdot_d:
            r = _meanline(theta, rpm=rpm, mdot=float(mm), frozen=frozen)
            pts.append(dict(
                mdot=float(mm), PR=r["PR"], eta_poly=r["eta_poly"],
                min_SM=r["min_SM"], M_rel_tip1=r["M_rel_tip1"],
                choke=bool(r["choke_flag"]),
                stall=bool(r["stall_flag"])))
        i_peak = int(np.argmax([p["PR"] for p in pts]))
        for j, p in enumerate(pts):
            p["unstable"] = bool(j < i_peak)
        valid = [p for p in pts if not (p["unstable"] or p["choke"] or p["stall"])]
        out["speedlines"].append(dict(
            rpm=float(rpm), speed_frac=float(fr), points=pts,
            mdot_surge=float(valid[0]["mdot"]) if valid else None,
            mdot_choke=float(valid[-1]["mdot"]) if valid else None))
    return out


# ==========================================================================
# 3. L1 — SCM de turbo-design (spool axial multi-fila, modo meanline)
#
# Hallazgos del spike M7 (TD3 1.4.2, documentados en docs/Science.md):
#   * num_streamlines=1 OBLIGATORIO: con >1 el ODE de equilibrio radial
#     (radeq/solve_ivp) colapsa el paso y se cuelga indefinidamente.
#   * El marching de balance de presión NO propaga el remolino de entrada
#     al trabajo de Euler del rotor (el rotor siempre ve entrada axial).
#     Por eso cada etapa se TRANSFORMA conservando el trabajo: entrada
#     axial, rotor con Cu2_eff = ψ·U (β2 equivalente) y estátor que
#     devuelve el flujo a axial. τ queda igual por construcción y TD3
#     aporta la marcha compresible + annulus real.
#   * La pérdida se impone como FixedPressureLoss (Yp referida a la carga
#     dinámica de la fila AGUAS ARRIBA en TD3): los ω̄ del L0 se reescalan
#     a esa referencia. La pérdida politrópica (patrón Phy-CC) NO es
#     alcanzable por el optimizador interno de TD3 en spools axiales.
#   * El bucle "Looping to converge massflow" puede colgarse → cada solve
#     corre en subproceso con timeout PHYAC_L1_TIMEOUT (def. 180 s).
# ==========================================================================
L1_TIMEOUT_S = float(os.environ.get("PHYAC_L1_TIMEOUT", "180"))


def _build_passage(theta: np.ndarray, record: dict):
    """Pasaje meridional AXIAL hub/shroud en metros, con estaciones en
    cada frontera de FILA (el área de salida del rotor importa: sin ella
    TD3 ve el área de entrada de la etapa y sobregira el flujo)."""
    M = _TD_MODS
    st = record["stage_table"]
    frozen = record.get("frozen") or {}
    areas_next = frozen.get("areas", [])
    z, areas = [0.0], [st[0]["A_in_m2"]]
    for i, s in enumerate(st):
        cr = s["chord_rotor_mm"] / 1000.0
        cs = s["chord_stator_mm"] / 1000.0
        z.append(z[-1] + cr)
        areas.append(s["A_rotor_exit_m2"])
        z.append(z[-1] + ROW_GAP_FRACTION * cr + (1 + ROW_GAP_FRACTION) * cs)
        if i + 1 < len(areas_next):
            areas.append(areas_next[i + 1])
        else:
            areas.append(s["A_rotor_exit_m2"])
    # r_mean por estación: la de su etapa (const_mean ⇒ constante)
    r_hub, r_tip = [], []
    r_means = [st[0]["r_mean_mm"] / 1000.0]
    for s in st:
        r_means += [s["r_mean_mm"] / 1000.0] * 2
    for A, r_m in zip(areas, r_means):
        h = A / (2 * math.pi * r_m)
        r_hub.append(r_m - 0.5 * h)
        r_tip.append(r_m + 0.5 * h)
    z = np.array(z)
    return M["Passage"](z, np.array(r_hub), z.copy(), np.array(r_tip),
                        passageType=M["PassageType"].Axial), z


def _scm_solve_direct(theta: np.ndarray) -> dict | None:
    """Corre el spool SCM axial EN ESTE PROCESO; None si diverge.

    Usar _scm_solve (wrapper con timeout) desde código de producto.
    """
    if not _try_load_turbodesign():
        return None
    M = _TD_MODS
    theta = np.asarray(theta, dtype=float)
    T0_in, P0_in, mdot = float(theta[10]), float(theta[11]), float(theta[12])
    RPM = float(theta[1])
    try:
        rec0 = _meanline(theta)
        st = rec0["stage_table"]
        _dbg = os.environ.get("PHYAC_DEBUG_L1")
        _co = contextlib.nullcontext() if _dbg else contextlib.redirect_stdout(io.StringIO())
        _ce = contextlib.nullcontext() if _dbg else contextlib.redirect_stderr(io.StringIO())
        with _co, _ce:
            n = len(st)
            # ángulos equivalentes iniciales (Cx de diseño); el lazo
            # exterior los corrige con el Cx REAL de TD3 (su densidad
            # local difiere de la del L0 y el acople Cx↔trabajo es fuerte
            # a stagger ~45°)
            b2_eff = []
            for s in st:
                U, Cx = s["U_m"], s["Cx"]
                b2_eff.append(math.degrees(math.atan2(U - s["psi"] * U, Cx)))

            PR = tau = None
            for _outer in range(6):
                passage, z_st = _build_passage(theta, rec0)
                z_total = float(z_st[-1])
                inlet = M["Inlet"](hub_location=0, alpha=[0])
                inlet.init_total(P0=P0_in, T0=T0_in, M=0.3)
                outlet = M["Outlet"](num_streamlines=1)
                outlet.init_static(P=0.9 * rec0["PR"] * P0_in,
                                   percent_radii=[0.5])
                rows = []
                for i, s in enumerate(st):
                    U, Cx = s["U_m"], s["Cx"]
                    W2_eff = math.hypot(Cx, U - s["psi"] * U)
                    # ω̄ del L0 reescalados a la referencia de TD3 (q de
                    # la fila aguas arriba): rotor ← q axial, estátor ←
                    # q relativa de salida del rotor.
                    om_r = (s["losses"]["omega_rotor"]
                            + s["losses"]["clearance"]
                            / max(0.5 * s["W1"] ** 2, 1e-3)) \
                        * (s["W1"] / max(Cx, 1e-3)) ** 2
                    om_s = s["losses"]["omega_stator"] * \
                        (s["C2"] / max(W2_eff, 1e-3)) ** 2
                    rows.append(M["make_rotor_row"](
                        hub_location=float(z_st[1 + 2 * i]) / z_total,
                        metal_exit_angle_deg=-b2_eff[i],   # rotores NEG.
                        P0_ratio=max(s["PR"], 1.05),       # gotcha: > 1
                        num_blades=s["n_blades_rotor"],
                        loss_function=M["FixedPressureLoss"](min(om_r, 0.9)),
                    ))
                    rows.append(M["make_stator_row"](
                        hub_location=float(z_st[2 + 2 * i]) / z_total,
                        metal_exit_angle_deg=0.0,          # → axial
                        num_blades=s["n_blades_stator"],
                        loss_function=M["FixedPressureLoss"](min(om_s, 0.9)),
                    ))
                spool = M["CompressorSpool"](
                    passage=passage, massflow=mdot, inlet=inlet,
                    outlet=outlet, rows=rows, rpm=RPM, num_streamlines=1,
                    fluid=M["Solution"]("air.yaml"),
                )
                spool.solve()
                allr = spool._all_rows()
                P0_out = float(np.mean(allr[-1].P0))
                P0_in_s = float(np.mean(allr[0].P0))
                T0_out = float(np.mean(allr[-1].T0))
                PR = P0_out / P0_in_s
                tau = T0_out / T0_in

                # corrección de trabajo por rotor: Cu2_logrado desde ΔT0
                # de la fila; nuevo ángulo con el Cx real implícito
                max_err = 0.0
                T_prev = T0_in
                for i, s in enumerate(st):
                    row_rot = allr[1 + 2 * i]
                    T_rot = float(np.mean(row_rot.T0))
                    dT = max(T_rot - T_prev, 1e-3)
                    T_prev = float(np.mean(allr[2 + 2 * i].T0))
                    U = s["U_m"]
                    Cu2_ach = dT * CP / U
                    Cu2_tgt = s["psi"] * U
                    tan_b2 = math.tan(math.radians(b2_eff[i]))
                    Cx_ach = max((U - Cu2_ach) / max(tan_b2, 1e-3), 10.0)
                    b2_new = math.degrees(math.atan2(U - Cu2_tgt, Cx_ach))
                    # amortiguado: el acople trabajo↔densidad↔Cx hace
                    # divergir el punto fijo sin relajación
                    b2_eff[i] = 0.5 * (b2_eff[i] + b2_new)
                    max_err = max(max_err,
                                  abs(Cu2_ach - Cu2_tgt) / max(Cu2_tgt, 1e-3))
                if max_err < 0.02:
                    break
        PR_L0 = rec0["PR"]
        if not (np.isfinite(PR) and 0.70 * PR_L0 <= PR <= 1.40 * PR_L0
                and 1.02 < PR < 40.0 and tau > 1.0):
            if os.environ.get("PHYAC_DEBUG_L1"):
                print(f"[L1] rechazado por ventana: PR={PR} (L0 {PR_L0}), "
                      f"tau={tau}")
            return None
        eta = float(np.clip(((GAMMA - 1) / GAMMA) * math.log(PR)
                            / max(math.log(tau), 1e-9), 0.05, 0.99))
        return dict(PR=PR, eta_poly=eta, T0_out=T0_out, source="scm_L1")
    except Exception:
        if os.environ.get("PHYAC_DEBUG_L1"):
            import traceback
            traceback.print_exc()
        return None


def _scm_worker_entry(theta_list, q):
    """Entry point del subproceso de _scm_solve (picklable, spawn)."""
    try:
        q.put(_scm_solve_direct(np.array(theta_list, dtype=float)))
    except Exception:
        q.put(None)


def _scm_solve(theta: np.ndarray) -> dict | None:
    """Spool SCM axial con timeout duro (subproceso).

    El bucle de convergencia de gasto de TD3 puede colgarse (spike M7);
    el subproceso se termina a los L1_TIMEOUT_S segundos y el punto
    degrada a L0 etiquetado.
    """
    if not _try_load_turbodesign():
        return None
    import multiprocessing as _mp
    ctx = _mp.get_context("spawn")
    q = ctx.Queue()
    p = ctx.Process(target=_scm_worker_entry,
                    args=(np.asarray(theta, dtype=float).tolist(), q))
    p.start()
    p.join(L1_TIMEOUT_S)
    if p.is_alive():
        p.terminate()
        p.join()
        return None
    try:
        return q.get_nowait()
    except Exception:
        return None


# ==========================================================================
# 4. Calibración de alta fidelidad (hook L2)
# ==========================================================================
class HiFiCalibration:
    """Corrección afín por salida: y_cal = a*y + b, ajustada con pares
    (y_modelo, y_hifi). Con <2 pares aplica identidad."""

    KEYS = ("PR", "eta_poly")

    def __init__(self):
        self.pairs: dict[str, list[tuple[float, float]]] = {k: [] for k in self.KEYS}
        self.coef: dict[str, tuple[float, float]] = {k: (1.0, 0.0) for k in self.KEYS}

    def register(self, y_model: dict, y_hifi: dict):
        for k in self.KEYS:
            if k in y_model and k in y_hifi:
                self.pairs[k].append((float(y_model[k]), float(y_hifi[k])))
        self._refit()

    def _refit(self):
        for k, p in self.pairs.items():
            if len(p) >= 2:
                x = np.array([q[0] for q in p])
                y = np.array([q[1] for q in p])
                a, b = np.polyfit(x, y, 1)
                self.coef[k] = (float(a), float(b))

    def apply(self, rec: dict) -> dict:
        rec = dict(rec)
        for k, (a, b) in self.coef.items():
            if k in rec and (a, b) != (1.0, 0.0):
                rec[k + "_raw"] = rec[k]
                rec[k] = a * rec[k] + b
                rec["calibrated"] = True
        return rec


CALIBRATION = HiFiCalibration()


def register_hifi_pair(theta: np.ndarray, y_hifi: dict, fidelity=Fidelity.L0):
    """Reinyecta un resultado de CFD/banco al loop (estilo Noyron)."""
    base = evaluate(theta, fidelity=fidelity, use_cache=True, calibrate=False)
    CALIBRATION.register(base, y_hifi)


# ==========================================================================
# 5. API pública: evaluate / evaluate_batch / caché / features
# ==========================================================================
_CACHE: dict[str, dict] = {}
CACHE_PATH = os.environ.get("PHYAC_CACHE", "")


def set_cache_path(path: str):
    """Activa (o cambia) el caché persistente y carga lo ya existente."""
    global CACHE_PATH
    CACHE_PATH = path or ""
    _cache_load()


def _cache_load():
    if CACHE_PATH and os.path.exists(CACHE_PATH):
        with open(CACHE_PATH, encoding="utf-8") as f:
            for line in f:
                try:
                    rec = json.loads(line)
                    _CACHE[rec.pop("_key")] = rec
                except Exception:
                    pass


def _cache_store(key: str, rec: dict):
    _CACHE[key] = rec
    # A disco solo lo caro (intentos L1+); el L0 se recomputa en ~1 ms.
    if CACHE_PATH and rec.get("fidelity", 0) >= int(Fidelity.L1):
        slim = {k: v for k, v in rec.items() if k != "frozen"}
        with open(CACHE_PATH, "a", encoding="utf-8") as f:
            f.write(json.dumps({"_key": key, **slim}) + "\n")


_cache_load()


def evaluate(theta: np.ndarray, fidelity: Fidelity = Fidelity.L1,
             use_cache: bool = True, calibrate: bool = True) -> dict:
    """Evalúa un punto de diseño físico (14-D, unidades físicas).

    L1 corre el SCM axial y, si diverge o no está instalado, degrada a L0
    etiquetando la fuente. Siempre devuelve el record completo del
    meanline (stage_table, restricciones, banderas) — el SCM solo
    sobreescribe PR/eta cuando converge.
    """
    theta = np.asarray(theta, dtype=float)
    key = _theta_key(theta, int(fidelity))
    if use_cache and key in _CACHE:
        rec = dict(_CACHE[key])
    else:
        rec = _meanline(theta)
        rec.update({n: float(v) for n, v in zip(VAR_NAMES, theta)})
        rec["n_stages"] = int(round(theta[0]))
        rec["fidelity"] = int(fidelity)
        if fidelity >= Fidelity.L1 and rec["feasible"]:
            scm = _scm_solve(theta)
            if scm:
                rec.update(scm)
            else:
                rec["source"] = rec["source"] + "(L1_unavailable_or_diverged)"
        if use_cache:
            _cache_store(key, rec)
    if calibrate:
        rec = CALIBRATION.apply(rec)
    return rec


def evaluate_batch(thetas: np.ndarray, fidelity: Fidelity = Fidelity.L1,
                   n_workers: int = 1) -> list[dict]:
    """Evalúa un batch. n_workers>1 paraleliza por procesos (útil en L1)."""
    thetas = np.atleast_2d(np.asarray(thetas, dtype=float))
    if n_workers <= 1 or fidelity == Fidelity.L0:
        return [evaluate(t, fidelity) for t in thetas]
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        return list(ex.map(_eval_worker, [(t.tolist(), int(fidelity))
                                          for t in thetas]))


def _eval_worker(args):
    theta, fid = args
    return evaluate(np.array(theta), Fidelity(fid))


# --- Physics embedding para el surrogate (capa 2) -------------------------
FEATURE_NAMES = VAR_NAMES + [
    "f_Utip", "f_Mu", "f_DFr", "f_DFs", "f_Mrel1", "f_SM",
    "f_HTRout", "f_hlast", "f_Mexit", "f_RxHub", "L0_logPR", "L0_eta",
]


def physics_features(theta: np.ndarray, record: dict | None = None) -> np.ndarray:
    """Vector de entrada del surrogate: theta normalizado + features
    físicos cerrados + predicciones L0 (residual learning).

    `record`: meanline L0 ya calculado para este mismo theta (DEBE ser de
    fidelidad L0 — un record L1 corrompería el residual learning).
    """
    theta = np.asarray(theta, dtype=float)
    rec = record if record is not None else _meanline(theta)
    st_last = rec["stage_table"][-1]
    htr_out = st_last["r_hub_mm"] / max(st_last["r_tip_mm"], 1e-6)
    feats = np.array([
        rec["U_tip"] / 550.0,
        rec["Mu"] / 1.8,
        min(rec["max_DF_rotor"], 1.2) / 1.2,
        min(rec["max_DF_stator"], 1.2) / 1.2,
        min(rec["M_rel_tip1"], 2.0) / 2.0,
        float(np.clip(rec["min_SM"], -1.0, 1.0)),
        htr_out,
        min(rec["h_blade_last_mm"], 120.0) / 120.0,
        min(rec["M_exit"], 1.2) / 1.2,
        float(np.clip(rec["Rx_hub_min"], -1.0, 1.5)) / 1.5,
        math.log(max(rec["PR"], 1.001)) / math.log(30.0),
        rec["eta_poly"],
    ])
    return np.concatenate([normalize(theta), feats])


N_FEATURES = NDIM + 12


if __name__ == "__main__":
    # Sanity check: compresor axial de 4 etapas representativo
    theta_ref = np.array([4.0, 12_500.0, 0.62, 0.55, 0.32, -0.10,
                          0.60, 1.20, 1.10, 2.20, 288.15, 101_325.0, 25.0])
    r = evaluate(theta_ref, fidelity=Fidelity.L1)
    print("Diseño axial de referencia (4 etapas):")
    for k in ["PR", "eta_poly", "eta_isen", "U_tip", "Mu", "r_tip_in_mm",
              "M_rel_tip1", "max_DF_rotor", "min_SM", "length_mm", "power_W",
              "feasible", "g", "source"]:
        print(f"  {k:16s} = {r[k]}")
