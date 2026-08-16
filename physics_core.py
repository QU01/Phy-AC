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
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from enum import IntEnum

import numpy as np

warnings.filterwarnings("ignore")

# ==========================================================================
# 0. Espacio de diseño (15-D) — distribuciones escalares + pendientes por
#    etapa (fase 8; ver docs/PLAN_Fases_7_8.md §8.1 y docs/Science.md §1.1)
# ==========================================================================
DESIGN_VARS = [
    # nombre          min        max          unidad
    ("n_stages",      1.0,       8.0),        # -    (se redondea a entero)
    ("RPM",       5_000.0,  25_000.0),        # rev/min
    ("HTR_in",        0.40,      0.80),       # -    hub-to-tip ratio entrada
    ("phi1",          0.35,      0.80),       # -    coef. de flujo Cx/Um (medio)
    ("psi_mid",       0.22,      0.45),       # -    carga media Δh0/Um²
    ("psi_slope",    -0.30,      0.30),       # -    gradiente de ψ frontal→trasero
    ("Rx_mean",       0.50,      0.85),       # -    grado de reacción medio
    ("sigma_r",       0.90,      1.60),       # -    solidez rotor (media)
    ("sigma_s",       0.80,      1.50),       # -    solidez estátor (media)
    ("AR",            1.20,      3.50),       # -    aspect ratio rotor h/c
    ("T0_in",       270.0,     320.0),        # K
    ("P0_in",    85_000.0, 110_000.0),        # Pa
    ("massflow",      5.0,      60.0),        # kg/s
    # --- fase 8: pendientes lineales por etapa (AL FINAL: los índices
    #     10-12 del punto de operación no se mueven — compatibilidad de
    #     checkpoints, fix_operating_point y CLI). φ cae hacia atrás y Rx
    #     sube en las máquinas reales (E³): φ_i = φ1·(1+s_φ·ξ_i),
    #     Rx_i = Rx·(1+s_Rx·ξ_i), ξ_i = 2i/(N−1)−1 (esquema de psi_slope).
    ("phi_slope",    -0.25,      0.25),       # -    gradiente de φ frontal→trasero
    ("Rx_slope",     -0.20,      0.20),       # -    gradiente de Rx frontal→trasero
]
# NOTA: r_tip de entrada NO es variable de diseño — queda DERIVADO de la
# continuidad (dados φ1, HTR, RPM, ṁ el annulus de entrada es único). Con
# r_tip libre el espacio quedaba sobredeterminado: la continuidad imponía
# un φ distinto del pedido y toda la parametrización de Smith se rompía.
VAR_NAMES = [v[0] for v in DESIGN_VARS]
BOUNDS_LO = np.array([v[1] for v in DESIGN_VARS])
BOUNDS_HI = np.array([v[2] for v in DESIGN_VARS])
NDIM = len(DESIGN_VARS)
NDIM_LEGACY = 13          # θ pre-fase-8 (sin pendientes) — ver pad_theta

GAMMA, CP, RGAS = 1.4, 1005.0, 287.0
#   GAMMA/CP se conservan como REFERENCIA (288 K) para el código que aún
#   los usa (geometría, estructura, tests). El stage-stacking usa desde la
#   fase 9 el gas caloríficamente imperfecto de §0.1 — a τ ≈ 2.4 (E³) el
#   cp del aire sube ~7% y γ cae a ~1.36; con cp constante el error entra
#   DIRECTO en PR (Dixon & Hall §5.9) y en cada Mach del stacking.
MU_AIR = 1.81e-5          # Pa·s
RE_REF = 1.0e6            # Re de cuerda nominal de las correlaciones de
#                           pérdida (Koch & Smith 1976: perfil liso a 10⁶)
RE_LAM = 2.0e5            # frontera de separación laminar (Wassell 1968)

# --- Parámetros del meanline (calibrables; ver docs/VALIDATION.md) --------
TIP_CLEARANCE_MM = float(os.environ.get("PHYAC_TIP_CLEARANCE_MM", "0.4"))
#   ε holgura de punta ABSOLUTA del rotor [mm]. En máquinas reales ε la
#   fijan tolerancias mecánicas y crecimiento térmico, NO la altura del
#   álabe — por eso ε/h CRECE hacia las etapas traseras (h cae con la
#   compresión), el efecto físico que un ε/h fijo borraba. 0.4 mm es
#   típico de compresor pequeño; la validación inyecta el ε publicado de
#   cada máquina (validation/machines.py) y el contrato de geometría
#   emite ESTE MISMO valor (annulus.tip_clearance_mm) — física y pieza
#   impresa comparten la holgura.
K_TIP_CLEARANCE = 1.8         # pendiente Δη/Δ(ε/h) del tramo lineal
EPS_H_KNEE = 0.008            # rodilla inferior: óptimo de holgura
EPS_H_HIGH = 0.034            # fin del tramo lineal (punta descargada)
EPS_H_MAX = 0.08              # clamp de ε/h para θ degenerado (g continuo)
BLOCKAGE_INIT = 0.98          # KB de la primera etapa
BLOCKAGE_PER_STAGE = 0.005    # OBSOLETO desde la fase 9: el bloqueo lo da
#                               ahora endwall_blockage() (Koch & Smith
#                               1976), que lo hace CRECER con la carga y
#                               con la holgura en vez de linealmente.
NU_OVER_DELTA = 0.48          # Σν/Σδ* — espesor de fuerza tangencial sobre
#                               espesor de desplazamiento (Koch & Smith
#                               1976, Fig. 10: línea única ≈ 0.48). Entra
#                               en el débito de rendimiento de pared, ec. 2.
ANNULUS_MODE = "const_mean"   # const_hub | const_mean | const_tip
VORTEX_N = float(os.environ.get("PHYAC_VORTEX_N", "-1.0"))
#   Exponente del torbellino: Cu(r) = Cu_m·(r/r_m)^n.
#     n = −1  vórtice LIBRE (r·Cu = cte) — el clásico y el default; único
#             con Cx radialmente constante.
#     n =  0  "constante": Cu igual en todo el span.
#     n = +1  vórtice FORZADO (sólido rígido) — reacción constante.
#   La práctica real de fans y de etapas frontales de HPC vive entre −1 y
#   0 ("controlled vortex"): el vórtice libre dispara Cu en el cubo y la
#   reacción de raíz se hunde, que es lo que hace INVIABLE un fan de HTR
#   bajo (Dixon & Hall cap. 6). Con n = −1 todo colapsa BIT A BIT al
#   comportamiento previo a esta fase.
#   Equilibrio radial simple (SRE), con la forma cerrada para Cu ∝ r^n:
#       Cx²(r) = Cx_m² − Cu_m²·((n+1)/n)·[(r/r_m)^{2n} − 1]
#   que en n = −1 da (n+1)/n = 0 ⇒ Cx constante, y en n = +1 reproduce el
#   resultado de cuerpo sólido.


def vortex_cu(cu_mean: float, r: float, r_mean: float,
              n: float = None) -> float:
    """Cu(r) del torbellino de exponente n."""
    n = VORTEX_N if n is None else n
    return cu_mean * (max(r, 1e-9) / max(r_mean, 1e-9)) ** n


def vortex_cx(cx_mean: float, cu_mean: float, r: float, r_mean: float,
              n: float = None) -> float:
    """Cx(r) por equilibrio radial simple para Cu ∝ r^n."""
    n = VORTEX_N if n is None else n
    if abs(n + 1.0) < 1e-12:
        return cx_mean                      # vórtice libre: Cx constante
    rr = max(r, 1e-9) / max(r_mean, 1e-9)
    if abs(n) < 1e-6:
        # límite n→0 (Cu constante en el span): (n+1)/n·[(r/r_m)^{2n} − 1]
        # → 2·ln(r/r_m). Rama continua, no un caso especial disfrazado.
        c2 = cx_mean ** 2 - 2.0 * cu_mean ** 2 * math.log(rr)
    else:
        c2 = cx_mean ** 2 - cu_mean ** 2 * ((n + 1.0) / n) * (rr ** (2.0 * n) - 1.0)
    return math.sqrt(max(c2, (0.25 * cx_mean) ** 2))


def vortex_reaction(rx_mean: float, r: float, r_mean: float,
                    n: float = None) -> float:
    """Rx(r) = 1 − (1−Rx_m)·(r/r_m)^{n−1}. En n = −1 recupera la forma
    de vórtice libre 1 − (1−Rx_m)·(r_m/r)²; en n = +1 es constante."""
    n = VORTEX_N if n is None else n
    return 1.0 - (1.0 - rx_mean) * \
        (max(r, 1e-9) / max(r_mean, 1e-9)) ** (n - 1.0)
STATOR_AR_FACTOR = 1.10       # AR_estator = 1.10 · AR_rotor
ROW_GAP_FRACTION = 0.25       # gap axial inter-fila / cuerda
CHORD_TAPER = 0.85            # c_tip/c_hub (usado por la capa 5a)
K_PROFILE = 1.24              # multiplicador global de pérdida de perfil.
#                               El ajuste θ/c de Lieblein es de CASCADA 2-D
#                               a incidencia de diseño y superficie lisa; en
#                               máquina real la pérdida de perfil es mayor
#                               (acabado, contracción del tubo de corriente,
#                               no estacionariedad). 1.24 es el factor global
#                               que, con el endwall de Koch & Smith ya
#                               físico, deja las cuatro máquinas NASA en
#                               tolerancia (docs/VALIDATION.md §3).
K_ENDWALL = 1.0               # multiplicador de pérdida secundaria/endwall.
#                               Recalibrado 2026-07-17 junto con el modelo
#                               de holgura por fila: el débito uniforme
#                               viejo (ε/h=1.5% fijo ≈ 3 pts) absorbía en
#                               silencio el endwall que el CDa simple de
#                               Howell subestima (Koch & Smith 1976
#                               atribuyen al endwall bastante más). Con la
#                               ε publicada de cada máquina, K=1.4 deja
#                               las 4 máquinas NASA en tolerancia
#                               (validation/RESULTS.md) sin ajuste por
#                               máquina individual.
K_MACH_PROFILE = 0.50         # pendiente del factor de Mach sobre la
#                               pérdida de perfil (Koch & Smith 1976,
#                               Fig. 6): f_M = 1 + K·(M²−0.0625).
K_SHOCK = 1.00                # multiplicador de pérdida de choque. El
#                               choque real del pasaje es OBLICUO; el modelo
#                               de choque normal al Mach de entrada
#                               sobreestima. Calibrado contra NASA Rotor
#                               37/67 (validation/RESULTS.md) manteniendo
#                               Stage 35 dentro de tolerancia.
MX_CHOKE = 0.78               # Mach axial máx. antes de declarar choke local

# --- Off-design (2026-07-17): bucket W(M), desviación y VSV ---------------
W_INC_LOW_DEG = 10.0    # semiancho del bucket de incidencia a M_rel <= 0.2
#                         (comportamiento clásico de baja velocidad)
W_INC_HIGH_DEG = 3.5    # semiancho a M_rel >= 0.8 — el rango de incidencia
#                         útil se ESTRECHA con el Mach (Aungier 2003,
#                         off-design performance; a M≈0.8 el bucket real es
#                         ±2-4°, no ±10°)
K_DEV_OFFDESIGN = 0.30  # dδ/di con incidencia POSITIVA (sub-giro creciente
#                         hacia stall — Creveling 1968, NASA CR-72427;
#                         práctica SP-36). En i<=0 no aplica.
VSV_GAIN_DEG = 50.0     # cierre de VSV por unidad de defecto de velocidad
#                         (schedule "auto" del mapa: Δ = GAIN·(1−N/Nd))
VSV_MAX_DEG = 35.0      # tope físico del restagger de estátores frontales

# Work-done factor de Howell (mediado por nº de etapa; tabla clásica)
_HOWELL_WDF = [0.982, 0.952, 0.929, 0.910, 0.895, 0.883, 0.875, 0.868]

# Límites de las restricciones aero
DF_MAX_R = 0.55           # factor de difusión de Lieblein, rotor
DF_MAX_S = 0.60           # factor de difusión, estátor
DEHALLER_MIN = 0.72       # W2/W1 (rotor) y C3/C2 (estátor) mínimos
M_REL_TIP_MAX = 1.35      # Mach relativo en punta del primer rotor
CH_STALL_MAX = 0.48       # OBSOLETO desde la fase 9: la capacidad de
#                           subida de presión al stall ya no es una
#                           constante sino la correlación de Koch 1981
#                           (koch_ch_stall) — depende de L/g₂, Re, ε/g,
#                           Δz/s y del tipo de triángulo (𝔉_ef). Se
#                           conserva el símbolo por compatibilidad.
KOCH_SM_MIN = 0.10        # margen mínimo de Koch en el PUNTO DE DISEÑO
M_EXIT_MAX = 0.55         # Mach absoluto de salida (amigable al difusor)
H_BLADE_MIN_MM = 8.0      # altura mínima del álabe de la última etapa
RX_HUB_MIN = 0.05         # reacción mínima en la raíz (fase 9). Con
#                           vórtice libre y HTR bajo la reacción de raíz
#                           se hunde y el rotor ACELERA el flujo en el
#                           cubo; por debajo de ~0.05 el diseño no es
#                           construible (Dixon & Hall cap. 6).
SM_FLOW_MIN = 0.15        # margen de bombeo MÍNIMO en gasto sobre la
#                           línea de trabajo (práctica: 15-25%). Es la
#                           forma en que la industria especifica el
#                           margen — ver surge_margin().

N_CONSTRAINTS = 9         # dim del vector g(theta) — solo aerodinámica


class Fidelity(IntEnum):
    L0 = 0   # stage-stacking meanline analítico
    L1 = 1   # SCM turbo-design (si está disponible)


# ==========================================================================
# 0.1 Gas caloríficamente imperfecto (fase 9)
# ==========================================================================
# cp(T) del aire seco por ajuste cuadrático de las tablas JANAF/NASA en
# 250-1000 K (<1% en todo el rango; anclado en 300/600/900 K):
#
#     cp(T) = A + B·T + C·T²      [J/(kg·K)]
#
# De ahí salen las tres funciones de estado que el stacking necesita:
#
#     h(T)   = ∫cp dT                         entalpía
#     phi(T) = ∫(cp/T) dT                     "función de entropía"
#     s2−s1  = phi(T2) − phi(T1) − R·ln(P2/P1)
#
# Con phi(T) los rendimientos dejan de necesitar γ constante y quedan
# EXACTOS para gas imperfecto (Dixon & Hall §1.11 generalizado):
#
#     η_poly = R·ln(PR) / [phi(T2) − phi(T1)]
#     η_isen = [h(T2s) − h(T1)] / [h(T2) − h(T1)],
#              con phi(T2s) = phi(T1) + R·ln(PR)
#
# γ(T) = cp/(cp−R) se usa localmente para Mach y relaciones isentrópicas
# de estación (donde la variación de T es pequeña y el gas perfecto local
# es exacto a segundo orden).
CP_A, CP_B, CP_C = 983.0, 0.0333333, 1.3333e-4
GAS_VARIABLE = os.environ.get("PHYAC_GAS", "variable").lower() != "const"
#   PHYAC_GAS=const restaura el gas perfecto (cp=1005, γ=1.4) — solo para
#   reproducir resultados pre-fase-9; NO es el modo de producto.


def cp_air(T: float) -> float:
    """cp(T) del aire [J/(kg·K)]. Clampeado a 200-1400 K."""
    if not GAS_VARIABLE:
        return CP
    t = min(max(float(T), 200.0), 1400.0)
    return CP_A + CP_B * t + CP_C * t * t


def gamma_air(T: float) -> float:
    """γ(T) = cp/(cp−R)."""
    if not GAS_VARIABLE:
        return GAMMA
    c = cp_air(T)
    return c / (c - RGAS)


def h_air(T: float) -> float:
    """Entalpía específica [J/kg] (referencia arbitraria, solo diferencias)."""
    if not GAS_VARIABLE:
        return CP * float(T)
    t = min(max(float(T), 200.0), 1400.0)
    return CP_A * t + 0.5 * CP_B * t * t + (CP_C / 3.0) * t ** 3


def phi_air(T: float) -> float:
    """Función de entropía ∫cp/T dT [J/(kg·K)]."""
    t = min(max(float(T), 200.0), 1400.0)
    if not GAS_VARIABLE:
        return CP * math.log(t)
    return CP_A * math.log(t) + CP_B * t + 0.5 * CP_C * t * t


def T_from_h(h: float, T_guess: float = 500.0) -> float:
    """Inversa de h_air por Newton (cp = dh/dT). 3 iteraciones bastan."""
    if not GAS_VARIABLE:
        return float(h) / CP
    T = min(max(float(T_guess), 200.0), 1400.0)
    for _ in range(6):
        f = h_air(T) - h
        d = cp_air(T)
        T_new = min(max(T - f / d, 200.0), 1400.0)
        if abs(T_new - T) < 1e-9:
            return T_new
        T = T_new
    return T


def T_from_phi(p: float, T_guess: float = 500.0) -> float:
    """Inversa de phi_air por Newton (dphi/dT = cp/T)."""
    T = min(max(float(T_guess), 200.0), 1400.0)
    for _ in range(8):
        f = phi_air(T) - p
        d = cp_air(T) / T
        T_new = min(max(T - f / d, 200.0), 1400.0)
        if abs(T_new - T) < 1e-9:
            return T_new
        T = T_new
    return T


def pr_from_work(T0_in: float, dh_ideal: float) -> float:
    """P02/P01 de un proceso con salto de entalpía ISENTRÓPICO dh_ideal
    desde T0_in (gas imperfecto). Equivale a (1+η·Δh0/(cp·T0))^(γ/(γ−1))
    en el límite de cp constante."""
    if dh_ideal <= 0.0:
        return 1.0
    T2s = T_from_h(h_air(T0_in) + dh_ideal, T0_in + dh_ideal / 1050.0)
    return math.exp((phi_air(T2s) - phi_air(T0_in)) / RGAS)


def pad_theta(theta: np.ndarray) -> np.ndarray:
    """θ legacy de 13 (pre-fase-8, sin pendientes) → 15-D con
    phi_slope = Rx_slope = 0.0.

    Con pendientes cero las distribuciones por etapa colapsan EXACTAMENTE
    a los escalares (φ_i = φ1·1.0, bit-exacto), así que un θ viejo produce
    el mismo record que antes de la fase 8 — las anclas de regresión y los
    checkpoints antiguos siguen siendo válidos sin re-congelar nada.
    """
    theta = np.asarray(theta, dtype=float).ravel()
    if theta.shape[0] == NDIM:
        return theta
    if theta.shape[0] == NDIM_LEGACY:
        return np.concatenate([theta, [0.0, 0.0]])
    raise ValueError(
        f"theta debe tener {NDIM_LEGACY} (legacy) o {NDIM} componentes — "
        f"llegaron {theta.shape[0]}")


def denormalize(u: np.ndarray) -> np.ndarray:
    """[0,1]^15 -> unidades físicas."""
    return BOUNDS_LO + np.asarray(u, dtype=float) * (BOUNDS_HI - BOUNDS_LO)


def normalize(theta: np.ndarray) -> np.ndarray:
    return (pad_theta(theta) - BOUNDS_LO) / (BOUNDS_HI - BOUNDS_LO)


def _physics_params_key() -> str:
    """Huella de los parámetros de MÓDULO que cambian el resultado del
    meanline sin cambiar θ. Sin esto el caché devuelve el record de otra
    configuración: la validación ya lo esquivaba con `use_cache=False` al
    inyectar la ε publicada de cada máquina, y el exponente de torbellino
    lo habría reproducido en silencio."""
    return f"{TIP_CLEARANCE_MM:.6g}|{VORTEX_N:.6g}|{int(GAS_VARIABLE)}"


def _theta_key(theta: np.ndarray, fidelity: int, extra: str = "") -> str:
    arr = np.round(np.asarray(theta, dtype=float), 9)
    return hashlib.sha1(arr.tobytes() + bytes([fidelity])
                        + extra.encode()
                        + _physics_params_key().encode()).hexdigest()


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
def _normal_shock_p0_ratio(M: float, g: float = GAMMA) -> float:
    """P02/P01 a través de un choque normal (M>1)."""
    if M <= 1.0:
        return 1.0
    t1 = ((g + 1) * M * M / ((g - 1) * M * M + 2)) ** (g / (g - 1))
    t2 = ((g + 1) / (2 * g * M * M - (g - 1))) ** (1 / (g - 1))
    return t1 * t2


def _re_correction(re_c: float | None) -> float:
    """Multiplicador de pérdida viscosa por número de Reynolds de cuerda.

    Las correlaciones de perfil/endwall son nominales a Re_c ≈ 10⁶ con
    superficie lisa (Koch & Smith 1976). Por debajo la pérdida crece como
    Re^-0.2 (capa límite turbulenta lisa; en máquina completa Wassell 1968
    y Schäffler 1980 miden (1−η) ∝ Re^-n con n ≈ 0.1–0.2) y como Re^-0.5
    bajo RE_LAM (separación laminar). Continua en ambas fronteras. SIN
    crédito por encima de RE_REF (conservador: mantiene la calibración
    NASA, cuyas máquinas corren a Re ≳ 10⁶). Se aplica a perfil y endwall
    (pérdidas de fricción), no al choque ni a la holgura de punta.
    """
    if re_c is None or not (re_c > 0.0):
        return 1.0
    if re_c >= RE_REF:
        return 1.0
    if re_c >= RE_LAM:
        return (RE_REF / re_c) ** 0.2
    return (RE_REF / RE_LAM) ** 0.2 * (RE_LAM / re_c) ** 0.5


def _clearance_loss_frac(eps_over_h: float) -> float:
    """Fracción del trabajo de etapa Δh0 perdida por fuga de punta (rotor).

    Mecanismo: el gasto de fuga cruza el hueco impulsado por la carga del
    álabe y se mezcla desalineado con el pasaje (Denton 1993, ζ_TL;
    Storer & Cumpsty 1991/94: pérdida ∝ gasto de fuga × desalineación).
    Regímenes medidos por Sakulkaew, Tan et al. (2013, J. Turbomach.):
    bajo ε/h ≈ 0.8% la sensibilidad cae (cizalladura de carcasa compite
    con la fuga — existe un óptimo de holgura), tramo LINEAL en
    0.8–3.4% con ~1.6 pts de η por 1% de ε/h, y sensibilidad decreciente
    por encima (la punta se descarga). Se modela como rampa a tramos con
    pendiente K_TIP_CLEARANCE en el tramo lineal y ½·K fuera de él —
    CONTINUA y monótona (la dominancia de Deb exige g continuo), con
    ε/h clampeado a EPS_H_MAX para θ degenerados.
    """
    k = K_TIP_CLEARANCE
    e = min(max(eps_over_h, 0.0), EPS_H_MAX)
    if e <= EPS_H_KNEE:
        return 0.5 * k * e
    f = 0.5 * k * EPS_H_KNEE
    if e <= EPS_H_HIGH:
        return f + k * (e - EPS_H_KNEE)
    f += k * (EPS_H_HIGH - EPS_H_KNEE)
    return f + 0.5 * k * (e - EPS_H_HIGH)


def _lieblein_theta_c(deq: float) -> float:
    """Espesor de momento de estela θ/c vs difusión equivalente (Lieblein
    1959, fit clásico). Diverge al acercarse Deq→e^(1/0.95): se acota."""
    deq = min(max(deq, 1.0), 2.6)
    denom = 1.0 - 0.95 * math.log(deq)
    return 0.0045 / max(denom, 0.08)


def _shock_omega(M: float, g: float = GAMMA) -> float:
    """ω̄ de choque normal referido a la presión dinámica de entrada."""
    if M <= 1.0:
        return 0.0
    dyn = 1.0 - (1.0 + 0.5 * (g - 1) * M ** 2) ** (-g / (g - 1))
    return (1.0 - _normal_shock_p0_ratio(M, g)) / max(dyn, 1e-3)


def _row_losses(beta1: float, beta2: float, W1: float, W2: float,
                sigma: float, h_over_c: float, M_tip: float,
                M_mean: float, re_c: float | None = None,
                gam: float = GAMMA) -> tuple[float, float, dict]:
    """Pérdidas de una fila en el marco relativo a la fila.

    beta1/beta2 en rad (ángulos de flujo entrada/salida respecto al eje,
    positivos), W1/W2 velocidades entrada/salida del marco de la fila,
    M_tip/M_mean Mach de entrada en punta y en la línea media (el choque
    se promedia entre ambos, práctica de Miller — el choque solo cubre el
    span exterior), re_c Reynolds de cuerda de la fila (None = sin
    corrección, equivale a Re ≥ 10⁶). Devuelve (omega_bar_total,
    dh_loss [J/kg], desglose). ω̄ está referida a la presión dinámica de
    ENTRADA de la fila.
    """
    cb1, cb2 = math.cos(beta1), math.cos(beta2)
    tb1, tb2 = math.tan(beta1), math.tan(beta2)
    f_re = _re_correction(re_c)
    # Difusión equivalente de Lieblein (forma de circulación, diseño)
    deq = (cb2 / max(cb1, 1e-3)) * (
        1.12 + 0.61 * (cb1 ** 2 / max(sigma, 1e-3)) * abs(tb1 - tb2))
    theta_c = _lieblein_theta_c(deq)
    # Efecto del Mach de entrada sobre la pérdida de PERFIL (Koch & Smith
    # 1976, §"Blade Profile Loss" y Fig. 6). La correlación de Lieblein es
    # incompresible: no sabe que al subir M el espesor de momento baja
    # pero el FACTOR DE FORMA sube, y que el balance neto casi DUPLICA el
    # coeficiente de pérdida entre M 0.1 y 1.5 (Fig. 6: ω̄ de rotor
    # 0.020 → 0.038). Ajuste de esa figura, anclado en M ≈ 0.25:
    f_mach = 1.0 + K_MACH_PROFILE * max(M_mean ** 2 - 0.0625, 0.0)
    # ω̄ de perfil (Lieblein): 2·(θ/c)·(σ/cosβ2)·(cosβ1/cosβ2)²
    om_profile = K_PROFILE * f_re * f_mach * 2.0 * theta_c * \
        (sigma / max(cb2, 1e-3)) * (cb1 / max(cb2, 1e-3)) ** 2
    # Secundarias (vórtice de pasaje, Howell): CDs = 0.018·CL². El término
    # de ANNULUS de Howell (CDa = 0.020·s/h) desaparece desde la fase 9: la
    # capa límite de pared la trata ahora el modelo de Koch & Smith 1976 a
    # nivel de ETAPA (espesores de desplazamiento y de fuerza tangencial,
    # ec. 2-3 del paper), que es física real en vez de un arrastre
    # equivalente — y por eso K_ENDWALL vuelve a 1.0 (el 1.4 existía para
    # tapar justo lo que el CDa subestimaba).
    beta_m = math.atan(0.5 * (tb1 + tb2))
    cbm = math.cos(beta_m)
    CL = 2.0 * (1.0 / max(sigma, 1e-3)) * cbm * abs(tb1 - tb2)
    cd_sec = 0.018 * CL ** 2
    om_endwall = K_ENDWALL * f_re * cd_sec * sigma * cb1 ** 2 / \
        max(cbm, 1e-3) ** 3
    # Choque: promedio punta/media (dos puntos del span, Miller-style);
    # sin f_re — la pérdida de choque no es de fricción
    om_shock = K_SHOCK * 0.5 * (_shock_omega(M_tip, gam)
                                + _shock_omega(M_mean, gam))
    om_total = om_profile + om_endwall + om_shock
    # ΔP0rel = ω̄·½ρW1² → pérdida de entalpía ≈ ΔP0rel/ρ = ω̄·½W1²
    dh = om_total * 0.5 * W1 ** 2
    detail = dict(profile=om_profile, endwall=om_endwall, shock=om_shock,
                  deq=deq, CL=CL, f_re=f_re, re_c=re_c)
    return om_total, dh, detail


def _solve_axial_mach(mdot: float, A: float, T0: float, P0: float,
                      kb: float, swirl_tan: float = 0.0,
                      gam: float | None = None,
                      cp: float | None = None) -> tuple[float, bool]:
    """Mach axial subsónico desde continuidad ṁ = ρ(M)·A·KB·Cx(M).

    swirl_tan = tanα del remolino en la estación (C = Cx·√(1+tan²α)).
    `gam`/`cp` locales (gas imperfecto, fase 9); None = valores de
    referencia. Devuelve (Mx, choked). Si no hay raíz subsónica devuelve
    el Mx del máximo de gasto (garganta) y choked=True — g queda continuo.
    """
    gam = gamma_air(T0) if gam is None else gam
    cp = cp_air(T0) if cp is None else cp
    a0 = math.sqrt(gam * RGAS * T0)
    rho0 = P0 / (RGAS * T0)
    k = 1.0 + swirl_tan ** 2

    def flow(mx: float) -> float:
        c = mx * a0  # aproximación: Mx referido a a0; corregido abajo
        # iteración corta densidad-estática
        T = T0
        for _ in range(4):
            C2 = k * c * c
            T = max(T0 - C2 / (2 * cp), 0.5 * T0)
            a = math.sqrt(gam * RGAS * T)
            c = mx * a
        rho = rho0 * (T / T0) ** (1.0 / (gam - 1))
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
                      swirl_tan: float = 0.0, gam: float | None = None
                      ) -> tuple[float, float, float, float]:
    """(T, P, rho, Cx) estáticos desde Mach axial + remolino."""
    gam = gamma_air(T0) if gam is None else gam
    k = 1.0 + swirl_tan ** 2
    T = T0 / (1.0 + 0.5 * (gam - 1) * k * mx * mx)
    a = math.sqrt(gam * RGAS * T)
    Cx = mx * a
    P = P0 * (T / T0) ** (gam / (gam - 1))
    rho = P / (RGAS * T)
    return T, P, rho, Cx


# ==========================================================================
# 2a. Capacidad de subida de presión al stall — Koch 1981 (fase 9)
#
# Sustituye a la constante CH_STALL_MAX = 0.48 del proxy anterior. Koch
# (ASME JEP 103, 1981) correlaciona el coeficiente de subida de presión
# estática EFECTIVO al stall de una etapa con la geometría de la cascada
# a través del parámetro de difusor L/g₂ (longitud de la línea media del
# arco circular / espaciado escalonado de salida), con correcciones por
# Reynolds, holgura de punta y espaciado axial, y con un factor de cabeza
# dinámica EFECTIVA 𝔉_ef que captura el efecto (grande) del tipo de
# triángulo de velocidades. La etapa entra en stall cuando
#
#     Ch_ef = Ch/𝔉_ef  ≥  Ch_ef,stall(L/g₂) · f_Re · f_ε · f_Δz
#
# Ajustes de las figuras del paper (documentados uno a uno abajo). El
# rango de validez declarado por Koch: L/g₂ ∈ [0.65, 3.5], radio de cubo
# 0.40-0.92, M de entrada 0.1-1.6, AR 0.35-5.0.
# ==========================================================================
KOCH_CH_A = 0.35          # Ch_ef,stall en L/g₂ = 1 (Fig. 15/18)
KOCH_CH_B = 0.145         # pendiente en ln(L/g₂) (Fig. 15/18)
KOCH_CH_MIN, KOCH_CH_MAX = 0.20, 0.58
KOCH_RE_REF = 1.3e5       # Re de referencia de la correlación (Fig. 4)
KOCH_EPS_G_REF = 0.055    # ε/g de referencia (Fig. 5)
KOCH_DZ_S_REF = 0.38      # Δz/s de referencia (Fig. 6)
# Ajuste de la rama de holgura CERO de Koch & Smith 1976 Fig. 8
# ((2δ*/g)|₀ vs Ch/Ch_max). La figura se lee con incertidumbre de ±0.03
# en ordenada; EW_A/EW_N son los dos parámetros del ajuste y se calibran
# globalmente contra las cuatro máquinas de la campaña (documentado en
# docs/VALIDATION.md §3). EW_X_MAX corta la extrapolación por encima del
# stall, donde la figura ya no tiene datos.
EW_A, EW_N, EW_X_MAX = 0.16, 3.0, 1.05


def koch_L_over_g2(camber_deg: float, sigma: float,
                   beta2_deg: float) -> float:
    """Parámetro de difusión L/g₂ de la cascada (Koch 1981, Fig. 13).

    L = c·(θ_rad/2)/sin(θ_rad/2)   longitud de la línea media del arco
                                   circular de comba θ
    g₂ = s·cos(χ₂) = (c/σ)·cos(β₂) espaciado escalonado en el borde de
                                   salida

    ⇒ L/g₂ = [(θ/2)/sin(θ/2)] · σ / cos(β₂)

    En la línea media se usa el GIRO DEL FLUJO como proxy de la comba (la
    comba real la resuelve la capa 5a con incidencia y desviación; el
    proxy la subestima unos grados, del lado conservador porque baja
    L/g₂ y por tanto la capacidad de subida de presión).
    """
    th = abs(math.radians(camber_deg))
    arc = 1.0 if th < 1e-6 else (0.5 * th) / math.sin(0.5 * th)
    cb2 = max(math.cos(math.radians(beta2_deg)), 0.15)
    return float(np.clip(arc * max(sigma, 0.2) / cb2, 0.30, 5.0))


def koch_f_ef(alpha_deg: float, beta_deg: float, U: float, V: float) -> float:
    """Factor de cabeza dinámica efectiva 𝔉_ef (Koch 1981, Fig. 13).

        𝔉_ef = (V² + 2.5·V_min² + 0.5·U²) / (4·V²)

    con V la velocidad de entrada de la fila en SU marco y V_min la
    mínima velocidad posible de entrada en presencia de estelas y capas
    límite de pared:

        V_min = V·sin(α+β)   si (α+β) ≤ 90°
        V_min = V            si (α+β) > 90°
        V_min = U            si el ángulo del OTRO marco es negativo

    Es el mecanismo que explica por qué las etapas de bajo stagger (y
    (α+β) pequeño) stallean antes: no pueden re-energizar el fluido de
    baja cantidad de movimiento que les llega. Rango medido: 0.78-1.11.
    """
    V = max(abs(V), 1e-3)
    ang = alpha_deg + beta_deg
    if min(alpha_deg, beta_deg) < 0.0:
        vmin2 = (U / V) ** 2
    elif ang > 90.0:
        vmin2 = 1.0
    else:
        vmin2 = math.sin(math.radians(ang)) ** 2
    return float(np.clip((1.0 + 2.5 * vmin2 + 0.5 * (U / V) ** 2) / 4.0,
                         0.55, 1.30))


def koch_ch_stall(L_over_g2: float, re_c: float | None = None,
                  eps_over_g: float = KOCH_EPS_G_REF,
                  dz_over_s: float = KOCH_DZ_S_REF) -> float:
    """Ch_ef,stall corregido (Koch 1981, Figs. 4-7 y 15/18)."""
    base = KOCH_CH_A + KOCH_CH_B * math.log(max(L_over_g2, 0.30))
    base = min(max(base, KOCH_CH_MIN), KOCH_CH_MAX)
    # Reynolds (Fig. 4): 0.80 a Re 1e4, 1.00 a 1.3e5, 1.05 a 1e6 (satura)
    if re_c and re_c > 0.0:
        lr = math.log10(max(re_c, 1e3))
        if lr <= math.log10(KOCH_RE_REF):
            f_re = 0.80 + 0.20 * (lr - 4.0) / (math.log10(KOCH_RE_REF) - 4.0)
        else:
            f_re = 1.0 + 0.05 * (lr - math.log10(KOCH_RE_REF)) / \
                (6.0 - math.log10(KOCH_RE_REF))
        f_re = min(max(f_re, 0.75), 1.05)
    else:
        f_re = 1.0
    # Holgura de punta (Fig. 5): 1.22 en ε/g=0, 1.00 en 0.055, 0.87 en 0.15
    e = min(max(eps_over_g, 0.0), 0.30)
    f_eps = 1.22 - 4.0 * e if e <= KOCH_EPS_G_REF \
        else 1.0 - 1.37 * (e - KOCH_EPS_G_REF)
    f_eps = min(max(f_eps, 0.55), 1.25)
    # Espaciado axial (Fig. 6): 1.08 en Δz/s=0.2, 1.00 en 0.38, 0.99 en 0.8
    dz = min(max(dz_over_s, 0.05), 1.0)
    f_dz = 1.0 + 0.44 * (KOCH_DZ_S_REF - dz) if dz <= KOCH_DZ_S_REF \
        else 1.0 - 0.024 * (dz - KOCH_DZ_S_REF)
    f_dz = min(max(f_dz, 0.95), 1.12)
    return float(base * f_re * f_eps * f_dz)


def endwall_blockage(x_load: float, eps_over_g: float) -> float:
    """Suma de espesores de desplazamiento de pared 2δ*/g (Koch & Smith
    1976, ec. 3 y Fig. 8).

        2δ*/g = (2δ*/g)|_{ε/g=0} + 2·(ε/g)·(Ch/Ch_max)

    con el término de holgura cero ajustado a la Fig. 8 como 0.25·x⁴
    (x = Ch/Ch_stall): plano lejos del stall y con la subida abrupta
    medida al acercarse a él. Ésta es la fuente del BLOQUEO del annulus,
    que así CRECE con la carga y con la holgura en vez de ser la recta
    inventada 0.98 − 0.005·i. Es también el mecanismo que hace que las
    etapas traseras de una máquina multietapa se descuelguen.
    """
    x = min(max(x_load, 0.0), EW_X_MAX)
    e = min(max(eps_over_g, 0.0), 0.30)
    return float(EW_A * x ** EW_N + 2.0 * e * x)


def _incidence_bucket(inc_deg: float, m1_rel: float) -> float:
    """Multiplicador de pérdida de perfil por incidencia fuera de diseño.

    Bucket parabólico con SEMIANCHO dependiente del Mach relativo de
    entrada: W = W_LOW → W_HIGH linealmente entre M 0.2 y 0.8 (Aungier
    2003 — el rango de incidencia útil se estrecha drásticamente con M;
    a bajo M reproduce el bucket clásico de ±10°). Rama de incidencia
    NEGATIVA 1.5× más tolerante (práctica estándar: el lado de choke del
    bucket es más ancho que el de stall). Continuo en inc y M; tope ×4
    (bucket saturado — el punto ya está profundamente separado).
    """
    m = min(max(m1_rel, 0.0), 1.2)
    f = min(max((m - 0.2) / 0.6, 0.0), 1.0)
    w = W_INC_LOW_DEG - (W_INC_LOW_DEG - W_INC_HIGH_DEG) * f
    if inc_deg < 0.0:
        w *= 1.5
    return min(1.0 + (inc_deg / w) ** 2, 4.0)


def _meanline(theta: np.ndarray, rpm: float | None = None,
              mdot: float | None = None,
              frozen: dict | None = None,
              vsv_deg: float = 0.0,
              per_stage: dict | None = None,
              bleed: dict | None = None) -> dict:
    """Stage-stacking 1-D en la línea media, termodinámicamente consistente.

    Punto de diseño: `_meanline(theta)` — los ángulos de flujo salen de
    (φᵢ, ψᵢ, Rxᵢ) y el álabe se talla a incidencia cero (convención del
    generador de geometría). Desde la fase 8 las distribuciones por etapa
    son pendientes lineales (θ[13]=phi_slope, θ[14]=Rx_slope, mismo
    esquema ξᵢ que psi_slope); con pendientes 0 colapsan bit-exacto a los
    escalares. `per_stage` (modo experto, nivel 8b): dict opcional
    {"phi": [...], "psi": [...], "Rx": [...]} de longitud n_stages que
    SOBREESCRIBE las distribuciones — no entra al espacio de búsqueda
    (se usa desde --eval-theta --per-stage y la validación E³).

    Fuera de diseño: pasar `rpm`/`mdot` y `frozen` (dict devuelto en
    record["frozen"] del diseño: ángulos de flujo y áreas congelados); la
    incidencia paga el bucket W(M) de `_incidence_bucket`, y el sub-giro
    por desviación crece con la incidencia positiva (K_DEV_OFFDESIGN).
    `vsv_deg` > 0 cierra los estátores variables frontales (más pre-swirl
    → menos incidencia del rotor a baja velocidad): ángulo pleno en la
    etapa 0 decayendo linealmente a 0 hacia la etapa n/2. En diseño
    (i = 0, vsv = 0) nada de esto altera el resultado.

    Convención de triángulos (Dixon, ángulos desde el eje, Cx constante
    en el diseño): Cu1/U = 1−Rx−ψ/2, Cu2/U = 1−Rx+ψ/2,
    tanβ = (U−Cu)/Cx, tanα = Cu/Cx.
    """
    theta = pad_theta(theta)
    (n_st_f, RPM_d, HTR, phi1, psi_mid, psi_slope, Rx,
     sigma_r, sigma_s, AR, T0_in, P0_in, mdot_d,
     phi_slope, Rx_slope) = [float(x) for x in theta]
    # clamp de seguridad holgado (los bounds del optimizador son 1-8; la
    # validación evalúa máquinas de hasta 10+ etapas fuera de bounds)
    n_st = int(round(min(max(n_st_f, 1.0), 16.0)))
    RPM = float(rpm) if rpm is not None else RPM_d
    mdot = float(mdot) if mdot is not None else mdot_d
    off_design = frozen is not None
    clipped: list[str] = []

    # --- Distribuciones por etapa (fase 8) --------------------------------
    # ξ_i = 2i/(N−1)−1; con pendiente 0, φ_i = φ1·1.0 (bit-exacto con el
    # espacio 13-D legacy). Clamps suaves para θ degenerado (g continuo).
    def _xi(i: int) -> float:
        return (2.0 * i / (n_st - 1) - 1.0) if n_st > 1 else 0.0

    phi_st = [max(phi1 * (1.0 + phi_slope * _xi(i)), 0.05)
              for i in range(n_st)]
    Rx_st = [Rx * (1.0 + Rx_slope * _xi(i)) for i in range(n_st)]
    psi_st = [max(psi_mid * (1.0 + psi_slope * _xi(i)), 0.05)
              for i in range(n_st)]
    if per_stage:
        for key, dst in (("phi", phi_st), ("psi", psi_st), ("Rx", Rx_st)):
            if key in per_stage:
                vals = list(per_stage[key])
                if len(vals) != n_st:
                    raise ValueError(
                        f"per_stage['{key}'] debe tener n_stages={n_st} "
                        f"valores — llegaron {len(vals)}")
                lo = 0.05 if key in ("phi", "psi") else -2.0
                dst[:] = [max(float(v), lo) for v in vals]

    # --- Sangrado por etapa (fase 9) --------------------------------------
    # Fracción del gasto de ENTRADA extraída AGUAS ABAJO de cada etapa.
    # Hasta la fase 9 el puerto de sangrado existía solo en la carcasa
    # impresa: la física pasaba el mismo ṁ de la entrada a la salida. El
    # sangrado intermedio es, junto con los VSV, el mecanismo por el que
    # una máquina de PR alto arranca y opera a velocidad parcial
    # (Dixon & Hall §5.9) — sin él el mapa de baja velocidad no existe.
    bleed_st = [0.0] * n_st
    if bleed:
        for k, v in bleed.items():
            i = int(k)
            if 0 <= i < n_st:
                bleed_st[i] = float(np.clip(float(v), 0.0, 0.25))

    omega = RPM * 2 * math.pi / 60.0
    cp_in = cp_air(T0_in)
    gam_in = gamma_air(T0_in)
    a0_in = math.sqrt(gam_in * RGAS * T0_in)
    rho0_in = P0_in / (RGAS * T0_in)

    # --- Annulus de entrada: r_tip DERIVADO de la continuidad --------------
    # Dados (φ_0, HTR, RPM, ṁ): Cx = φ_0·ω·r_mean y A = π·r_tip²(1−HTR²)
    # deben transportar ṁ — punto fijo en r_tip (único y contractivo). La
    # etapa 0 manda: su φ define el annulus de entrada (con pendientes, la
    # etapa frontal corre a φ1·(1−s_φ)).
    tan_a1_first = (1.0 - Rx_st[0] - psi_st[0] / 2.0) / phi_st[0]
    r_tip = 0.20
    for _ in range(80):
        r_mean_i = 0.5 * (1.0 + HTR) * r_tip
        Cx_i = phi_st[0] * omega * r_mean_i
        C1_i = Cx_i * math.sqrt(1.0 + tan_a1_first ** 2)
        if C1_i > 0.95 * a0_in:
            C1_i = 0.95 * a0_in
        T1_i = max(T0_in - C1_i ** 2 / (2 * cp_in), 0.5 * T0_in)
        rho1_i = rho0_in * (T1_i / T0_in) ** (1.0 / (gam_in - 1))
        A_req = mdot / max(rho1_i * Cx_i * BLOCKAGE_INIT, 1e-6)
        r_tip_new = math.sqrt(A_req / (math.pi * max(1.0 - HTR ** 2, 1e-3)))
        if abs(r_tip_new - r_tip) < 1e-9:
            r_tip = r_tip_new
            break
        r_tip = 0.5 * (r_tip + r_tip_new)
    if phi_st[0] * omega * 0.5 * (1.0 + HTR) * r_tip * \
            math.sqrt(1.0 + tan_a1_first ** 2) > 0.95 * a0_in:
        clipped.append("inlet_C1_sonic")
    r_hub = HTR * r_tip
    r_mean = 0.5 * (r_hub + r_tip)          # línea media (const en const_mean)
    A1 = math.pi * (r_tip ** 2 - r_hub ** 2)
    if off_design:
        A1 = frozen["areas"][0]

    T0, P0 = T0_in, P0_in
    kb = BLOCKAGE_INIT
    mdot_i = mdot                           # gasto que ve la etapa i
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
    ds_total = 0.0

    # Ángulo de entrada de la primera etapa (IGV implícito = pre-swirl del
    # triángulo repetitivo de la etapa 1; su pérdida se carga a la etapa 1).
    for i in range(n_st):
        U = omega * (frozen["r_m"][i] if off_design and i < len(frozen["r_m"])
                     else r_mean)
        wdf = _HOWELL_WDF[min(i, len(_HOWELL_WDF) - 1)]
        psi_i = psi_st[i]
        Rx_i = Rx_st[i]

        if off_design:
            tan_a1 = frozen["alpha1"][i]
            # VSV: cerrar los estátores frontales añade pre-swirl al rotor
            # siguiente (menos incidencia a baja N). Ángulo pleno en la
            # etapa 0, decayendo linealmente a 0 hacia la etapa n/2.
            if vsv_deg:
                n_half = max(n_st // 2, 1)
                fade = max(1.0 - i / n_half, 0.0)
                if fade > 0.0:
                    a1_new = math.atan(tan_a1) \
                        + math.radians(vsv_deg * fade)
                    tan_a1 = math.tan(min(a1_new, 1.45))
            tan_b2_fr = frozen["beta2"][i]
            A_in = frozen["areas"][i]
        else:
            tan_a1 = (1.0 - Rx_i - psi_i / 2.0) / phi_st[i]
            tan_b2_fr = None
            A_in = frozen_out["areas"][i]

        # Propiedades locales del gas (fase 9): el aire de la etapa 10 de
        # un HPC está a ~700 K, donde cp es 7% mayor y γ 3% menor.
        cp_i = cp_air(T0)
        gam_i = gamma_air(T0)

        # Continuidad en la entrada de la etapa (con remolino α1)
        mx, chk = _solve_axial_mach(mdot_i, A_in, T0, P0, kb,
                                    swirl_tan=tan_a1, gam=gam_i, cp=cp_i)
        if chk or mx > MX_CHOKE:
            choke_any = True
            g_soft_pen += max(mx / MX_CHOKE - 1.0, 0.0) + (0.5 if chk else 0.0)
            clipped.append(f"choke_st{i}")
        T1, P1, rho1, Cx = _static_from_mach(mx, T0, P0, swirl_tan=tan_a1,
                                             gam=gam_i)
        phi = Cx / max(U, 1e-3)

        inc_deg = 0.0
        if off_design:
            # geometría congelada: α1 y β2 metálicos fijos → ψ emergente.
            # Desviación fuera de diseño (Creveling 1968 / SP-36): con
            # incidencia POSITIVA el flujo sub-gira Δδ = k·i⁺ y la ψ
            # lograda cae — la pendiente de la speedline cerca de surge
            # se aplana de forma realista. En diseño i=0 ⇒ Δδ=0 (punto
            # de diseño invariante). Clamps para θ degenerados: Δδ ≤ 10°,
            # β2 ≤ 83° — g queda continuo.
            tan_b1_fr = (U - Cx * tan_a1) / max(Cx, 1e-3)
            inc_deg = math.degrees(math.atan(tan_b1_fr)) \
                - frozen["beta1_deg"][i]
            ddev = min(K_DEV_OFFDESIGN * max(inc_deg, 0.0), 10.0)
            b2_eff = min(math.atan(tan_b2_fr) + math.radians(ddev), 1.45)
            tan_b2 = math.tan(b2_eff)
            psi_eff = 1.0 - phi * (tan_a1 + tan_b2)
            psi_eff = float(np.clip(psi_eff, -0.2, 0.9))
            Cu1 = Cx * tan_a1
            Cu2 = U - Cx * tan_b2
        else:
            Cu1 = (1.0 - Rx_i - psi_i / 2.0) * U
            Cu2 = (1.0 - Rx_i + psi_i / 2.0) * U
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

        a1_sound = math.sqrt(gam_i * RGAS * T1)
        # Punta local de esta etapa; Mach relativo con free-vortex Cu·r=const
        h_blade = A_in / (2 * math.pi * r_mean)
        r_tip_i = r_mean + 0.5 * h_blade
        r_hub_i = r_mean - 0.5 * h_blade
        # Punta y raíz con el torbellino de exponente VORTEX_N: la
        # velocidad axial deja de ser radialmente constante salvo en el
        # vórtice libre (equilibrio radial simple, forma cerrada).
        U_t = omega * r_tip_i
        Cu1_t = vortex_cu(Cu1, r_tip_i, r_mean)
        Cx_t = vortex_cx(Cx, Cu1, r_tip_i, r_mean)
        W1_tip = math.sqrt(Cx_t ** 2 + (U_t - Cu1_t) ** 2)
        M_rel_tip = W1_tip / a1_sound
        if i == 0:
            M_rel_tip1 = M_rel_tip
            U_tip1 = U_t
        M1_rel_mean = W1 / a1_sound

        # Reacción de raíz — la restricción dura g[8]
        Rx_hub = vortex_reaction(1.0 - (Cu1 + Cu2) / max(2.0 * U, 1e-6),
                                 r_hub_i, r_mean)

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
        # Solidez REAL de la pieza: el nº de álabes es entero, así que la σ
        # que se fabrica no es la σ continua que optimiza el buscador. Se
        # reporta el delta (fase 9, H9) en vez de dejarlo silencioso.
        sig_r_act = n_bl_r * c_rotor / max(2 * math.pi * r_mean, 1e-9)
        sig_s_act = n_bl_s * c_stator / max(2 * math.pi * r_mean, 1e-9)

        h_c_r = h_blade / max(c_rotor, 1e-4)
        h_c_s = h_blade / max(c_stator, 1e-4)

        # Pérdidas por fila. Re de cuerda con el estado de ENTRADA de la
        # etapa (ρ1 también para el estátor: ρ sube fila a fila, así que
        # subestima ligeramente su Re — lado conservador).
        Re_r = rho1 * W1 * c_rotor / MU_AIR
        Re_s = rho1 * C2 * c_stator / MU_AIR
        om_r, dh_r, det_r = _row_losses(b1, b2r, W1, W2, sigma_r, h_c_r,
                                        M_rel_tip, M1_rel_mean, re_c=Re_r,
                                        gam=gam_i)
        M2_abs = C2 / a1_sound
        om_s, dh_s, det_s = _row_losses(a2r, a1r, C2, C1, sigma_s, h_c_s,
                                        M2_abs, M2_abs, re_c=Re_s, gam=gam_i)
        # Fuera de diseño: bucket de incidencia sobre el rotor con
        # semiancho dependiente del Mach (Aungier 2003) — la incidencia
        # ya se calculó arriba (misma β1: tanβ1 = (U−Cu1)/Cx).
        f_inc = _incidence_bucket(inc_deg, M1_rel_mean) if off_design else 1.0
        om_r *= f_inc
        # Holgura de punta (solo rotor): ε ABSOLUTA en mm ⇒ ε/h crece
        # hacia las etapas traseras (h cae) — pérdida por fila, no débito
        # uniforme. Modelo por tramos de _clearance_loss_frac.
        eps_h = min(TIP_CLEARANCE_MM / max(h_blade * 1000.0, 1e-3),
                    EPS_H_MAX)
        dh_cl = _clearance_loss_frac(eps_h) * dh0 if dh0 > 0 else 0.0

        # ================= MARCHA TERMODINÁMICA EXACTA DE LA ETAPA ========
        # (fase 9) Antes: Δh_pérdida = ω̄·½W₁² = ΔP₀/ρ a la ENTRADA de la
        # fila, y PR_i = (1+η·Δh₀/(cp·T₀))^(γ/(γ−1)). Eso subestima la
        # pérdida de las etapas calientes en el factor T₀₃/T₀_fila (15-25%
        # por etapa, acumulativo). Dixon & Hall §5.5 ec. (5.4)-(5.9): el
        # equivalente en TRABAJO de una pérdida es T₀₃·Δs, con Δs de los
        # coeficientes Yp de cada fila. Aquí se marcha fila a fila con las
        # presiones totales REALES (relativa en el rotor, absoluta en el
        # estátor), así que PR sale de la marcha en vez de una fórmula.
        kx = gam_i / (gam_i - 1.0)
        # -- Rotor: total relativo a la entrada. ω̄ está referida a la
        #    cabeza dinámica COMPRESIBLE (P₀ − p) — la definición de Koch
        #    & Smith 1976 (nomenclatura: ω̄ = ΔP/(P−p)₁). Usar ½ρW² en su
        #    lugar subestima la pérdida un 60% a M_rel ≈ 1.4, justo donde
        #    viven los rotores transónicos de la campaña de validación.
        T01rel = T1 + W1 ** 2 / (2.0 * cp_i)
        P01rel = P1 * (T01rel / max(T1, 1.0)) ** kx
        q1rel = max(P01rel - P1, 1e-3)
        dP0rel = om_r * q1rel
        P02rel = max(P01rel - dP0rel, 0.05 * P01rel)
        ds_r = RGAS * math.log(P01rel / P02rel)
        # -- Salida del rotor (estación 2): T₀₂ del trabajo de Euler con
        #    cp variable; el total RELATIVO se conserva a radio constante
        T02_i = T_from_h(h_air(T0) + dh0, T0 + dh0 / max(cp_i, 1.0))
        T2_st = max(T02_i - C2 ** 2 / (2.0 * cp_i), 0.5 * T0)
        P2_st = P02rel * (T2_st / max(T01rel, 1.0)) ** kx
        P02_i = P2_st * (T02_i / max(T2_st, 1.0)) ** kx
        rho2_st = max(P2_st / (RGAS * T2_st), 1e-3)
        # -- Estátor: sin trabajo, T₀₃ = T₀₂ (misma referencia compresible)
        dP0_s = om_s * max(P02_i - P2_st, 1e-3)
        P03_i = max(P02_i - dP0_s, 0.05 * P02_i)
        ds_s = RGAS * math.log(P02_i / P03_i)
        T0_out_i = T02_i
        # -- Fuga de punta: se contabiliza como generación de entropía a
        #    la temperatura de salida de la etapa (mismo criterio)
        ds_cl = dh_cl / max(T0_out_i, 1.0)
        P0_fs = P03_i * math.exp(-ds_cl / RGAS)     # aún sin endwall
        ds_fs = ds_r + ds_s + ds_cl
        eta_fs = float(np.clip((dh0 - T0_out_i * ds_fs) / max(dh0, 1e-3),
                               0.05, 1.0)) if dh0 > 1.0 else 0.05
        if dh0 <= 1.0:
            clipped.append(f"psi_collapse_st{i}")

        # ================= CAPACIDAD DE SUBIDA DE PRESIÓN (Koch 1981) =====
        # Ch = subida de presión estática ENTALPÍA-EQUIVALENTE (isentrópica
        # de p₁ a p₃) sobre la suma de cabezas dinámicas de entrada de las
        # dos filas. Antes se usaba Δh₀ (el trabajo REAL) en el numerador,
        # que sobreestima Ch en 1/η ≈ 1.12 y por tanto inventa stall.
        T3_st = max(T0_out_i - C1 ** 2 / (2.0 * cp_i), 0.5 * T0)
        p3_st = P0_fs * (T3_st / max(T0_out_i, 1.0)) ** kx
        Ts_id = T_from_phi(phi_air(T1) + RGAS * math.log(
            max(p3_st, 1.0) / max(P1, 1.0)), T1 + 20.0)
        dh_static_id = max(h_air(Ts_id) - h_air(T1), 0.0)
        q_r = 0.5 * rho1 * W1 ** 2
        q_s = 0.5 * rho2_st * C2 ** 2
        Ch = dh_static_id / max(0.5 * (W1 ** 2 + C2 ** 2), 1e-3)
        # Cabeza dinámica EFECTIVA: media de las dos filas pesada por su
        # cabeza dinámica de entrada (criterio de promediado de Koch)
        fef_r = koch_f_ef(math.degrees(a1r), math.degrees(b1), U, W1)
        fef_s = koch_f_ef(math.degrees(a2r), math.degrees(b2r), U, C2)
        w_sum = max(q_r + q_s, 1e-9)
        fef = (fef_r * q_r + fef_s * q_s) / w_sum
        # Parámetro de difusión de la cascada, mismo promediado
        lg_r = koch_L_over_g2(math.degrees(b1 - b2r), sig_r_act,
                              math.degrees(b2r))
        lg_s = koch_L_over_g2(math.degrees(a2r - a1r), sig_s_act,
                              math.degrees(a1r))
        L_g2 = (lg_r * q_r + lg_s * q_s) / w_sum
        # ε/g y Δz/s con la geometría real de la fila de rotor
        s_pitch = 2.0 * math.pi * r_mean / max(n_bl_r, 1)
        beta_m_r = 0.5 * (b1 + b2r)
        g_stag = max(s_pitch * math.cos(beta_m_r), 1e-4)
        eps_g = (TIP_CLEARANCE_MM / 1000.0) / g_stag
        dz_s = ROW_GAP_FRACTION * c_rotor / max(s_pitch, 1e-6)
        Ch_stall = koch_ch_stall(L_g2, re_c=Re_r, eps_over_g=eps_g,
                                 dz_over_s=dz_s)
        Ch_ef = Ch / max(fef, 0.3)
        x_load = Ch_ef / max(Ch_stall, 1e-3)
        SM_i = 1.0 - x_load

        # ================= BLOQUEO DE PARED (Koch & Smith 1976) ===========
        # El bloqueo del annulus deja de ser la recta 0.98 − 0.005·i y pasa
        # a EMERGER de la carga: 2δ*/g crece como x⁴ al acercarse al stall
        # y linealmente con la holgura (ec. 3 y Fig. 8 del paper). Es el
        # mecanismo que descuelga a las etapas traseras de un multietapa.
        d_star_2 = endwall_blockage(x_load, eps_g) * g_stag
        a_ew = float(np.clip(d_star_2 / max(h_blade, 1e-4), 0.0, 0.35))
        kb_next = float(np.clip(1.0 - a_ew, 0.80, 0.995))
        # ...y su DÉBITO DE RENDIMIENTO, ec. (2) del mismo paper:
        #     η = η̃·(1 − Σδ*/h)/(1 − Σν/h),   Σν ≈ 0.48·Σδ*  (Fig. 10)
        # η̃ es el rendimiento "de flujo libre" (perfil + secundaria +
        # choque + fuga) que acaba de salir de la marcha entrópica. Esto
        # sustituye al término de annulus de Howell CDa = 0.020·s/h, que
        # era un arrastre equivalente calibrado con K_ENDWALL = 1.4.
        f_ew = (1.0 - a_ew) / max(1.0 - NU_OVER_DELTA * a_ew, 1e-3)
        eta_tt = float(np.clip(eta_fs * f_ew, 0.05, 1.0))
        dh_ew = max(dh0 * (eta_fs - eta_tt), 0.0)
        ds_ew = dh_ew / max(T0_out_i, 1.0)
        P0_out_i = P0_fs * math.exp(-ds_ew / RGAS)
        PR_i = P0_out_i / max(P0, 1.0)
        ds_stage = ds_fs + ds_ew
        dh_loss_i = T0_out_i * ds_stage        # equivalente en trabajo
        dh_r = T0_out_i * ds_r                 # desglose coherente
        dh_s = T0_out_i * ds_s

        # Estado estático y área en la SALIDA DEL ROTOR (estación 2) — lo
        # necesita el pasaje L1 para contraer el annulus fila a fila (sin
        # esto TD3 ve el área de entrada en el rotor y sobregira el flujo).
        A_rotor_exit = mdot_i / max(rho2_st * Cx * kb, 1e-6)

        stage_table.append(dict(
            stage=i, phi=phi, psi=psi_eff, Rx=Rx_i, Rx_hub=Rx_hub, U_m=U,
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
            # --- fase 9: Koch / Koch & Smith / gas variable / sangrado ---
            Ch_ef=Ch_ef, Ch_stall=Ch_stall, load_ratio=x_load,
            L_over_g2=L_g2, f_ef=fef, eps_over_g=eps_g,
            blockage_in=kb, blockage_out=kb_next,
            delta_star_sum_mm=d_star_2 * 1000.0,
            cp_J_kgK=cp_i, gamma=gam_i,
            mdot_in_kgs=mdot_i, bleed_frac=bleed_st[i],
            sigma_rotor_actual=sig_r_act, sigma_stator_actual=sig_s_act,
            ds_stage_J_kgK=ds_stage,
            losses=dict(rotor=dh_r, stator=dh_s, clearance=dh_cl,
                        eps_over_h=eps_h,
                        omega_rotor=om_r, omega_stator=om_s,
                        ds_rotor=ds_r, ds_stator=ds_s, ds_clearance=ds_cl,
                        deq_rotor=det_r["deq"], deq_stator=det_s["deq"],
                        Re_rotor=Re_r, Re_stator=Re_s,
                        f_re_rotor=det_r["f_re"], f_re_stator=det_s["f_re"]),
        ))

        dh_shaft_total += dh0 * mdot_i
        dh_loss_total += dh_loss_i * mdot_i
        ds_total += ds_stage
        T0, P0 = T0_out_i, P0_out_i
        kb = kb_next
        mdot_i *= (1.0 - bleed_st[i])

        # Annulus de la siguiente estación por continuidad (Cx de diseño
        # constante; en const_mean, r_mean fijo y h se ajusta)
        if not off_design:
            # Área de la siguiente estación que MANTIENE el Cx de diseño
            # de la etapa que la recibe (φ_{i+1}·Um; con pendientes el Cx
            # de diseño varía etapa a etapa), evaluada con el estado
            # (T0,P0) tras la etapa. Para la última estación (salida tras
            # el OGV, remolino ya quitado) manda el φ de la última etapa.
            tan_a_next = tan_a1 if i < n_st - 1 else 0.0
            Cx_dsg = phi_st[min(i + 1, n_st - 1)] * omega * r_mean
            C_n = Cx_dsg * math.sqrt(1.0 + tan_a_next ** 2)
            cp_n, gam_n = cp_air(T0), gamma_air(T0)
            T_n = max(T0 - C_n ** 2 / (2 * cp_n), 0.5 * T0)
            rho_n = (P0 / (RGAS * T0)) * (T_n / T0) ** (1.0 / (gam_n - 1))
            A_next = mdot_i / max(rho_n * Cx_dsg * kb, 1e-6)
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

    # --- IGV: fila que da el pre-swirl α₁ de la etapa 1 (fase 9) ----------
    # La capa 5a la construye desde la fase 3 y la máquina impresa la
    # lleva, pero hasta ahora NO pagaba pérdida ni longitud: la física
    # contaba n rotores + n estátores + OGV. Es una fila aceleradora
    # (axial → α₁), barata pero no gratis; sin ella el η y la longitud
    # describían una máquina distinta de la que sale por la impresora.
    s_first = stage_table[0]
    a_igv = math.radians(s_first["alpha1_deg"])
    igv = dict(present=False, omega_bar=0.0, dP0_Pa=0.0, length_mm=0.0)
    if abs(a_igv) > 1e-4:
        Cx_igv = s_first["Cx"]
        C_igv_out = s_first["C1"]
        c_igv_m = s_first["chord_stator_mm"] / 1000.0
        T_igv, P_igv, rho_igv, _ = _static_from_mach(
            s_first["Mx"], T0_in, P0_in, swirl_tan=0.0, gam=gam_in)
        Re_igv = rho_igv * C_igv_out * c_igv_m / MU_AIR
        M_igv = C_igv_out / math.sqrt(gam_in * RGAS * max(T_igv, 1.0))
        om_igv, _, _ = _row_losses(
            0.0, a_igv, Cx_igv, C_igv_out, sigma_s,
            s_first["h_blade_mm"] / max(s_first["chord_stator_mm"], 1e-3),
            M_igv, M_igv, re_c=Re_igv, gam=gam_in)
        dP0_igv = om_igv * 0.5 * rho_igv * Cx_igv ** 2
        P0_igv_out = max(P0_in - dP0_igv, 0.5 * P0_in)
        ds_igv = RGAS * math.log(P0_in / P0_igv_out)
        # La caída se aplica al nivel de presión de TODA la máquina (la
        # marcha de etapas partió de P0_in) y su entropía entra al total.
        P0 = P0 * (P0_igv_out / P0_in)
        ds_total += ds_igv
        length_m += (1.0 + ROW_GAP_FRACTION) * c_igv_m
        igv = dict(present=True, omega_bar=float(om_igv),
                   dP0_Pa=float(dP0_igv), length_mm=c_igv_m * 1000.0,
                   ds_J_kgK=float(ds_igv), alpha_out_deg=math.degrees(a_igv))

    # --- OGV: quita el remolino residual α1 → axial (fila de estátor) -----
    s_last = stage_table[-1]
    a_ogv = math.radians(s_last["alpha1_deg"])
    C_ogv_in = s_last["C1"]
    Cx_last = s_last["Cx"]
    cp_e, gam_e = cp_air(T0), gamma_air(T0)
    T_st, _, rho_st, _ = _static_from_mach(s_last["Mx"], T0, P0,
                                           swirl_tan=math.tan(a_ogv),
                                           gam=gam_e)
    Re_ogv = rho_st * C_ogv_in * s_last["chord_stator_mm"] / 1000.0 / MU_AIR
    om_ogv, _, _ = _row_losses(a_ogv, 0.0, C_ogv_in, Cx_last, sigma_s,
                               s_last["h_blade_mm"] /
                               max(s_last["chord_stator_mm"], 1e-3),
                               0.5, 0.5, re_c=Re_ogv, gam=gam_e)
    P0_pre_ogv = P0
    P0 = max(P0 - om_ogv * 0.5 * rho_st * C_ogv_in ** 2, 0.5 * P0)
    ds_total += RGAS * math.log(P0_pre_ogv / P0)
    length_m += (1.0 + ROW_GAP_FRACTION) * s_last["chord_stator_mm"] / 1000.0
    # estado de salida tras OGV (flujo axial)
    mx_e, chk_e = _solve_axial_mach(mdot_i,
                                    frozen_out["areas"][-1] if not off_design
                                    else frozen["areas"][-1], T0, P0, kb, 0.0,
                                    gam=gam_e, cp=cp_e)
    if chk_e:
        choke_any = True
        clipped.append("choke_exit")
    T_e, P_e, rho_e, Cx_e = _static_from_mach(mx_e, T0, P0, 0.0, gam=gam_e)
    M_exit = Cx_e / math.sqrt(gam_e * RGAS * T_e)

    # --- Totales de máquina ------------------------------------------------
    # Rendimientos EXACTOS para gas imperfecto vía la función de entropía
    # phi(T) = ∫cp/T dT (fase 9). En el límite cp = cte estas expresiones
    # colapsan a ((γ−1)/γ)·lnPR/lnτ y (PR^((γ−1)/γ)−1)/(τ−1).
    PR = P0 / P0_in
    tau = T0 / T0_in
    if PR > 1.0 + 1e-9 and tau > 1.0 + 1e-9:
        dphi = phi_air(T0) - phi_air(T0_in)
        eta_poly = RGAS * math.log(PR) / max(dphi, 1e-9)
        T_isen = T_from_phi(phi_air(T0_in) + RGAS * math.log(PR), T0)
        eta_isen = (h_air(T_isen) - h_air(T0_in)) / \
            max(h_air(T0) - h_air(T0_in), 1e-6)
    else:
        eta_poly, eta_isen = 0.05, 0.05
    eta_poly = float(np.clip(eta_poly, 0.05, 1.0))
    eta_isen = float(np.clip(eta_isen, 0.05, 1.0))
    # Potencia con el gasto REAL de cada etapa (el sangrado se lleva su
    # parte del trabajo ya hecho): dh_shaft_total viene acumulado en W/(kg/s)
    power = dh_shaft_total
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
    Rx_hub_min = min(s["Rx_hub"] for s in st)
    stage_stall_first = int(min(st, key=lambda s: s["SM"])["stage"])

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
        # fase 9: reacción de raíz. Con vórtice libre y HTR bajo, Cu·r=cte
        # dispara Cu en el cubo y la reacción se hunde (flujo que ACELERA
        # en el rotor del hub → separación garantizada). Se calculaba, se
        # reportaba y se le pasaba al surrogate desde v0.1 — pero no
        # restringía nada.
        (RX_HUB_MIN - Rx_hub_min) / max(RX_HUB_MIN, 1e-3),
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
        eps_tip_mm=TIP_CLEARANCE_MM,
        M_rel_tip1=_safe(M_rel_tip1, 9.9), M_exit=_safe(M_exit, 9.9),
        max_DF_rotor=_safe(max_DF_r, 9.9), max_DF_stator=_safe(max_DF_s, 9.9),
        min_dehaller_rotor=_safe(min_dh_r, 0.0),
        min_dehaller_stator=_safe(min_dh_s, 0.0),
        min_SM=_safe(min_SM, -9.9), h_blade_last_mm=_safe(h_last, 0.0),
        Rx_hub_min=_safe(Rx_hub_min, -9.9),
        # --- fase 9: trazabilidad del nuevo modelo ------------------------
        stage_stall_first=stage_stall_first,
        ds_machine_J_kgK=_safe(ds_total, 0.0),
        blockage_exit=_safe(kb, 0.9),
        igv=igv,
        bleed=dict(per_stage=list(bleed_st),
                   total_frac=float(1.0 - mdot_i / max(mdot, 1e-9)),
                   mdot_exit_kgs=float(mdot_i)),
        gas_model="variable_cp" if GAS_VARIABLE else "perfect",
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
              frozen: dict | None = None, vsv_deg: float = 0.0) -> dict:
    """Evalúa el diseño theta en otro punto (RPM, ṁ) con la geometría
    congelada del diseño (ángulos de flujo y áreas). `vsv_deg` cierra los
    estátores variables frontales (ver _meanline)."""
    theta = pad_theta(theta)
    if frozen is None:
        frozen = _meanline(theta)["frozen"]
    frozen = dict(frozen)
    frozen.setdefault("phi_d", float(theta[3]))
    return _meanline(theta, rpm=rpm, mdot=mdot, frozen=frozen,
                     vsv_deg=vsv_deg)


def working_line_mdot(PR: float, PR_d: float, mdot_d: float,
                      eta_poly: float, gam: float = GAMMA) -> float:
    """Gasto de la LÍNEA DE TRABAJO para una relación de presión dada.

    Dixon & Hall §5.9, ec. (5.26b): con una tobera/estrangulador de salida
    de área fija y AHOGADA aguas abajo (el caso de un turborreactor y el
    de casi todo banco de ensayo), el gasto adimensional de entrada queda
    ligado al PR:

        ṁ√(c_p T₀₁)/(D²p₀₁) = C · (p₀ₑ/p₀₁)^(1 − (γ−1)/(2γη_p))

    A condiciones de entrada y diámetro fijos eso es ṁ ∝ PR^k. La
    constante C se fija exigiendo que la línea pase por el PUNTO DE
    DISEÑO. Sin esta línea la frase "margen de bombeo" no tiene
    denominador: el margen se mide SIEMPRE respecto al punto de trabajo,
    no respecto a un punto cualquiera de la speedline.
    """
    k = 1.0 - (gam - 1.0) / (2.0 * gam * max(eta_poly, 0.3))
    return float(mdot_d * (max(PR, 1.0) / max(PR_d, 1.0)) ** k)


def _choke_mdot(theta: np.ndarray, rpm: float, frozen: dict,
                mdot_d: float, vsv_deg: float = 0.0,
                hi_frac: float = 1.6) -> float:
    """Gasto máximo que la máquina PUEDE pasar a esa velocidad.

    Antes el meanline seguía calculando con más gasto del que el annulus
    admite y devolvía relaciones de presión que caían por debajo de 1
    (una máquina que "expande"). Físicamente la speedline se hace
    VERTICAL en el choke: el gasto no puede crecer más. Se localiza por
    bisección sobre la bandera de choke.
    """
    def chokes(m: float) -> bool:
        r = _meanline(theta, rpm=rpm, mdot=m, frozen=frozen, vsv_deg=vsv_deg)
        return bool(r["choke_flag"])

    lo = 0.25 * mdot_d
    hi = hi_frac * mdot_d
    if not chokes(hi):
        return hi
    if chokes(lo):
        return lo
    for _ in range(28):
        mid = 0.5 * (lo + hi)
        if chokes(mid):
            hi = mid
        else:
            lo = mid
    return lo


SURGE_SCAN_N = 36


def _surge_mdot(theta: np.ndarray, rpm: float, frozen: dict,
                mdot_d: float, mdot_choke: float,
                vsv_deg: float = 0.0) -> tuple[float, dict | None]:
    """Gasto de bombeo de una speedline. DOS criterios, gana el que se
    encuentra primero viniendo del choke hacia abajo:

      (a) Capacidad de subida de presión agotada — alguna etapa alcanza su
          Ch de stall de Koch 1981 (min_SM ≤ 0).
      (b) Pendiente de la speedline nula o positiva — a la izquierda del
          pico de PR la operación es inestable frente a la resistencia del
          estrangulador (Dixon & Hall §5.11, Fig. 5.14: "surge appears to
          occur when this slope is zero or even a little negative").

    El criterio (b) es imprescindible: el crecimiento de pérdida por
    incidencia puede aplanar la línea ANTES de que Ch llegue a su límite,
    y entonces (a) solo nunca dispara.
    """
    lo = 0.18 * mdot_d
    hi = max(mdot_choke, lo * 1.05)
    ms = np.linspace(hi, lo, SURGE_SCAN_N)      # de choke hacia abajo
    recs = [_meanline(theta, rpm=rpm, mdot=float(m), frozen=frozen,
                      vsv_deg=vsv_deg) for m in ms]
    prs = [r["PR"] for r in recs]
    i_koch = next((j for j, r in enumerate(recs) if r["min_SM"] <= 0.0), None)
    # pico de PR recorriendo hacia gastos decrecientes: el primer punto a
    # partir del cual PR deja de crecer
    i_peak = None
    for j in range(1, len(prs)):
        if prs[j] <= prs[j - 1]:
            i_peak = j - 1
            break
    cands = [j for j in (i_koch, i_peak) if j is not None]
    if not cands:
        return float(ms[-1]), recs[-1]
    j = min(cands)                              # el de MAYOR gasto
    return float(ms[j]), recs[j]


def surge_margin(theta: np.ndarray, record: dict | None = None,
                 vsv_deg: float = 0.0) -> dict:
    """Margen de bombeo del diseño EN LA LÍNEA DE TRABAJO.

    Definición industrial (la que se escribe en una espec):

        SM = (PR_surge/ṁ_surge) / (PR_wl/ṁ_wl) − 1

    evaluada a velocidad de diseño. Se reporta también el margen en
    GASTO puro (ṁ_wl − ṁ_surge)/ṁ_wl, que es el número que la gente usa
    para hablar. La práctica pide 15-25%; SM_FLOW_MIN es el mínimo que el
    optimizador exige como restricción dura.
    """
    theta = pad_theta(theta)
    rec = record if record is not None else _meanline(theta)
    frozen = dict(rec["frozen"]) if rec.get("frozen") else \
        dict(_meanline(theta)["frozen"])
    frozen.setdefault("phi_d", float(theta[3]))
    RPM_d, mdot_d = float(theta[1]), float(theta[12])
    PR_d, eta_d = rec["PR"], rec["eta_poly"]
    gam_d = gamma_air(float(theta[10]))

    m_choke = _choke_mdot(theta, RPM_d, frozen, mdot_d, vsv_deg)
    m_surge, rec_s = _surge_mdot(theta, RPM_d, frozen, mdot_d, m_choke,
                                 vsv_deg)
    if rec_s is None:
        return dict(ok=False, sm_flow=0.0, sm_koch=0.0,
                    mdot_surge=float(m_surge), mdot_choke=float(m_choke),
                    mdot_wl=float(mdot_d), PR_surge=float(PR_d),
                    PR_wl=float(PR_d), stage_stall_first=-1,
                    reason="sin punto de bombeo localizable en el rango")

    # Punto de trabajo: intersección de la línea de trabajo con esta
    # speedline (punto fijo en ṁ; 12 iteraciones bastan y es monótono).
    m_wl = mdot_d
    PR_wl = PR_d
    for _ in range(12):
        r = _meanline(theta, rpm=RPM_d, mdot=m_wl, frozen=frozen,
                      vsv_deg=vsv_deg)
        PR_wl = r["PR"]
        m_new = working_line_mdot(PR_wl, PR_d, mdot_d, eta_d, gam_d)
        m_new = float(np.clip(m_new, m_surge, m_choke))
        if abs(m_new - m_wl) < 1e-4 * mdot_d:
            m_wl = m_new
            break
        m_wl = 0.5 * (m_wl + m_new)

    sm_koch = (rec_s["PR"] / max(m_surge, 1e-9)) / \
        max(PR_wl / max(m_wl, 1e-9), 1e-9) - 1.0
    sm_flow = (m_wl - m_surge) / max(m_wl, 1e-9)
    return dict(ok=True, sm_flow=float(sm_flow), sm_koch=float(sm_koch),
                mdot_surge=float(m_surge), mdot_choke=float(m_choke),
                mdot_wl=float(m_wl), PR_surge=float(rec_s["PR"]),
                PR_wl=float(PR_wl),
                stage_stall_first=int(rec_s["stage_stall_first"]),
                stage_SM=[float(s["SM"]) for s in rec_s["stage_table"]])


def compressor_map(theta: np.ndarray,
                   speed_fracs=(0.7, 0.8, 0.9, 1.0, 1.05),
                   mdot_range=(0.55, 1.25), n_points=17,
                   vsv: str = "none") -> dict:
    """Mapa de operación L0: speedlines PR(ṁ) con marcadores de límite.

    Surge: SM de Koch <= 0 o rama de pendiente positiva a la izquierda del
    pico de PR. Choke: continuidad sin raíz subsónica en alguna estación.
    Da la TENDENCIA; recalibrar contra CFD/banco (docs/VALIDATION.md).

    `vsv`: "none" = geometría fija (default). Para PR ≳ 4 las speedlines
    de baja velocidad con geometría fija son EXTRAPOLACIÓN NO FÍSICA —
    las etapas frontales están en stall profundo; es la razón histórica
    de los estátores variables (J79). Con "auto", cada speedline cierra
    los VSV frontales Δ = VSV_GAIN_DEG·(1−N/Nd) (tope VSV_MAX_DEG, solo
    a N < Nd) — schedule lineal realista de primer orden.
    """
    theta = pad_theta(theta)
    base = _meanline(theta)
    frozen = dict(base["frozen"])
    frozen.setdefault("phi_d", float(theta[3]))
    RPM_d, mdot_d = float(theta[1]), float(theta[12])
    PR_d, eta_d = base["PR"], base["eta_poly"]
    gam_d = gamma_air(float(theta[10]))
    out = dict(design=dict(RPM=RPM_d, mdot=mdot_d, PR=PR_d,
                           eta_poly=eta_d), speedlines=[])
    for fr in speed_fracs:
        rpm = fr * RPM_d
        vsv_deg = 0.0
        if vsv == "auto":
            vsv_deg = min(VSV_GAIN_DEG * max(1.0 - fr, 0.0), VSV_MAX_DEG)
        # Límites REALES de esta speedline: el choke acota el gasto por
        # arriba (la línea se hace vertical, no se desploma) y el bombeo
        # por abajo. Antes se barría un rango fijo 0.55-1.25·ṁ_d y se
        # marcaban banderas, con lo que casi todos los puntos quedaban
        # fuera y el "rango estable" era un único punto.
        m_choke = _choke_mdot(theta, rpm, frozen, mdot_d, vsv_deg)
        m_surge, rec_s = _surge_mdot(theta, rpm, frozen, mdot_d, m_choke,
                                     vsv_deg)
        lo = max(m_surge, mdot_range[0] * mdot_d * 0.5)
        hi = m_choke
        pts = []
        if hi > lo * 1.0005:
            for mm in np.linspace(lo, hi, n_points):
                r = _meanline(theta, rpm=rpm, mdot=float(mm), frozen=frozen,
                              vsv_deg=vsv_deg)
                pts.append(dict(
                    mdot=float(mm), PR=r["PR"], eta_poly=r["eta_poly"],
                    min_SM=r["min_SM"], M_rel_tip1=r["M_rel_tip1"],
                    choke=bool(r["choke_flag"]),
                    stall=bool(r["stall_flag"]),
                    stage_stall_first=int(r["stage_stall_first"]),
                    unstable=False))
        # Punto de trabajo de esta velocidad (intersección con la línea de
        # trabajo que pasa por el diseño)
        wl = None
        if pts:
            m_wl = float(np.clip(mdot_d * fr, lo, hi))
            for _ in range(12):
                r = _meanline(theta, rpm=rpm, mdot=m_wl, frozen=frozen,
                              vsv_deg=vsv_deg)
                m_new = float(np.clip(
                    working_line_mdot(r["PR"], PR_d, mdot_d, eta_d, gam_d),
                    lo, hi))
                if abs(m_new - m_wl) < 1e-4 * mdot_d:
                    m_wl = m_new
                    break
                m_wl = 0.5 * (m_wl + m_new)
            r = _meanline(theta, rpm=rpm, mdot=m_wl, frozen=frozen,
                          vsv_deg=vsv_deg)
            sm_flow = (m_wl - m_surge) / max(m_wl, 1e-9)
            sm_koch = ((rec_s["PR"] / max(m_surge, 1e-9)) /
                       max(r["PR"] / max(m_wl, 1e-9), 1e-9) - 1.0) \
                if rec_s else 0.0
            wl = dict(mdot=float(m_wl), PR=float(r["PR"]),
                      eta_poly=float(r["eta_poly"]),
                      sm_flow=float(sm_flow), sm_koch=float(sm_koch),
                      stage_stall_first=int(r["stage_stall_first"]),
                      stage_SM=[float(s["SM"]) for s in r["stage_table"]])
        out["speedlines"].append(dict(
            rpm=float(rpm), speed_frac=float(fr), vsv_deg=float(vsv_deg),
            points=pts,
            mdot_surge=float(m_surge) if pts else None,
            mdot_choke=float(m_choke) if pts else None,
            PR_surge=float(rec_s["PR"]) if rec_s else None,
            working_point=wl))
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


_L1_MARK = "<<<PHYAC_L1>>>"

_L1_CHILD_SRC = (
    "import json,sys;"
    "sys.path.insert(0, sys.argv[2]);"
    "import numpy as np, physics_core as pc;"
    "r = pc._scm_solve_direct(np.array(json.loads(sys.argv[1]), "
    "dtype=float));"
    "sys.stdout.write({mark!r} + json.dumps(r, default=float))"
).format(mark=_L1_MARK)


def _scm_solve(theta: np.ndarray) -> dict | None:
    """Spool SCM axial con timeout duro (proceso hijo).

    El bucle de convergencia de gasto de TD3 puede colgarse (spike M7);
    el hijo se termina a los L1_TIMEOUT_S segundos y el punto degrada a L0
    etiquetado.

    Se lanza con `python -c`, NO con multiprocessing. En Windows (método
    `spawn`) el hijo de multiprocessing reimporta el módulo `__main__` del
    padre para poder deserializar el target: cualquier script que llame a
    la física sin envolver su cuerpo en `if __name__ == "__main__"` se
    reejecuta entero dentro del hijo, agota el timeout y degrada a L0 en
    SILENCIO. Así se descubrió: la suite de verificación —que no tiene
    guard, porque es un script— creía correr en L1 desde siempre. El
    proceso hijo explícito no importa nada del padre y elimina la
    dependencia de cómo esté escrito el que llama.
    """
    if not _try_load_turbodesign():
        return None
    import subprocess
    root = os.path.dirname(os.path.abspath(__file__))
    try:
        out = subprocess.run(
            [sys.executable, "-c", _L1_CHILD_SRC,
             json.dumps(np.asarray(theta, dtype=float).tolist()), root],
            capture_output=True, text=True, timeout=L1_TIMEOUT_S,
            cwd=root)
    except Exception:      # TimeoutExpired incluido: el hijo ya está muerto
        return None
    if out.returncode != 0 or _L1_MARK not in out.stdout:
        if os.environ.get("PHYAC_DEBUG_L1"):
            print(f"[L1] hijo rc={out.returncode}\n{out.stderr[-2000:]}")
        return None
    try:
        rec = json.loads(out.stdout.split(_L1_MARK, 1)[1])
    except Exception:
        return None
    return rec if isinstance(rec, dict) else None


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
             use_cache: bool = True, calibrate: bool = True,
             per_stage: dict | None = None,
             with_surge_margin: bool = False) -> dict:
    """Evalúa un punto de diseño físico (15-D; acepta el θ legacy de 13
    vía pad_theta — pendientes 0, bit-exacto con el espacio pre-fase-8).

    L1 corre el SCM axial y, si diverge o no está instalado, degrada a L0
    etiquetando la fuente. Siempre devuelve el record completo del
    meanline (stage_table, restricciones, banderas) — el SCM solo
    sobreescribe PR/eta cuando converge.

    `per_stage` (modo experto, fase 8b): dict {"phi"/"psi"/"Rx": [...]}
    que sobreescribe las distribuciones por etapa (ver _meanline). Entra
    a la clave de caché y fuerza la ruta L0 pura (la transformación del
    spool L1 no conoce el override).
    """
    theta = pad_theta(theta)
    extra = json.dumps(per_stage, sort_keys=True) if per_stage else ""
    key = _theta_key(theta, int(fidelity), extra)
    if use_cache and key in _CACHE:
        rec = dict(_CACHE[key])
    else:
        rec = _meanline(theta, per_stage=per_stage)
        rec.update({n: float(v) for n, v in zip(VAR_NAMES, theta)})
        rec["n_stages"] = int(round(theta[0]))
        rec["fidelity"] = int(fidelity)
        if fidelity >= Fidelity.L1 and rec["feasible"] and not per_stage:
            scm = _scm_solve(theta)
            if scm:
                rec.update(scm)
            else:
                rec["source"] = rec["source"] + "(L1_unavailable_or_diverged)"
        if use_cache:
            _cache_store(key, rec)
    if with_surge_margin and "sm_flow" not in rec:
        # ~80 llamadas al meanline (≈80 ms). Se pide SOLO en las
        # evaluaciones físicas del orquestador, nunca en el lazo interno
        # del NSGA-II (que corre 5 760 individuos por ronda).
        try:
            sm = surge_margin(theta, rec if rec.get("frozen") else None)
            rec["sm_flow"] = float(sm["sm_flow"])
            rec["sm_koch_wl"] = float(sm["sm_koch"])
            rec["surge"] = sm
        except Exception:
            rec["sm_flow"] = 0.0
            rec["sm_koch_wl"] = 0.0
    if calibrate:
        rec = CALIBRATION.apply(rec)
    return rec


def evaluate_batch(thetas: np.ndarray, fidelity: Fidelity = Fidelity.L1,
                   n_workers: int = 1,
                   with_surge_margin: bool = False) -> list[dict]:
    """Evalúa un batch. n_workers>1 paraleliza por procesos (útil en L1)."""
    thetas = np.atleast_2d(np.asarray(thetas, dtype=float))
    if n_workers <= 1 or fidelity == Fidelity.L0:
        return [evaluate(t, fidelity, with_surge_margin=with_surge_margin)
                for t in thetas]
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        return list(ex.map(_eval_worker,
                           [(t.tolist(), int(fidelity), with_surge_margin)
                            for t in thetas]))


def _eval_worker(args):
    theta, fid, sm = args
    return evaluate(np.array(theta), Fidelity(fid), with_surge_margin=sm)


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
