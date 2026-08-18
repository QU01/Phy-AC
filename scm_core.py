"""
QUASAR Phy-AC · scm_core.py — FIDELIDAD L1: through-flow por curvatura de
líneas de corriente (SCM)
=========================================================================
Segundo peldaño REAL de la escalera de fidelidad. Resuelve el flujo
axisimétrico sobre N líneas de corriente (5-11) con la ecuación de
equilibrio radial COMPLETA —incluido el término de curvatura meridional—
y continuidad por tubo de corriente, en vez de una sola línea media.

Por qué existe
--------------
Hasta la fase 11 «L1» era `turbo-design` 1.4.2 con tres parches, corriendo
obligatoriamente con UNA línea de corriente porque con más el ODE de
equilibrio radial de esa librería divergía. Es decir: L1 también era
meanline. Con un solo peldaño real, el ensemble residual de la capa 2 no
tiene residual que aprender (se cortocircuitaba en L0) y la calibración
afín L2 era una API esperando usuario. Todo lo que se construye encima
—geometría, estructura, optimizador— se apoyaba en un 1-D con
correlaciones.

Qué resuelve
------------
El álabe es un objeto FIJO EN EL ESPACIO: sus ángulos metálicos β₂(r) y
α_salida(r) salen de la ley de torbellino de diseño y se congelan. El
solver busca el campo (Cm, Cu, ρ) que satisface a la vez:

  1. equilibrio radial en cada estación,
  2. continuidad global y por tubo de corriente,
  3. el ángulo que impone el álabe en cada borde de salida,
  4. Euler POR LÍNEA DE CORRIENTE — el trabajo deja de ser uniforme en el
     span, que es justo lo que el meanline no puede ver.

De ahí sale un resultado que NO es el del meanline: el trabajo se
redistribuye radialmente, la densidad y el bloqueo varían con el radio, y
las pérdidas se calculan donde ocurren (perfil en todo el span, choque
solo donde el Mach relativo lo justifica, holgura en la banda de punta,
capa límite de pared en las bandas de pared).

Ecuación de equilibrio radial
-----------------------------
Sobre una cuasi-ortogonal radial, con h₀ la entalpía de parada rotálpica
del marco absoluto, s la entropía y r_c el radio de curvatura meridional:

    Cm·∂Cm/∂r = ∂h₀/∂r − T·∂s/∂r − (Cu/r)·∂(r·Cu)/∂r − Cm²·cos(γ)/r_c

que sale de combinar dh = T ds + dp/ρ con el equilibrio normal
(1/ρ)∂p/∂r = Cu²/r + (Cm²/r_c)·cos γ. El signo: 1/r_c > 0 cuando el
centro de curvatura está a MENOR radio (línea cóncava hacia el eje), es
decir 1/r_c = −r''/(1+r'²)^{3/2}, y cos γ = 1/√(1+r'²).

Con curvatura nula y Cu ∝ r^n se recupera exactamente la forma cerrada de
`physics_core.vortex_cx` de la fase 9.1 — el test T21 lo comprueba, y esa
compatibilidad es lo que hace que este módulo sea una EXTENSIÓN del
modelo anterior y no un modelo paralelo.

Numérica
--------
Marcha estación a estación (borde de ataque y de salida de cada fila),
con un lazo exterior que recoloca las líneas de corriente por fracción de
gasto y actualiza la curvatura. El nivel de Cm en cada estación sale de
resolver la continuidad global por bisección sobre Cm en el cubo; el
PERFIL de Cm lo fija el equilibrio radial. El término de curvatura va
retrasado una iteración (práctica estándar; es el que desestabiliza si se
mete implícito).

Sin dependencias externas: solo numpy y el resto de physics_core.
"""

from __future__ import annotations

import math
import os

import numpy as np

import physics_core as pc

# ---------------------------------------------------------------------------
# Parámetros del solver
# ---------------------------------------------------------------------------
N_STREAMLINES = int(os.environ.get("PHYAC_SCM_STREAMLINES", "9"))
#   Líneas de corriente hub→punta, repartidas por FRACCIÓN DE GASTO (no
#   por radio): así los tubos de corriente llevan el mismo caudal y la
#   resolución se concentra sola donde la densidad de gasto es alta.
#   5 es el mínimo con el que la derivada radial tiene sentido; por encima
#   de 11 el coste sube sin que el resultado se mueva (test T21).

SCM_MAX_ITER = int(os.environ.get("PHYAC_SCM_ITER", "60"))
SCM_TOL = 1e-5                 # convergencia en radio relativo
SCM_RELAX_R = 0.35             # relajación al recolocar líneas
SCM_RELAX_CURV = 0.30          # relajación de la curvatura meridional
CURV_MAX = 8.0                 # |1/r_c| máx [1/m] — el término de
#                                curvatura es pequeño en un axial y un
#                                pico numérico en r'' lo puede disparar
CM_FLOOR_FRAC = 0.25           # Cm mínimo como fracción del Cm de la línea
#                                media (evita raíces negativas con leyes
#                                de torbellino agresivas)
CM_SPREAD_MIN = 0.30           # límites del PERFIL de Cm respecto a su
CM_SPREAD_MAX = 2.40           # media. Es un limitador NUMÉRICO, no un
#                                modelo: el cierre del álabe en la punta de
#                                un fan es rígido (tanβ₂ ≈ 2, así que 10
#                                m/s de error en Cm son 20 en Cu y ~7 kJ/kg
#                                de trabajo) y sin él una línea se colapsa.
#                                Si sigue ACTIVO al converger, el resultado
#                                no vale: `solve` lanza SCMDiverged en vez
#                                de devolver un campo amasado por el
#                                limitador.
BLADE_RELAX = 0.35             # relajación del cierre Cm↔Cu del álabe
PR_WINDOW = 0.15               # |PR_L1/PR_L0 − 1| máximo admisible.
#   NO es una afirmación de física: es la guarda contra un problema
#   ESTRUCTURAL de acoplar los dos niveles. El annulus lo dimensiona L0
#   con SU Cx uniforme y SU densidad media; L1 resuelve un perfil, así
#   que su Cm medio no coincide exactamente, el álabe —de ángulo fijo—
#   convierte esa diferencia en trabajo, el trabajo cambia la densidad y
#   la siguiente estación arranca con más diferencia. En 1-4 etapas es
#   ruido; en 7-8 el lazo se compone hasta ±27% de PR, y eso no es
#   fidelidad, es deriva. Medido en validation/BENCH_SCM.md.
#   La cura de fondo es que el annulus salga del MISMO solver que lo usa
#   (o que L0 lo dimensione con el perfil de L1); mientras tanto, el
#   punto se rechaza y se degrada a L0 etiquetado en vez de devolver un
#   número que nadie debería usar.
TIP_BAND_FRAC = 0.25           # fracción de span sobre la que se reparte
#                                la pérdida de holgura de punta
WALL_BAND_FRAC = 0.30          # ídem para la capa límite de pared
SPAN_MIX_PER_ROW = 0.20        # mezcla radial EFECTIVA entre tubos de
#   corriente vecinos, aplicada una vez por fila a s y h₀. Sin ella el
#   método no tiene mecanismo que reparta lo que la turbulencia y los
#   flujos secundarios reparten en la máquina real (Adkins & Smith 1982,
#   Gallimore & Cumpsty 1986): las bandas de pared acumulan la entropía
#   de TODAS las filas en las mismas líneas, el gradiente T·∂s/∂r
#   distorsiona el equilibrio radial, el álabe de ángulo fijo convierte
#   la distorsión en trabajo no uniforme y el lazo se compone. Medido en
#   el E³ (10 etapas) antes de esta constante: ds de pared 105-135
#   J/kg·K contra 20 en medio span, h₀ del cubo 212 kJ/kg contra 117 en
#   punta y CONTRA-remolino de -224 m/s en una salida de rotor — hasta
#   bloquear la marcha en la etapa 5. En monoetapa apenas actúa (una o
#   dos aplicaciones); su efecto medido en los anclajes L1 es <0.1% de
#   PR. El intercambio es SIMÉTRICO entre tubos de igual gasto, así que
#   conserva la media másica de h₀ y s exactamente (la generación de
#   entropía del propio mezclado es de segundo orden y se desprecia).


class SCMDiverged(RuntimeError):
    """El solver no convergió. Nunca se traga: el que llama decide si
    degrada a L0 y lo ETIQUETA."""


# ---------------------------------------------------------------------------
# Utilidades radiales
# ---------------------------------------------------------------------------
def _ddr(f: np.ndarray, r: np.ndarray) -> np.ndarray:
    """Derivada radial por diferencias centradas sobre malla no uniforme."""
    d = np.empty_like(f)
    d[1:-1] = ((f[2:] - f[:-2]) / np.maximum(r[2:] - r[:-2], 1e-12))
    d[0] = (f[1] - f[0]) / max(r[1] - r[0], 1e-12)
    d[-1] = (f[-1] - f[-2]) / max(r[-1] - r[-2], 1e-12)
    return d


def _curvature(z: np.ndarray, r: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(1/r_c, cos γ) de una línea de corriente muestreada en (z, r).

    1/r_c = −r''/(1+r'²)^{3/2} con la convención del encabezado: positivo
    cuando el centro de curvatura está a MENOR radio.
    """
    n = len(z)
    rp = np.zeros(n)
    rpp = np.zeros(n)
    if n >= 3:
        rp[1:-1] = (r[2:] - r[:-2]) / np.maximum(z[2:] - z[:-2], 1e-12)
        rp[0] = (r[1] - r[0]) / max(z[1] - z[0], 1e-12)
        rp[-1] = (r[-1] - r[-2]) / max(z[-1] - z[-2], 1e-12)
        for k in range(1, n - 1):
            h0 = max(z[k] - z[k - 1], 1e-9)
            h1 = max(z[k + 1] - z[k], 1e-9)
            rpp[k] = 2.0 * (h0 * r[k + 1] - (h0 + h1) * r[k] + h1 * r[k - 1]) \
                / (h0 * h1 * (h0 + h1))
    cos_g = 1.0 / np.sqrt(1.0 + rp ** 2)
    kappa = -rpp / np.power(1.0 + rp ** 2, 1.5)
    return np.clip(kappa, -CURV_MAX, CURV_MAX), cos_g


# ---------------------------------------------------------------------------
# Estado termodinámico por línea de corriente
# ---------------------------------------------------------------------------
def _cp_v(T: np.ndarray) -> np.ndarray:
    """cp(T) VECTORIZADO. Misma ley que physics_core.cp_air."""
    if not pc.GAS_VARIABLE:
        return np.full_like(T, pc.CP)
    t = np.clip(T, 200.0, 1400.0)
    return pc.CP_A + pc.CP_B * t + pc.CP_C * t * t


def _h_v(T: np.ndarray) -> np.ndarray:
    if not pc.GAS_VARIABLE:
        return pc.CP * T
    t = np.clip(T, 200.0, 1400.0)
    return pc.CP_A * t + 0.5 * pc.CP_B * t * t + (pc.CP_C / 3.0) * t ** 3


def _phi_v(T: np.ndarray) -> np.ndarray:
    t = np.clip(T, 200.0, 1400.0)
    if not pc.GAS_VARIABLE:
        return pc.CP * np.log(t)
    return pc.CP_A * np.log(t) + pc.CP_B * t + 0.5 * pc.CP_C * t * t


def _T_from_h_v(h: np.ndarray, T_guess: np.ndarray) -> np.ndarray:
    """Inversa de h(T) por Newton VECTORIZADO.

    Es el cuello de botella real del solver: la bisección de continuidad
    la llama ~90 veces por estación y por iteración interna. En bucle
    Python sobre las líneas de corriente el SCM tardaba 12 s por máquina;
    vectorizado baja a menos de 1 s sin cambiar un dígito del resultado.
    """
    if not pc.GAS_VARIABLE:
        return h / pc.CP
    # Newton SIN clip dentro del lazo: sobre vectores de 9 elementos el
    # overhead de np.clip domina el coste (1.9 M llamadas por máquina en
    # el perfilado). El rango se acota una sola vez al salir; el suelo de
    # entalpía de _static_state ya impide entrar fuera de rango.
    a, b, c = pc.CP_A, pc.CP_B, pc.CP_C
    T = T_guess
    for _ in range(5):
        T = T - (a * T + 0.5 * b * T * T + (c / 3.0) * T * T * T - h)             / (a + b * T + c * T * T)
    return np.clip(T, 200.0, 1400.0)


def _static_state(h0: np.ndarray, c2: np.ndarray, ds: np.ndarray,
                  T0_in: float, P0_in: float
                  ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(T, p, ρ) estáticos desde entalpía de parada, energía cinética y
    la entropía GENERADA respecto a la entrada.

    Gas caloríficamente imperfecto, igual que el meanline: la presión sale
    de la función de entropía, s − s_in = φ(T) − φ(T_in) − R·ln(p/p_in).
    """
    phi_in = pc.phi_air(T0_in)
    # Suelo de temperatura estática: con Cm de prueba grandes la bisección
    # explora estados donde h₀ − c²/2 se vuelve negativo y la inversa de
    # h devuelve basura. Medio T₀ es un suelo holgado (M ≈ 2.2).
    h_st = np.maximum(h0 - 0.5 * c2, pc.h_air(0.5 * T0_in))
    T = _T_from_h_v(h_st, np.full(len(h0), max(T0_in, 200.0)))
    ln_p = np.clip((_phi_v(T) - phi_in - ds) / pc.RGAS, -20.0, 20.0)
    p = P0_in * np.exp(ln_p)
    rho = p / (pc.RGAS * np.maximum(T, 1.0))
    return T, p, rho


def _total_pressure(T: np.ndarray, p: np.ndarray,
                    h0: np.ndarray) -> np.ndarray:
    """P₀ por línea de corriente, isentrópico desde el estático local."""
    T0 = _T_from_h_v(h0, T + 50.0)
    return p * np.exp((_phi_v(T0) - _phi_v(T)) / pc.RGAS)


# ---------------------------------------------------------------------------
# Geometría de estaciones
# ---------------------------------------------------------------------------
def _stations(record: dict) -> list[dict]:
    """Estaciones de cálculo desde el record meanline.

    2n+1 estaciones: borde de ataque del rotor i (= salida del estátor
    i−1), borde de salida del rotor i, borde de salida del estátor i.
    La geometría (r_hub, r_tip, z) es la MISMA que emite la capa 5a, así
    que el SCM resuelve el flujo en la máquina que se fabrica.
    """
    st = record["stage_table"]
    out: list[dict] = []
    z = 0.0
    for i, s in enumerate(st):
        r_m = s["r_mean_mm"] / 1000.0
        c_r = s["chord_rotor_mm"] / 1000.0
        c_s = s["chord_stator_mm"] / 1000.0
        h_in = s["A_in_m2"] / (2.0 * math.pi * r_m)
        h_re = s["A_rotor_exit_m2"] / (2.0 * math.pi * r_m)
        if i == 0:
            out.append(dict(kind="le", stage=i, z=z, r_m=r_m,
                            r_hub=r_m - 0.5 * h_in, r_tip=r_m + 0.5 * h_in,
                            area=s["A_in_m2"], kb=s["blockage_in"]))
        # el rotor ocupa su cuerda axial (proyección por el stagger medio)
        ax_r = c_r * math.cos(math.radians(
            0.5 * (s["beta1_deg"] + s["beta2_deg"])))
        z += max(ax_r, 1e-4)
        out.append(dict(kind="rte", stage=i, z=z, r_m=r_m,
                        r_hub=r_m - 0.5 * h_re, r_tip=r_m + 0.5 * h_re,
                        area=s["A_rotor_exit_m2"], kb=s["blockage_in"]))
        ax_s = c_s * math.cos(math.radians(
            0.5 * (s["alpha1_deg"] + s["alpha2_deg"])))
        z += max(ax_s, 1e-4) * (1.0 + pc.ROW_GAP_FRACTION)
        nxt = st[i + 1] if i + 1 < len(st) else None
        if nxt is not None:
            r_m2 = nxt["r_mean_mm"] / 1000.0
            h2 = nxt["A_in_m2"] / (2.0 * math.pi * r_m2)
            out.append(dict(kind="le", stage=i + 1, z=z, r_m=r_m2,
                            r_hub=r_m2 - 0.5 * h2, r_tip=r_m2 + 0.5 * h2,
                            area=nxt["A_in_m2"], kb=nxt["blockage_in"]))
        else:
            # Salida de la última etapa. El area la da la marcha del
            # meanline (frozen["areas"][-1]): con el gas ya comprimido, el
            # annulus que mantiene el Cx de diseño es MENOR que el de
            # entrada de esa etapa, y usar A_in dejaba una Cm de salida
            # 25% baja que falseaba la energía cinética del promediado.
            fr = record.get("frozen") or {}
            a_ex = float((fr.get("areas") or [s["A_in_m2"]])[-1])
            h2 = a_ex / (2.0 * math.pi * r_m)
            out.append(dict(kind="ste", stage=i, z=z, r_m=r_m,
                            r_hub=r_m - 0.5 * h2, r_tip=r_m + 0.5 * h2,
                            area=a_ex, kb=s["blockage_out"]))
    return out


def _blade_angles(record: dict, vortex_n: float | None = None) -> list[dict]:
    """Ángulos METÁLICOS de cada fila como función del radio.

    Se congelan desde la ley de torbellino de diseño: el álabe es un
    objeto fijo en el espacio y el solver NO puede moverlo. Es la
    diferencia entre resolver el flujo y volver a postular el triángulo.
    """
    n = pc.VORTEX_N if vortex_n is None else vortex_n
    st = record["stage_table"]
    rows = []
    for i, s in enumerate(st):
        r_m = s["r_mean_mm"] / 1000.0
        Cx = s["Cx"]
        U = s["U_m"]
        Cu1 = Cx * math.tan(math.radians(s["alpha1_deg"]))
        Cu2 = Cx * math.tan(math.radians(s["alpha2_deg"]))
        # El estátor de la etapa i entrega el pre-swirl que necesita la
        # etapa i+1, NO el de la suya: con pendientes por etapa (psi_slope,
        # Rx_slope) esos dos ángulos son distintos, y usar el propio hacía
        # que todas las etapas vieran la misma entrada y el trabajo no
        # cayera hacia atrás como manda el diseño (medido: +13.6% de
        # trabajo total). En la última etapa el remolino residual lo quita
        # el OGV, así que el estátor entrega el alpha1 repetitivo.
        s_next = st[min(i + 1, len(st) - 1)]
        Cu1_next = s_next["Cx"] * math.tan(math.radians(s_next["alpha1_deg"]))
        rows.append(dict(stage=i, r_m=r_m, Cx=Cx, U_m=U, Cu1_m=Cu1,
                         Cu2_m=Cu2, n=n,
                         Cu1_next_m=Cu1_next,
                         Cx_next=s_next["Cx"],
                         r_m_next=s_next["r_mean_mm"] / 1000.0,
                         sigma_r=s["sigma_rotor_actual"],
                         sigma_s=s["sigma_stator_actual"],
                         c_r=s["chord_rotor_mm"] / 1000.0,
                         c_s=s["chord_stator_mm"] / 1000.0,
                         h_blade=s["h_blade_mm"] / 1000.0,
                         camber_r=s["beta1_deg"] - s["beta2_deg"],
                         camber_s=s["alpha2_deg"] - s["alpha1_deg"]))
    return rows


def _beta2_metal(row: dict, r: np.ndarray, omega: float) -> np.ndarray:
    """β₂(r) del rotor: el ángulo relativo de salida que talla el álabe."""
    out = np.empty(len(r))
    for j, rj in enumerate(r):
        cu2 = pc.vortex_cu(row["Cu2_m"], rj, row["r_m"], row["n"])
        cx = pc.vortex_cx(row["Cx"], row["Cu2_m"], rj, row["r_m"], row["n"])
        out[j] = math.atan2(omega * rj - cu2, max(cx, 1e-3))
    return out


def _alpha_out_metal(row: dict, r: np.ndarray) -> np.ndarray:
    """α_salida(r) del estátor: devuelve el flujo al pre-swirl de la etapa
    SIGUIENTE (o al alpha1 repetitivo en la última, donde el remolino
    residual lo quita el OGV)."""
    cu_t, cx_t, r_t = row["Cu1_next_m"], row["Cx_next"], row["r_m_next"]
    out = np.empty(len(r))
    for j, rj in enumerate(r):
        cu1 = pc.vortex_cu(cu_t, rj, r_t, row["n"])
        cx = pc.vortex_cx(cx_t, cu_t, rj, r_t, row["n"])
        out[j] = math.atan2(cu1, max(cx, 1e-3))
    return out


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------
def _span_mix(f: np.ndarray,
              alpha: float = SPAN_MIX_PER_ROW) -> np.ndarray:
    """Una pasada de mezcla radial entre tubos de corriente vecinos.

    Difusión explícita en el índice de tubo (los tubos llevan el mismo
    gasto, así que el intercambio simétrico conserva la media másica
    exactamente), con frontera de Neumann: la pared no aporta ni quita.
    α ≤ 0.5 por estabilidad del esquema explícito.
    """
    g = f.copy()
    g[1:-1] += alpha * (f[2:] - 2.0 * f[1:-1] + f[:-2])
    g[0] += alpha * (f[1] - f[0])
    g[-1] += alpha * (f[-2] - f[-1])
    return g


def _solve_station_cm(r: np.ndarray, cu: np.ndarray, h0: np.ndarray,
                      ds: np.ndarray, T: np.ndarray, cm_prev: np.ndarray,
                      kappa: np.ndarray, cos_g: np.ndarray,
                      mdot: float, kb: float, T0_in: float, P0_in: float,
                      cm_ref_guess: float) -> tuple[np.ndarray, np.ndarray]:
    """Perfil de Cm(r) en una estación: equilibrio radial + continuidad.

    El equilibrio radial fija la FORMA (Cm² relativo al cubo) y la
    continuidad global el NIVEL. Se resuelve el nivel por bisección: el
    gasto crece de forma monótona con Cm mientras el flujo sea subsónico,
    así que la bisección es robusta donde Newton no lo es.
    """
    rcu = r * cu
    # F(r) = dh0/dr − T·ds/dr − (Cu/r)·d(rCu)/dr
    F = _ddr(h0, r) - T * _ddr(ds, r) - (cu / np.maximum(r, 1e-9)) * _ddr(rcu, r)
    # término de curvatura, retrasado una iteración
    G = F - cm_prev ** 2 * kappa * cos_g
    # Δ_j = Cm²_j − Cm²_0 por trapecio
    dlt = np.zeros(len(r))
    for j in range(1, len(r)):
        dlt[j] = dlt[j - 1] + (r[j] - r[j - 1]) * (G[j] + G[j - 1])

    trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

    limiter = {"hit": False}

    def mass(cm_hub: float):
        cm2 = cm_hub ** 2 + dlt
        floor = (CM_FLOOR_FRAC * max(cm_hub, 1.0)) ** 2
        cm = np.sqrt(np.maximum(cm2, floor))
        mean = max(float(np.mean(cm)), 1e-6)
        lo_c, hi_c = CM_SPREAD_MIN * mean, CM_SPREAD_MAX * mean
        if float(np.min(cm)) < lo_c or float(np.max(cm)) > hi_c:
            limiter["hit"] = True
            cm = np.clip(cm, lo_c, hi_c)
        c2 = cm ** 2 + cu ** 2
        _T, _p, rho = _static_state(h0, c2, ds, T0_in, P0_in)
        a = np.sqrt(pc.GAMMA * pc.RGAS * np.maximum(_T, 1.0))
        # Mach MERIDIONAL, no absoluto. Dentro de esta bisección Cu está
        # FIJO, así que d(ρ·Cm)/dCm = ρ·(1 − Cm²/a²): el gasto satura
        # cuando el Mach MERIDIONAL llega a 1, sin que importe el remolino.
        # Con la guarda sobre el Mach absoluto (el bug que esto arregla),
        # el cubo de un vórtice libre con HTR bajo —Cu_hub ∝ r_m/r_hub—
        # tocaba M_abs 0.98 con Cm aún pequeño y la estación se declaraba
        # BLOQUEADA al ~55% de su capacidad real: el R67 (HTR 0.375) y el
        # E³ enteros se perdían por esto. Flujo localmente supersónico en
        # ABSOLUTO con Cm subsónico es un estado normal detrás del cubo de
        # un rotor transónico; el estátor de detrás paga su choque en el
        # modelo de pérdidas, no aquí.
        mach = float(np.max(cm / a))
        f = rho * cm * cos_g * 2.0 * math.pi * r * kb
        m = 0.5 * float(np.sum((f[1:] + f[:-1]) * (r[1:] - r[:-1])))
        return m, cm, _T, mach

    # El gasto SOLO crece con Cm en la rama subsónica; pasado el punto
    # sónico vuelve a caer, y una bisección ciega se va a la rama
    # supersónica y devuelve un estado sin sentido (medido: PR = 0.19 con
    # T₀ de salida por DEBAJO de la de entrada). Se acota el intervalo por
    # el límite sónico antes de bisecar.
    lo = 1.0
    hi = max(0.5 * cm_ref_guess, 10.0)
    m_hi, _, _, mach_hi = mass(hi)
    for _ in range(90):
        if m_hi >= mdot or mach_hi >= 0.98:
            break
        hi *= 1.25
        m_hi, _, _, mach_hi = mass(hi)
    if m_hi < mdot:
        raise SCMDiverged("estación BLOQUEADA: ni en el límite sónico pasa "
                          f"el gasto ({m_hi:.3f} < {mdot:.3f} kg/s)")
    # 34 bisecciones dan 6e-11 de precisión relativa en el intervalo — de
    # sobra, y la bisección es el 60% del coste del solver (se llama una
    # vez por iteración del cierre del álabe, por fila y por pasada).
    for _ in range(34):
        mid = 0.5 * (lo + hi)
        m_mid, cm, T_new, _ = mass(mid)
        if abs(m_mid - mdot) < 1e-9 * max(mdot, 1.0):
            break
        if m_mid < mdot:
            lo = mid
        else:
            hi = mid
    _, cm, T_new, _ = mass(0.5 * (lo + hi))
    return cm, T_new, limiter["hit"]


def _reposition(r_hub: float, r_tip: float, r: np.ndarray, rho: np.ndarray,
                cm: np.ndarray, cos_g: np.ndarray, kb: float,
                frac: np.ndarray) -> np.ndarray:
    """Recoloca las líneas de corriente para que cada tubo lleve su
    fracción de gasto. Es la condición que cierra el problema: sin ella el
    solver resolvería el equilibrio radial sobre una malla arbitraria."""
    f = rho * cm * cos_g * 2.0 * math.pi * r * kb
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1])
                                           * np.diff(r))])
    tot = cum[-1]
    if tot <= 0:
        raise SCMDiverged("gasto acumulado nulo al recolocar líneas")
    r_new = np.interp(frac * tot, cum, r)
    r_new[0], r_new[-1] = r_hub, r_tip      # las paredes no se mueven
    return r_new


def solve(theta: np.ndarray, record: dict,
          n_streamlines: int | None = None) -> dict:
    """Resuelve el through-flow L1 sobre el diseño que describe `record`.

    Devuelve un dict con PR, eta_poly, T0_out y los perfiles radiales.
    Lanza SCMDiverged si no converge — el que llama decide si degrada.
    """
    theta = pc.pad_theta(np.asarray(theta, dtype=float))
    T0_in, P0_in, mdot = float(theta[10]), float(theta[11]), float(theta[12])
    omega = float(theta[1]) * 2.0 * math.pi / 60.0
    n_sl = int(n_streamlines or N_STREAMLINES)
    if n_sl < 5:
        raise ValueError("el equilibrio radial necesita ≥ 5 líneas")

    stations = _stations(record)
    rows = _blade_angles(record)
    n_stat = len(stations)
    frac = np.linspace(0.0, 1.0, n_sl)

    # --- inicialización: líneas equiespaciadas en ÁREA en cada estación ---
    R = np.zeros((n_stat, n_sl))
    for k, stn in enumerate(stations):
        a = np.linspace(0.0, 1.0, n_sl)
        R[k] = np.sqrt(stn["r_hub"] ** 2
                       + a * (stn["r_tip"] ** 2 - stn["r_hub"] ** 2))
    # Cm inicial = el Cx de DISEÑO de la etapa de cada estación. La
    # versión anterior era ṁ/(A·kb·1.2), con 1.2 = densidad AMBIENTE:
    # válida en la etapa 0 y un disparate creciente hacia atrás (en la
    # etapa 4 del E³, ρ real ≈ 4.7 → Cm inicial 552 m/s, meridional
    # supersónico; el cierre del álabe arrancaba en Cu = ωr − 552·tanβ₂ ≈
    # −200 m/s y la primera bisección declaraba la estación bloqueada).
    # Era LA razón de que la cobertura del banco cayera con el número de
    # etapas: no un límite del acoplamiento L0↔L1, un arranque en frío
    # con la densidad equivocada. El meanline ya calculó el Cx de cada
    # etapa con su densidad real — se arranca de ahí.
    CM = np.zeros((n_stat, n_sl))
    st_tab = record["stage_table"]
    for k, stn in enumerate(stations):
        CM[k] = float(st_tab[min(stn["stage"], len(st_tab) - 1)]["Cx"])
    KAPPA = np.zeros((n_stat, n_sl))
    COSG = np.ones((n_stat, n_sl))
    z_stat = np.array([s["z"] for s in stations])

    CU = np.zeros((n_stat, n_sl))
    H0 = np.zeros((n_stat, n_sl))
    DS = np.zeros((n_stat, n_sl))
    TS = np.full((n_stat, n_sl), T0_in)
    rows_out: list[dict] = []

    h0_in = pc.h_air(T0_in)
    converged = False
    lim_hits = 0
    for it in range(SCM_MAX_ITER):
        R_old = R.copy()
        rows_out = []
        lim_hits = 0

        # ---- estación 0: entrada, con el pre-swirl del IGV ---------------
        row0 = rows[0]
        r0 = R[0]
        cu0 = np.array([pc.vortex_cu(row0["Cu1_m"], rj, row0["r_m"],
                                     row0["n"]) for rj in r0])
        CU[0] = cu0
        H0[0] = h0_in
        DS[0] = 0.0
        CM[0], TS[0], _lim = _solve_station_cm(
            r0, cu0, H0[0], DS[0], TS[0], CM[0], KAPPA[0], COSG[0],
            mdot, stations[0]["kb"], T0_in, P0_in, float(np.mean(CM[0])))

        # ---- marcha por filas -------------------------------------------
        for i, row in enumerate(rows):
            k_le, k_te, k_se = 2 * i, 2 * i + 1, 2 * i + 2

            # ===== ROTOR: LE → TE ========================================
            r_le, r_te = R[k_le], R[k_te]
            cm_le, cu_le = CM[k_le], CU[k_le]
            b2m = _beta2_metal(row, r_te, omega)
            # Cu₂ lo impone el álabe: Cu₂ = U − Cm·tanβ₂. Cm sale del
            # equilibrio radial, que a su vez depende de h₀ (Euler con
            # Cu₂): se itera el par.
            # Arranque en caliente desde la pasada exterior anterior: el
            # lazo del cierre del álabe converge en 3-4 iteraciones en vez
            # de 25 cuando las líneas ya casi no se mueven.
            cm_te = CM[k_te].copy()
            ds_te = DS[k_te].copy() if it else DS[k_le].copy()
            ds_rot = np.zeros(n_sl)
            t0_te = np.full(n_sl, T0_in)
            for _ in range(40):
                cu_te = omega * r_te - cm_te * np.tan(b2m)
                h0_te = H0[k_le] + omega * (r_te * cu_te - r_le * cu_le)
                cm_new, _T, _lim = _solve_station_cm(
                    r_te, cu_te, h0_te, ds_te, TS[k_te], cm_te,
                    KAPPA[k_te], COSG[k_te], mdot, stations[k_te]["kb"],
                    T0_in, P0_in, float(np.mean(cm_te)))
                lim_hits += int(_lim)
                # --- perdidas del ROTOR, linea a linea, con ESE estado ---
                w1 = np.sqrt(cm_le ** 2 + (omega * r_le - cu_le) ** 2)
                w2 = np.sqrt(cm_new ** 2 + (omega * r_te - cu_te) ** 2)
                b1 = np.arctan2(omega * r_le - cu_le,
                                np.maximum(cm_le, 1e-3))
                _, p_le, rho_le = _static_state(
                    H0[k_le], cm_le ** 2 + cu_le ** 2, DS[k_le],
                    T0_in, P0_in)
                a_le = np.sqrt(pc.GAMMA * pc.RGAS
                               * np.maximum(TS[k_le], 1.0))
                a_te = np.sqrt(pc.GAMMA * pc.RGAS
                               * np.maximum(TS[k_te], 1.0))
                ds_rot = _row_loss_span(
                    b1, np.abs(b2m), w1, w2, w1 / a_le,
                    row["sigma_r"] * row["r_m"] / np.maximum(r_le, 1e-6),
                    row["h_blade"] / max(row["c_r"], 1e-4), rho_le,
                    row["c_r"], p_le, TS[k_le], h0_te,
                    r_le, w2 / a_te, _streamtube_ratios(r_le, r_te),
                    int(record["stage_table"][i]["n_blades_rotor"]),
                    "rotor")
                t0_te = _T_from_h_v(h0_te, np.full(n_sl, T0_in + 100.0))
                ds_rot = ds_rot + _tip_clearance_span(
                    r_te, h0_te - H0[k_le], row["h_blade"], t0_te)
                d_cm = float(np.max(np.abs(cm_new - cm_te)))
                cm_te = cm_te + BLADE_RELAX * (cm_new - cm_te)
                ds_te = ds_te + BLADE_RELAX * (DS[k_le] + ds_rot - ds_te)
                if d_cm < 1e-4 * max(float(np.mean(cm_new)), 1.0):
                    break
            cu_te = omega * r_te - cm_te * np.tan(b2m)
            h0_te = H0[k_le] + omega * (r_te * cu_te - r_le * cu_le)
            # mezcla radial efectiva en el hueco tras la fila (ver
            # SPAN_MIX_PER_ROW): sin ella la entropía de pared se acumula
            # fila a fila en las mismas líneas y en multietapa el
            # equilibrio radial se distorsiona hasta bloquear la marcha
            h0_te = _span_mix(h0_te)
            ds_te = _span_mix(ds_te)
            CU[k_te], H0[k_te], DS[k_te], CM[k_te] = (cu_te, h0_te,
                                                      ds_te, cm_te)
            _, TS[k_te], _lim = _solve_station_cm(
                r_te, cu_te, h0_te, ds_te, TS[k_te], cm_te, KAPPA[k_te],
                COSG[k_te], mdot, stations[k_te]["kb"], T0_in, P0_in,
                float(np.mean(cm_te)))

            # ===== ESTATOR: TE -> salida =================================
            if k_se >= n_stat:
                break
            r_se = R[k_se]
            a_out = _alpha_out_metal(row, r_se)
            cm_se = CM[k_se].copy()
            ds_se = DS[k_se].copy() if it else DS[k_te].copy()
            ds_sta = np.zeros(n_sl)
            for _ in range(40):
                cu_se = cm_se * np.tan(a_out)
                cm_new, _T, _lim = _solve_station_cm(
                    r_se, cu_se, H0[k_te], ds_se, TS[k_se], cm_se,
                    KAPPA[k_se], COSG[k_se], mdot, stations[k_se]["kb"],
                    T0_in, P0_in, float(np.mean(cm_se)))
                lim_hits += int(_lim)
                c2v = np.sqrt(CM[k_te] ** 2 + CU[k_te] ** 2)
                c3v = np.sqrt(cm_new ** 2 + cu_se ** 2)
                a2v = np.arctan2(CU[k_te], np.maximum(CM[k_te], 1e-3))
                _, p_te, rho_te = _static_state(
                    H0[k_te], c2v ** 2, DS[k_te], T0_in, P0_in)
                a_te = np.sqrt(pc.GAMMA * pc.RGAS
                               * np.maximum(TS[k_te], 1.0))
                a_se = np.sqrt(pc.GAMMA * pc.RGAS
                               * np.maximum(TS[k_se], 1.0))
                ds_sta = _row_loss_span(
                    np.abs(a2v), np.abs(a_out), c2v, c3v, c2v / a_te,
                    row["sigma_s"] * row["r_m"] / np.maximum(r_te, 1e-6),
                    row["h_blade"] / max(row["c_s"], 1e-4), rho_te,
                    row["c_s"], p_te, TS[k_te], H0[k_te],
                    r_te, c3v / a_se, _streamtube_ratios(r_te, r_se),
                    int(record["stage_table"][i]["n_blades_stator"]),
                    "stator")
                ds_new = DS[k_te] + ds_sta + _endwall_span(
                    r_se, record["stage_table"][i],
                    H0[k_te] - H0[k_le], t0_te)
                d_cm = float(np.max(np.abs(cm_new - cm_se)))
                cm_se = cm_se + BLADE_RELAX * (cm_new - cm_se)
                ds_se = ds_se + BLADE_RELAX * (ds_new - ds_se)
                if d_cm < 1e-4 * max(float(np.mean(cm_new)), 1.0):
                    break
            cu_se = cm_se * np.tan(a_out)
            # misma mezcla radial en el hueco tras el estátor
            h0_se = _span_mix(H0[k_te])
            ds_se = _span_mix(ds_se)
            CU[k_se], H0[k_se], DS[k_se], CM[k_se] = (cu_se, h0_se,
                                                      ds_se, cm_se)
            _, TS[k_se], _lim = _solve_station_cm(
                r_se, cu_se, H0[k_se], ds_se, TS[k_se], cm_se, KAPPA[k_se],
                COSG[k_se], mdot, stations[k_se]["kb"], T0_in, P0_in,
                float(np.mean(cm_se)))

            rows_out.append(dict(stage=i, beta2_metal_deg=np.degrees(b2m),
                                 alpha_out_deg=np.degrees(a_out),
                                 ds_rotor=ds_rot, ds_stator=ds_sta))

        # ---- recolocar líneas y actualizar curvatura ---------------------
        for k, stn in enumerate(stations):
            c2v = CM[k] ** 2 + CU[k] ** 2
            _, _, rho = _static_state(H0[k], c2v, DS[k], T0_in, P0_in)
            r_new = _reposition(stn["r_hub"], stn["r_tip"], R[k], rho,
                                CM[k], COSG[k], stn["kb"], frac)
            R[k] = R[k] + SCM_RELAX_R * (r_new - R[k])
        for j in range(n_sl):
            kap, cg = _curvature(z_stat, R[:, j])
            KAPPA[:, j] += SCM_RELAX_CURV * (kap - KAPPA[:, j])
            COSG[:, j] = cg

        move = float(np.max(np.abs(R - R_old))
                     / max(float(np.mean(R)), 1e-9))
        if move < SCM_TOL:
            converged = True
            break

    if not converged:
        raise SCMDiverged(f"sin converger en {SCM_MAX_ITER} iteraciones "
                          f"(movimiento residual {move:.2e})")
    # ¿El limitador de perfil AMASÓ el campo final? La versión anterior
    # rechazaba si el limitador se había tocado en CUALQUIER llamada de la
    # última pasada — pero el cierre del álabe pasa por transitorios que
    # el limitador estabiliza (para eso está) y que luego convergen a un
    # campo sano lejos del límite: el E³ entero se rechazaba por golpes en
    # iteraciones intermedias con un campo final impecable. El criterio
    # honesto es el CAMPO CONVERGIDO: si alguna estación tiene el perfil
    # clavado en los límites, el número que saldría es del limitador y no
    # de la física — se degrada a L0 etiquetado. Si el campo final vive
    # dentro de los límites, los transitorios son historia del solver.
    lim_final = 0
    for k in range(n_stat):
        mean = max(float(np.mean(CM[k])), 1e-6)
        if (float(np.min(CM[k])) <= CM_SPREAD_MIN * mean * 1.001
                or float(np.max(CM[k])) >= CM_SPREAD_MAX * mean * 0.999):
            lim_final += 1
    if lim_final:
        raise SCMDiverged(
            f"el limitador de perfil de Cm sigue activo al converger "
            f"({lim_final} estaciones con el perfil clavado en "
            f"[{CM_SPREAD_MIN}, {CM_SPREAD_MAX}]x la media)")

    # ---- IGV y OGV: las dos filas que el SCM no modela --------------------
    # Son filas de guía sin trabajo que ocupan todo el span; el meanline ya
    # calcula su pérdida (el IGV al entrar, el OGV al salir quitando el
    # remolino residual). Su entropía se recupera como el resto de la de
    # máquina que no está en las etapas, y se aplica UNIFORME en el span
    # sobre la última estación, junto con la condición de salida AXIAL que
    # deja el OGV. Sin este paso el PR del SCM se compara contra un P₀ que
    # todavía no ha pagado esas dos filas.
    ds_stages = sum(float(s.get("ds_stage_J_kgK", 0.0))
                    for s in record["stage_table"])
    ds_guide = max(float(record.get("ds_machine_J_kgK", ds_stages))
                   - ds_stages, 0.0)
    k = n_stat - 1
    DS[k] = DS[k] + ds_guide
    CU[k] = np.zeros(n_sl)
    CM[k], TS[k], _lim_exit = _solve_station_cm(
        R[k], CU[k], H0[k], DS[k], TS[k], CM[k], KAPPA[k], COSG[k],
        mdot, stations[k]["kb"], T0_in, P0_in, float(np.mean(CM[k])))

    # ---- promediado por gasto a la salida --------------------------------
    trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

    def _plane(q: int) -> tuple[float, float]:
        """(P₀, h₀) promediados por GASTO en la estación q."""
        rq, cmq, cuq = R[q], CM[q], CU[q]
        Tq, pq, rhoq = _static_state(H0[q], cmq ** 2 + cuq ** 2, DS[q],
                                     T0_in, P0_in)
        wq = (rhoq * cmq * COSG[q] * 2.0 * math.pi * rq
              * stations[q]["kb"])
        mq = float(trapz(wq, rq))
        if mq <= 0:
            raise SCMDiverged("gasto nulo en el promediado")
        return (float(trapz(wq * _total_pressure(Tq, pq, H0[q]), rq)) / mq,
                float(trapz(wq * H0[q], rq)) / mq)

    r, cm, cu = R[k], CM[k], CU[k]
    _T_ex, _p_ex, _ = _static_state(H0[k], cm ** 2 + cu ** 2, DS[k],
                                    T0_in, P0_in)
    P0_sl = _total_pressure(_T_ex, _p_ex, H0[k])
    P0_out, h0_out = _plane(k)
    T0_out = pc.T_from_h(h0_out, T0_in + 100.0)

    # Plano de SALIDA DEL ROTOR de la etapa 1 (estación 1). El Rotor 37 y
    # el 67 están medidos ahí —rotor aislado, sin estátor detrás— y sin
    # emitir este plano no se pueden calificar contra L1 aunque el solver
    # resuelva: lo único que salía era la salida de máquina, que es otro
    # plano. Emitirlo es lo que convierte esas dos máquinas en anclajes.
    rotor1: dict = {}
    if n_stat > 1:
        try:
            P0_r1, h0_r1 = _plane(1)
            T0_r1 = pc.T_from_h(h0_r1, T0_in + 100.0)
            pr_r1 = P0_r1 / P0_in
            dphi_r1 = pc.phi_air(T0_r1) - pc.phi_air(T0_in)
            if pr_r1 > 1.001 and dphi_r1 > 1e-6:
                T_id = pc.T_from_phi(
                    pc.phi_air(T0_in) + pc.RGAS * math.log(pr_r1), T0_r1)
                rotor1 = dict(
                    PR=float(pr_r1),
                    T0_out=float(T0_r1),
                    eta_poly=float(np.clip(
                        pc.RGAS * math.log(pr_r1) / dphi_r1, 0.05, 0.999)),
                    eta_isen=float(np.clip(
                        (pc.h_air(T_id) - h0_in)
                        / max(h0_r1 - h0_in, 1e-6), 0.05, 0.999)))
        except SCMDiverged:
            rotor1 = {}

    PR = P0_out / P0_in
    dphi = pc.phi_air(T0_out) - pc.phi_air(T0_in)
    if not (np.isfinite(PR) and PR > 1.001 and dphi > 1e-6):
        raise SCMDiverged(f"salida no física (PR={PR}, dphi={dphi})")
    pr_l0 = float(record.get("PR", 0.0) or 0.0)
    if pr_l0 > 1.0 and abs(PR / pr_l0 - 1.0) > PR_WINDOW:
        raise SCMDiverged(
            f"DERIVA respecto a L0: PR {pr_l0:.3f} → {PR:.3f} "
            f"({100 * (PR / pr_l0 - 1):+.1f}%, ventana "
            f"±{100 * PR_WINDOW:.0f}%). El annulus lo dimensionó L0 con su "
            f"Cx uniforme; sobre {len(record['stage_table'])} etapas la "
            "diferencia se compone a través del álabe de ángulo fijo")
    eta_poly = float(np.clip(pc.RGAS * math.log(PR) / dphi, 0.05, 0.999))

    # eficiencia isentrópica exacta de gas imperfecto
    T_id = pc.T_from_phi(pc.phi_air(T0_in) + pc.RGAS * math.log(PR),
                         T0_out)
    eta_isen = float(np.clip(
        (pc.h_air(T_id) - h0_in) / max(h0_out - h0_in, 1e-6), 0.05, 0.999))

    span = (r - r[0]) / max(r[-1] - r[0], 1e-9)
    return dict(
        PR=float(PR), eta_poly=eta_poly, eta_isen=eta_isen,
        T0_out=float(T0_out), source="scm_L1",
        scm=dict(
            rotor1=rotor1,
            n_streamlines=n_sl, n_stations=n_stat, iterations=it + 1,
            residual=move,
            span_frac=span.tolist(),
            r_exit_mm=(r * 1000.0).tolist(),
            cm_exit=cm.tolist(), cu_exit=cu.tolist(),
            P0_exit=P0_sl.tolist(),
            dh0_span=(H0[k] - h0_in).tolist(),
            ds_span=DS[k].tolist(),
            work_spread=float((np.max(H0[k]) - np.min(H0[k]))
                              / max(h0_out - h0_in, 1e-6)),
            cm_mean_stations=[float(np.mean(CM[q])) for q in range(n_stat)],
            cu_mean_stations=[float(np.mean(CU[q])) for q in range(n_stat)],
            h0_mean_stations=[float(np.mean(H0[q])) for q in range(n_stat)],
            ds_mean_stations=[float(np.mean(DS[q])) for q in range(n_stat)],
            r_hub_mm=[s["r_hub"] * 1000.0 for s in stations],
            r_tip_mm=[s["r_tip"] * 1000.0 for s in stations],
            z_mm=(z_stat * 1000.0).tolist(),
            streamlines_mm=(R * 1000.0).tolist(),
        ))


# ---------------------------------------------------------------------------
# Pérdidas resueltas en el span
# ---------------------------------------------------------------------------
def _blade_thickness_law(kind: str) -> tuple[float, float, float, float]:
    """(t/c en cubo, t/c en punta, taper de cuerda, semi-espesor del BA).

    Se leen de la capa 5a: es la ley de espesor del álabe que REALMENTE se
    fabrica, no una constante paralela. Que L1 pague la pérdida de choque
    con el espesor de la sección que sale en el STEP es la mitad de la
    razón de tener un nivel L1.
    """
    import geometry_generator as gg          # diferido: 5a importa 1
    tc = ((gg.TC_ROOT_R, gg.TC_TIP_R) if kind == "rotor"
          else (gg.TC_ROOT_S, gg.TC_TIP_S))
    # El generador DCA de blade_profiles acota el semi-espesor a 0.02·t/c
    # en el borde: ese es el t_LE del perfil biconvexo que se voxeliza.
    return tc[0], tc[1], gg.CHORD_TAPER, 0.02


def _streamtube_ratios(r1: np.ndarray,
                       r2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(A₂/A₁, h₁/h₂) del TUBO DE CORRIENTE de cada línea en la fila.

    Koch & Smith corrigen tanto el pico de succión (ec. 21-22) como el
    espesor de momento (Fig. 4a) por la contracción del tubo, con una
    relación de alturas que en un meanline hay que estimar del annulus.
    Aquí es exacta: la altura del tubo es la separación entre líneas
    vecinas, y el área va con r·h porque el tubo es un anillo.
    """
    def h(r: np.ndarray) -> np.ndarray:
        d = np.empty(len(r))
        d[1:-1] = 0.5 * (r[2:] - r[:-2])
        d[0] = r[1] - r[0]
        d[-1] = r[-1] - r[-2]
        return np.maximum(np.abs(d), 1e-6)
    h1, h2 = h(r1), h(r2)
    a1 = np.maximum(r1, 1e-6) * h1
    a2 = np.maximum(r2, 1e-6) * h2
    return (np.clip(a2 / np.maximum(a1, 1e-12), 0.4, 2.0),
            np.clip(h1 / np.maximum(h2, 1e-12), 0.6, 2.0))


def _wall_band_weight(r: np.ndarray, frac: float) -> np.ndarray:
    """Peso que concentra una pérdida en las bandas de CUBO y PUNTA.

    Normalizado a media 1: redistribuye sin cambiar el total, así que
    mueve el PERFIL radial de la pérdida sin tocar su nivel.
    """
    span = (r - r[0]) / max(r[-1] - r[0], 1e-9)
    w = np.maximum(np.clip((frac - span) / frac, 0.0, 1.0),
                   np.clip((span - (1.0 - frac)) / frac, 0.0, 1.0))
    s = float(np.sum(w))
    if s <= 0:
        return np.ones(len(r))
    return w * len(r) / s


def _row_loss_span(b1: np.ndarray, b2: np.ndarray, w1: np.ndarray,
                   w2: np.ndarray, m1: np.ndarray, sigma: np.ndarray,
                   h_over_c: float, rho1: np.ndarray, chord: float,
                   p1: np.ndarray, T1: np.ndarray,
                   h0_out: np.ndarray, r1: np.ndarray, m2: np.ndarray,
                   ratios: tuple[np.ndarray, np.ndarray], n_blades: int,
                   kind: str) -> np.ndarray:
    """Δs de una fila, línea de corriente a línea de corriente.

    Es la ganancia real de resolver el span: el meanline calcula UNA
    pérdida con el triángulo medio y la reparte por igual. Aquí cada línea
    paga su propia difusión, su propio Reynolds y su propio Mach.

    El CHOQUE va por el modelo de Koch & Smith 1976 (physics_core §1b), y
    esto es lo que solo se puede hacer con el span resuelto: el Mach
    representativo del choque de pasaje sale del pico de succión de ESA
    sección, que depende de su solidez, su espesor, su ángulo y la
    contracción de SU tubo de corriente — cuatro cosas que en el meanline
    solo existen en la línea media. Y el romo del borde de ataque paga
    contra el espaciado tangencial local, que va con el radio.

    La SECUNDARIA de Howell se redistribuye a las bandas de pared. Koch &
    Smith suman perfil y choque sección a sección y tratan el end-wall
    aparte, a nivel de etapa, porque el vórtice de pasaje es un fenómeno
    de pared: untarlo plano en el span —como hacía esto— contradice la
    razón de existir del módulo. La media se conserva, así que cambia el
    PERFIL de la pérdida, no su nivel; y ese perfil entra en el equilibrio
    radial por el término T·∂s/∂r.
    """
    n = len(b1)
    area_ratio, h_ratio = ratios
    tc_hub, tc_tip, taper, le_frac = _blade_thickness_law(kind)
    span = (r1 - r1[0]) / max(r1[-1] - r1[0], 1e-9)
    om_pp = np.empty(n)     # perfil + choque, se quedan donde ocurren
    om_sec = np.empty(n)    # secundaria, se redistribuye a las paredes
    ds_le = np.empty(n)
    for j in range(n):
        re_c = float(rho1[j] * w1[j] * chord / pc.MU_AIR)
        t_c = tc_hub + (tc_tip - tc_hub) * float(span[j])
        c_loc = chord * (1.0 + (taper - 1.0) * float(span[j]))
        b_tan = 2.0 * math.pi * float(r1[j]) / max(n_blades, 1)
        _, _, det = pc._row_losses(
            float(abs(b1[j])), float(abs(b2[j])), float(w1[j]),
            float(w2[j]), float(max(sigma[j], 0.2)), h_over_c,
            float(m1[j]), float(m1[j]), re_c=re_c,
            gam=pc.gamma_air(float(T1[j])),
            ks=dict(M2=float(m2[j]), tmax_c=t_c,
                    area_ratio=float(area_ratio[j]), T1=float(T1[j]),
                    h1_over_h2=float(h_ratio[j]),
                    t_le_over_b=2.0 * le_frac * t_c * c_loc
                    / max(b_tan, 1e-6)))
        om_pp[j] = det["profile"] + det["shock"]
        om_sec[j] = det["endwall"]
        ds_le[j] = det["ds_le"]
    om = om_pp + float(np.mean(om_sec)) * _wall_band_weight(
        r1, WALL_BAND_FRAC)

    ds = np.empty(n)
    for j in range(n):
        # ω̄ referida a la cabeza dinámica COMPRESIBLE (P₀−p), igual que el
        # meanline desde la fase 9
        kx = pc.gamma_air(float(T1[j])) / (pc.gamma_air(float(T1[j])) - 1.0)
        T01 = float(T1[j]) + float(w1[j]) ** 2 / (
            2.0 * pc.cp_air(float(T1[j])))
        P01 = float(p1[j]) * (T01 / max(float(T1[j]), 1.0)) ** kx
        q1 = max(P01 - float(p1[j]), 1e-3)
        P02 = max(P01 - om[j] * q1, 0.05 * P01)
        ds[j] = pc.RGAS * math.log(P01 / P02) + ds_le[j]
    return ds


def _tip_clearance_span(r: np.ndarray, dh0: np.ndarray, h_blade: float,
                        t0: np.ndarray) -> np.ndarray:
    """Pérdida de holgura, concentrada en la banda de PUNTA.

    El meanline la reparte por igual sobre toda la etapa; el vórtice de
    holgura vive en el 20-30% exterior del span. Repartirla donde ocurre
    cambia el perfil de rendimiento, que es la información que el
    meanline no tiene.
    """
    eps_h = min(pc.TIP_CLEARANCE_MM / max(h_blade * 1000.0, 1e-3),
                pc.EPS_H_MAX)
    frac = pc._clearance_loss_frac(eps_h)
    span = (r - r[0]) / max(r[-1] - r[0], 1e-9)
    w = np.clip((span - (1.0 - TIP_BAND_FRAC)) / TIP_BAND_FRAC, 0.0, 1.0)
    if float(np.sum(w)) <= 0:
        return np.zeros(len(r))
    # se reparte de forma que la MEDIA en el span sea la del meanline
    w = w * len(r) / float(np.sum(w))
    dh = frac * np.maximum(dh0, 0.0) * w
    return dh / np.maximum(t0, 1.0)


def _endwall_span(r: np.ndarray, stage: dict, dh0: np.ndarray,
                  t0: np.ndarray) -> np.ndarray:
    """Débito de capa límite de pared (Koch & Smith), repartido sobre las
    bandas de cubo y punta en vez de uniformemente."""
    a_ew = float(np.clip(stage.get("delta_star_sum_mm", 0.0)
                         / max(stage.get("h_blade_mm", 1.0), 1e-3),
                         0.0, 0.35))
    if a_ew <= 0.0:
        return np.zeros(len(r))
    f_ew = (1.0 - a_ew) / max(1.0 - pc.NU_OVER_DELTA * a_ew, 1e-3)
    w = _wall_band_weight(r, WALL_BAND_FRAC)
    dh = np.maximum(dh0, 0.0) * max(1.0 - f_ew, 0.0) * w
    return dh / np.maximum(t0, 1.0)
