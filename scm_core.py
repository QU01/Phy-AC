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
TIP_BAND_FRAC = 0.25           # fracción de span sobre la que se reparte
#                                la pérdida de holgura de punta
WALL_BAND_FRAC = 0.30          # ídem para la capa límite de pared


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
        mach = float(np.max(np.sqrt(c2) / a))
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
    CM = np.zeros((n_stat, n_sl))
    for k, stn in enumerate(stations):
        CM[k] = mdot / max(stn["area"] * stn["kb"] * 1.2, 1e-6)
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
                ds_rot = _row_loss_span(
                    b1, np.abs(b2m), w1, w2, w1 / a_le,
                    row["sigma_r"] * row["r_m"] / np.maximum(r_le, 1e-6),
                    row["h_blade"] / max(row["c_r"], 1e-4), rho_le,
                    row["c_r"], p_le, TS[k_le], h0_te)
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
                ds_sta = _row_loss_span(
                    np.abs(a2v), np.abs(a_out), c2v, c3v, c2v / a_te,
                    row["sigma_s"] * row["r_m"] / np.maximum(r_te, 1e-6),
                    row["h_blade"] / max(row["c_s"], 1e-4), rho_te,
                    row["c_s"], p_te, TS[k_te], H0[k_te])
                ds_new = DS[k_te] + ds_sta + _endwall_span(
                    r_se, record["stage_table"][i],
                    H0[k_te] - H0[k_le], t0_te)
                d_cm = float(np.max(np.abs(cm_new - cm_se)))
                cm_se = cm_se + BLADE_RELAX * (cm_new - cm_se)
                ds_se = ds_se + BLADE_RELAX * (ds_new - ds_se)
                if d_cm < 1e-4 * max(float(np.mean(cm_new)), 1.0):
                    break
            cu_se = cm_se * np.tan(a_out)
            CU[k_se], H0[k_se], DS[k_se], CM[k_se] = (cu_se, H0[k_te],
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
    if lim_hits:
        # El limitador de perfil sigue activo en la última pasada: el campo
        # que saldría está amasado por él, no por la física. Se degrada a
        # L0 ETIQUETADO en vez de devolver un número con buena pinta.
        raise SCMDiverged(
            f"el limitador de perfil de Cm sigue activo al converger "
            f"({lim_hits} estaciones): el equilibrio radial pide un perfil "
            f"fuera de [{CM_SPREAD_MIN}, {CM_SPREAD_MAX}]x la media")

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
    r, cm, cu = R[k], CM[k], CU[k]
    c2v = cm ** 2 + cu ** 2
    T, p, rho = _static_state(H0[k], c2v, DS[k], T0_in, P0_in)
    P0_sl = _total_pressure(T, p, H0[k])
    w = rho * cm * COSG[k] * 2.0 * math.pi * r * stations[k]["kb"]
    trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    m_tot = float(trapz(w, r))
    if m_tot <= 0:
        raise SCMDiverged("gasto de salida nulo")
    P0_out = float(trapz(w * P0_sl, r)) / m_tot
    h0_out = float(trapz(w * H0[k], r)) / m_tot
    T0_out = pc.T_from_h(h0_out, T0_in + 100.0)

    PR = P0_out / P0_in
    dphi = pc.phi_air(T0_out) - pc.phi_air(T0_in)
    if not (np.isfinite(PR) and PR > 1.001 and dphi > 1e-6):
        raise SCMDiverged(f"salida no física (PR={PR}, dphi={dphi})")
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
def _row_loss_span(b1: np.ndarray, b2: np.ndarray, w1: np.ndarray,
                   w2: np.ndarray, m1: np.ndarray, sigma: np.ndarray,
                   h_over_c: float, rho1: np.ndarray, chord: float,
                   p1: np.ndarray, T1: np.ndarray,
                   h0_out: np.ndarray) -> np.ndarray:
    """Δs de una fila, línea de corriente a línea de corriente.

    Es la ganancia real de resolver el span: el meanline calcula UNA
    pérdida con el triángulo medio y la reparte por igual. Aquí cada línea
    paga su propia difusión, su propio Reynolds y su propio Mach — y el
    choque solo lo pagan las líneas donde el Mach relativo lo justifica,
    que es la mitad exterior de un rotor transónico.
    """
    n = len(b1)
    ds = np.empty(n)
    for j in range(n):
        re_c = float(rho1[j] * w1[j] * chord / pc.MU_AIR)
        om, _, _ = pc._row_losses(
            float(abs(b1[j])), float(abs(b2[j])), float(w1[j]),
            float(w2[j]), float(max(sigma[j], 0.2)), h_over_c,
            float(m1[j]), float(m1[j]), re_c=re_c,
            gam=pc.gamma_air(float(T1[j])))
        # ω̄ referida a la cabeza dinámica COMPRESIBLE (P₀−p), igual que el
        # meanline desde la fase 9
        kx = pc.gamma_air(float(T1[j])) / (pc.gamma_air(float(T1[j])) - 1.0)
        T01 = float(T1[j]) + float(w1[j]) ** 2 / (
            2.0 * pc.cp_air(float(T1[j])))
        P01 = float(p1[j]) * (T01 / max(float(T1[j]), 1.0)) ** kx
        q1 = max(P01 - float(p1[j]), 1e-3)
        P02 = max(P01 - om * q1, 0.05 * P01)
        ds[j] = pc.RGAS * math.log(P01 / P02)
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
    span = (r - r[0]) / max(r[-1] - r[0], 1e-9)
    w = np.maximum(
        np.clip((WALL_BAND_FRAC - span) / WALL_BAND_FRAC, 0.0, 1.0),
        np.clip((span - (1.0 - WALL_BAND_FRAC)) / WALL_BAND_FRAC, 0.0, 1.0))
    if float(np.sum(w)) <= 0:
        return np.zeros(len(r))
    w = w * len(r) / float(np.sum(w))
    dh = np.maximum(dh0, 0.0) * max(1.0 - f_ew, 0.0) * w
    return dh / np.maximum(t0, 1.0)
