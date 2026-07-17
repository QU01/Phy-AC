"""
QUASAR Phy-AC · blade_profiles.py
=================================
Perfiles de compresor axial y correlaciones de ángulo metálico (capa 5a,
inspirado en turbodesigner/MetalAngles de OpenOrion pero NumPy puro).

Contenido:
  * naca65_profile / dca_profile — secciones 2D de cuerda unitaria sobre
    línea de comba de ARCO CIRCULAR: NACA 65 usa la distribución de espesor
    65-010 tabulada (escalada); DCA usa espesor biconvexo parabólico
    (aprox. estándar a nivel meanline). Contrato: polígono CERRADO de
    exactamente N puntos, CCW, empezando en el borde de ataque, sin
    autointersección — la capa C# hace loft entre secciones y EXIGE conteo
    igual y winding consistente.
  * incidence_design — incidencia de diseño i* (fits de Aungier 2003 de las
    cartas de Lieblein, NASA SP-36 cap. VI).
  * carter_deviation — regla de Carter δ = m·θ/√σ, m = 0.23(2a/c)² + χ₂/500
    (arco circular, a/c=0.5).
  * metal_angles — punto fijo χ₁ = β₁−i*, χ₂ = β₂−δ(θ), θ = χ₁−χ₂.

Convención de ángulos: grados, medidos desde el eje de la máquina, en el
marco de la fila (β para rotor, α para estátor).
"""

from __future__ import annotations

import math

import numpy as np

# Distribución de espesor NACA 65-010 (x/c, t/2 en % de cuerda) — tabla
# clásica (Abbott & von Doenhoff / NASA SP-36); se escala por t/c.
_NACA65_XT = np.array([
    0.000, 0.005, 0.0075, 0.0125, 0.025, 0.050, 0.075, 0.100, 0.150,
    0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600,
    0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 1.000])
_NACA65_HT = np.array([
    0.000, 0.772, 0.932, 1.169, 1.574, 2.177, 2.647, 3.040, 3.666,
    4.143, 4.503, 4.760, 4.924, 4.996, 4.963, 4.812, 4.530, 4.146,
    3.682, 3.156, 2.584, 1.987, 1.385, 0.810, 0.306, 0.000]) / 100.0

N_POINTS = 60   # puntos por sección (contrato con la capa C#)


def _circular_arc_camber(camber_deg: float, x: np.ndarray):
    """Línea de comba de arco circular en cuerda unitaria.

    Devuelve (yc, dyc/dx). Ángulo de comba θ = χ1−χ2; el arco subtiende θ
    con la cuerda: yc(x) máx en x=0.5, pendiente = tan(θ/2) en el LE.
    """
    th = math.radians(max(abs(camber_deg), 1e-4))
    R = 1.0 / (2.0 * math.sin(th / 2.0))
    y0 = R * math.cos(th / 2.0)
    yc = np.sqrt(np.maximum(R * R - (x - 0.5) ** 2, 0.0)) - y0
    dyc = -(x - 0.5) / np.sqrt(np.maximum(R * R - (x - 0.5) ** 2, 1e-9))
    if camber_deg < 0:
        yc, dyc = -yc, -dyc
    return yc, dyc


def _thickness(x: np.ndarray, t_over_c: float, profile: str) -> np.ndarray:
    """Semi-espesor t(x)/2 en cuerda unitaria."""
    if profile.upper() == "NACA65":
        ht = np.interp(x, _NACA65_XT, _NACA65_HT) * (t_over_c / 0.10)
    else:  # DCA — biconvexo (parabólico, aprox. del arco circular)
        ht = 0.5 * t_over_c * (1.0 - (2.0 * x - 1.0) ** 2)
        # bordes redondeados mínimos para robustez del voxelizado
        ht = np.maximum(ht, 0.02 * t_over_c)
        ht[x <= 0.0] = 0.0
        ht[x >= 1.0] = 0.0
    return ht


def _profile(camber_deg: float, t_over_c: float, profile: str,
             n: int = N_POINTS) -> np.ndarray:
    """Polígono cerrado (n, 2) CCW empezando en el LE, cuerda unitaria."""
    m = n // 2 + 1
    # muestreo coseno: densifica LE/TE
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, math.pi, m)))
    yc, dyc = _circular_arc_camber(camber_deg, x)
    ht = _thickness(x, t_over_c, profile)
    th = np.arctan(dyc)
    xu = x - ht * np.sin(th)
    yu = yc + ht * np.cos(th)
    xl = x + ht * np.sin(th)
    yl = yc - ht * np.cos(th)
    # CCW desde el LE: lado de presión (inferior) LE→TE, succión TE→LE
    pts = np.concatenate([
        np.column_stack([xl, yl]),          # LE → TE por abajo
        np.column_stack([xu, yu])[::-1][1:],  # TE → LE por arriba
    ])
    # exactamente n puntos (sin repetir el LE al final)
    if len(pts) > n:
        pts = pts[:n]
    while len(pts) < n:
        pts = np.vstack([pts, pts[-1] + (pts[0] - pts[-1]) * 0.5])
    return pts


def naca65_profile(camber_deg: float, t_over_c: float,
                   n: int = N_POINTS) -> np.ndarray:
    """Sección NACA 65 (comba de arco circular + espesor 65-010)."""
    return _profile(camber_deg, t_over_c, "NACA65", n)


def dca_profile(camber_deg: float, t_over_c: float,
                n: int = N_POINTS) -> np.ndarray:
    """Sección de doble arco circular (biconvexa)."""
    return _profile(camber_deg, t_over_c, "DCA", n)


# ---------------------------------------------------------------------------
# Correlaciones de incidencia y desviación
# ---------------------------------------------------------------------------
def incidence_design(beta1_deg: float, sigma: float, t_over_c: float,
                     camber_deg: float, profile: str = "NACA65") -> float:
    """Incidencia de diseño i* [deg] — fits de Aungier (2003) de las cartas
    de Lieblein (NASA SP-36 cap. VI): i* = Ksh·Kti·i0_10 + n·θ."""
    b1 = min(max(abs(beta1_deg), 0.0), 80.0)
    sig = max(sigma, 0.4)
    p = 0.914 + sig ** 3 / 160.0
    i0_10 = b1 ** p / (5.0 + 46.0 * math.exp(-2.3 * sig)) \
        - 0.1 * sig ** 3 * math.exp((b1 - 70.0) / 4.0)
    n_slope = 0.025 * sig - 0.06 - (b1 / 90.0) ** (1.0 + 1.2 * sig) / \
        (1.5 + 0.43 * sig)
    q = 0.28 / (0.1 + (t_over_c) ** 0.3)
    k_ti = (10.0 * t_over_c) ** q
    k_sh = 1.0 if profile.upper() == "NACA65" else 0.7   # DCA
    return k_sh * k_ti * i0_10 + n_slope * abs(camber_deg)


def carter_m(stagger_deg: float, a_over_c: float = 0.5) -> float:
    """Coeficiente m de Carter (comba de arco circular)."""
    return 0.23 * (2.0 * a_over_c) ** 2 + abs(stagger_deg) / 500.0


def carter_deviation(camber_deg: float, stagger_deg: float,
                     sigma: float) -> float:
    """Desviación de Carter δ = m·θ/√σ [deg]."""
    return carter_m(stagger_deg) * abs(camber_deg) / \
        math.sqrt(max(sigma, 0.2))


def metal_angles(beta1_flow: float, beta2_flow: float, sigma: float,
                 t_over_c: float, profile: str = "NACA65",
                 iters: int = 6) -> dict:
    """Ángulos metálicos por punto fijo (grados, desde el eje).

    χ1 = β1 − i*(θ),  χ2 = β2 − sgn(θ)·δ(θ, γ),  θ = χ1 − χ2,
    γ = (χ1+χ2)/2. Converge en 3-4 iteraciones para filas convencionales.

    El signo de la desviación sigue al del camber: el flujo siempre gira
    MENOS que el metal, así que en filas difusoras (θ>0, rotor/estátor)
    χ2 = β2 − δ y en filas aceleradoras (θ<0, IGV) χ2 = β2 + δ — el
    metal debe girar MÁS allá del ángulo de flujo objetivo en ambos casos.
    """
    chi1, chi2 = beta1_flow, beta2_flow
    i_star = dev = 0.0
    for _ in range(iters):
        camber = chi1 - chi2
        stagger = 0.5 * (chi1 + chi2)
        i_star = incidence_design(beta1_flow, sigma, t_over_c, camber,
                                  profile)
        dev = carter_deviation(camber, stagger, sigma)
        chi1_new = beta1_flow - i_star
        chi2_new = beta2_flow - math.copysign(dev, camber)
        if abs(chi1_new - chi1) < 1e-4 and abs(chi2_new - chi2) < 1e-4:
            chi1, chi2 = chi1_new, chi2_new
            break
        chi1 = 0.5 * (chi1 + chi1_new)
        chi2 = 0.5 * (chi2 + chi2_new)
    return dict(chi1_deg=chi1, chi2_deg=chi2, camber_deg=chi1 - chi2,
                stagger_deg=0.5 * (chi1 + chi2), incidence_deg=i_star,
                deviation_deg=dev)


if __name__ == "__main__":
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ma = metal_angles(55.0, 35.0, 1.2, 0.08)
    print("metal_angles(55→35, σ=1.2):",
          {k: round(v, 2) for k, v in ma.items()})
    for name, fn in (("NACA65", naca65_profile), ("DCA", dca_profile)):
        p = fn(ma["camber_deg"], 0.08)
        closed = np.linalg.norm(p[0] - p[-1]) < 0.2
        print(f"{name}: {len(p)} pts, x∈[{p[:,0].min():.3f},"
              f"{p[:,0].max():.3f}], t_max={2*p[:,1].max():.3f}")
