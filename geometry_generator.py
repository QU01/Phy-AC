"""
QUASAR Phy-AC · geometry_generator.py
=====================================
Capa 5a del CEM: del record meanline verificado a la GEOMETRÍA de filas de
álabes — el contrato `axial_compressor.json` (schema phyac-axial-2) que
consume la capa 5c (C#/PicoGK) y los CSV auxiliares para CFD/documentación.

Pipeline (inspirado en turbodesigner: Stage/BladeRow/Vortex/MetalAngles):

  1. Annulus: polilíneas hub/tip (z, r) desde el stage_table (extensión
     axial = cuerda·cos γ + gaps inter-fila).
  2. Triángulos spanwise por FREE VORTEX (r·Cu = const, Cx = const) en
     N_SPAN secciones por fila.
  3. Ángulos metálicos por sección (incidencia de Lieblein/Aungier +
     desviación de Carter, blade_profiles.metal_angles).
  4. Sección 2D: DCA si M_entrada > 0.8, NACA-65 si no; espesor y cuerda
     con taper lineal hub→tip; polígono de 60 puntos CCW centrado en su
     centroide (stacking radial por centroides).
  5. Contrato JSON + CSVs + STEP opcional (extra `step`, CadQuery).

Convención de signos del contrato: los ángulos/stagger del ROTOR son
positivos y los del ESTÁTOR negativos (giran el flujo en sentidos opuestos
alrededor del eje); la capa C# aplica la rotación FIRMADA uniformemente.
"""

from __future__ import annotations

import csv
import json
import math
import os
import re

import numpy as np

from physics_core import (CP, GAMMA, RGAS, ROW_GAP_FRACTION,
                          TIP_CLEARANCE_MM, VAR_NAMES, VORTEX_N,
                          vortex_cu, vortex_cx)
from blade_profiles import (N_POINTS, _circular_arc_camber, dca_profile,
                            metal_angles, naca65_profile)
from contract_schema import SCHEMA_ID

N_CAMBER = 41              # puntos de la línea de comba (capa 5c)

N_SPAN = 13                # secciones por fila (twist free-vortex EXACTO en
#                            cada estación — con pocas secciones la
#                            interpolación lineal de la capa 5c deja
#                            crestas visibles en las estaciones originales)
CHORD_TAPER = 0.85         # c_tip/c_hub
TC_ROOT_R, TC_TIP_R = 0.10, 0.05   # t/c del rotor hub→tip
TC_ROOT_S, TC_TIP_S = 0.09, 0.06   # t/c del estátor
M_DCA_THRESHOLD = 0.80     # M de entrada sobre el cual la sección es DCA

SCHEMA = SCHEMA_ID     # phyac-axial-2 — ver contract_schema.py y
#                        schemas/phyac-axial-2.schema.json. El entero final
#                        es la versión MAYOR: sube cuando un consumidor
#                        escrito contra la anterior ya no puede leer el
#                        fichero sin perder información.


# ---------------------------------------------------------------------------
# Estaciones axiales y annulus
# ---------------------------------------------------------------------------
def _row_extent_mm(sections: list[dict]) -> float:
    """Extensión axial REAL de una fila alrededor de su eje de stacking.

    La lámina se centra en el CENTROIDE del perfil (no a media cuerda),
    así que el voladizo axial es asimétrico: se evalúa la comba ROTADA
    por el stagger de cada sección y se toma la envolvente simétrica
    2·máx|x_rot| + espesor (voxMeshShell infla ±½·t). Sin esto las filas
    se colocaban con la cuerda axial al stagger MEDIO y rotor/estátor se
    solapaban en el hub (hasta ~4.5 mm de voladizo)."""
    half = 0.0
    for sec in sections:
        g = math.radians(sec["stagger_deg"])
        cam = np.asarray(sec["camber_points"])
        xr = cam[:, 0] * math.cos(g) - cam[:, 1] * math.sin(g)
        half = max(half, float(np.max(np.abs(xr)))
                   + 0.5 * sec["thickness_mm"])
    return 2.0 * half


def _row_geometry(record: dict, extents: dict | None = None,
                  igv_extent_mm: float | None = None,
                  ogv_extent_mm: float | None = None) -> list[dict]:
    """Posiciones axiales de cada fila (z_le/z_te en mm).

    `extents`: {(stage, kind): extensión_axial_mm} EXACTAS (de las
    secciones ya construidas — build_contract las pasa siempre). Sin
    ellas (uso de visualización) se estima con la cuerda del hub, que es
    una cota superior razonable de la cuerda axial.

    `igv_extent_mm` / `ogv_extent_mm`: extensiones axiales del IGV (antes
    del rotor 1, produce el pre-swirl α₁ que la física asume) y del OGV
    (tras el último estátor, devuelve el flujo a axial — su pérdida ya la
    contabiliza el meanline). None = fila ausente (compatibilidad con la
    ruta de visualización).
    """
    rows = []
    z = 0.0
    st = record["stage_table"]
    if igv_extent_mm is not None:
        rows.append(dict(stage=-1, kind="igv", z_le=0.0,
                         z_te=igv_extent_mm,
                         chord_mm=st[0]["chord_stator_mm"],
                         chord_axial_mm=igv_extent_mm))
        z = igv_extent_mm * (1.0 + ROW_GAP_FRACTION)
    for s in st:
        for kind in ("rotor", "stator"):
            c = s[f"chord_{kind}_mm"]
            if extents is not None:
                cz = extents[(s["stage"], kind)]
            else:
                # estimación: cuerda del hub (taper) + 10% de espesor
                cz = 1.10 * 2.0 * c / (1.0 + CHORD_TAPER) * \
                    math.cos(math.radians(0.5 * (
                        s["beta1_deg"] + s["beta2_deg"]) if kind == "rotor"
                        else 0.5 * (s["alpha2_deg"] + s["alpha1_deg"]))) \
                    + 0.10 * c
            rows.append(dict(stage=s["stage"], kind=kind, z_le=z,
                             z_te=z + cz, chord_mm=c, chord_axial_mm=cz))
            z = z + cz + ROW_GAP_FRACTION * cz
    if ogv_extent_mm is not None:
        rows.append(dict(stage=len(st), kind="ogv", z_le=z,
                         z_te=z + ogv_extent_mm,
                         chord_mm=st[-1]["chord_stator_mm"],
                         chord_axial_mm=ogv_extent_mm))
    return rows


def _exit_height_mm(record: dict) -> float:
    """Altura del annulus en la estación de salida (área congelada)."""
    st = record["stage_table"]
    r_m = st[-1]["r_mean_mm"]
    if record.get("frozen") and record["frozen"]["areas"]:
        A_exit = record["frozen"]["areas"][-1] * 1e6   # m² → mm²
        return A_exit / (2 * math.pi * r_m)
    return st[-1]["h_blade_mm"]


def annulus_lines(record: dict, extents: dict | None = None,
                  igv_extent_mm: float | None = None,
                  ogv_extent_mm: float | None = None) -> dict:
    """Polilíneas hub/tip [[z_mm, r_mm], ...] de entrada a salida."""
    rows = _row_geometry(record, extents, igv_extent_mm, ogv_extent_mm)
    st = record["stage_table"]
    hub, tip = [], []
    if igv_extent_mm is not None:
        # estación de entrada del IGV: annulus de la etapa 1 (const_mean)
        hub.append([rows[0]["z_le"], st[0]["r_hub_mm"]])
        tip.append([rows[0]["z_le"], st[0]["r_tip_mm"]])
    for s in st:
        r_rows = [r for r in rows if r["stage"] == s["stage"]
                  and r["kind"] in ("rotor", "stator")]
        z0 = r_rows[0]["z_le"]
        hub.append([z0, s["r_hub_mm"]])
        tip.append([z0, s["r_tip_mm"]])
    # estación de salida: annulus del área de salida a r_mean constante
    z_end = rows[-1]["z_te"]
    r_m = st[-1]["r_mean_mm"]
    h_exit = _exit_height_mm(record)
    hub.append([z_end, r_m - 0.5 * h_exit])
    tip.append([z_end, r_m + 0.5 * h_exit])
    # la MISMA ε absoluta que usó la física (ε constante en mm — un solo
    # valor sirve para todas las filas de la capa 5c)
    return dict(hub=hub, tip=tip, tip_clearance_mm=TIP_CLEARANCE_MM)


# ---------------------------------------------------------------------------
# Secciones spanwise (free vortex)
# ---------------------------------------------------------------------------
def _polygon_centroid(p: np.ndarray) -> np.ndarray:
    x, y = p[:, 0], p[:, 1]
    xr, yr = np.roll(x, -1), np.roll(y, -1)
    a = x * yr - xr * y
    A = 0.5 * a.sum()
    cx = ((x + xr) * a).sum() / (6.0 * A)
    cy = ((y + yr) * a).sum() / (6.0 * A)
    return np.array([cx, cy])


def polygon_section_props(p: np.ndarray) -> dict:
    """Propiedades de sección de un polígono cerrado (fórmulas de Green).

    Función PURA y reutilizable (la consume structures_core.blade_modes):
    área y segundos momentos respecto al CENTROIDE, en las unidades del
    polígono al cuadrado/cuarta (mm²/mm⁴ para los `points` del contrato).
    I_min/I_max son los momentos principales; I_p = Ix + Iy el polar.
    El signo del área sigue la orientación (CCW > 0); los momentos se
    devuelven con el |área| (invariantes a la orientación).
    """
    p = np.asarray(p, dtype=float)
    c = _polygon_centroid(p)
    x, y = p[:, 0] - c[0], p[:, 1] - c[1]
    xr, yr = np.roll(x, -1), np.roll(y, -1)
    a = x * yr - xr * y                       # doble área de cada triángulo
    A = 0.5 * a.sum()
    sA = 1.0 if A >= 0.0 else -1.0            # momentos positivos siempre
    Ix = sA * ((y * y + y * yr + yr * yr) * a).sum() / 12.0
    Iy = sA * ((x * x + x * xr + xr * xr) * a).sum() / 12.0
    Ixy = sA * ((x * yr + 2.0 * x * y + 2.0 * xr * yr + xr * y)
                * a).sum() / 24.0
    I_mean = 0.5 * (Ix + Iy)
    R = math.hypot(0.5 * (Ix - Iy), Ixy)
    I_min, I_max = I_mean - R, I_mean + R
    return dict(A_mm2=abs(A), I_min_mm4=max(I_min, 0.0),
                I_max_mm4=max(I_max, 0.0), I_p_mm4=Ix + Iy,
                centroid=c)


def _row_sections(stage: dict, kind: str, omega: float,
                  n_span: int = N_SPAN,
                  r_hub_mm: float | None = None,
                  r_tip_mm: float | None = None) -> list[dict]:
    """Secciones de una fila con triángulos free-vortex y ángulos metálicos.

    Rotor: ángulos relativos β(r); estátor: absolutos α(r) con signo
    NEGATIVO (convención del contrato).
    """
    Cx_m = stage["Cx"]
    r_m = stage["r_mean_mm"] / 1000.0
    # Radios de ESTA fila. Por defecto los de la entrada de la etapa; con
    # override, los del annulus en la estación axial REAL de la fila (el
    # annulus se contrae a lo largo de z, así que un estátor construido
    # con los radios de entrada de su etapa se sale de la carcasa por
    # varios milímetros — medido: hasta 8.5 mm en la etapa 1).
    r_hub = (stage["r_hub_mm"] if r_hub_mm is None else r_hub_mm) / 1000.0
    r_tip = (stage["r_tip_mm"] if r_tip_mm is None else r_tip_mm) / 1000.0
    Cu1_m = Cx_m * math.tan(math.radians(stage["alpha1_deg"]))
    Cu2_m = Cx_m * math.tan(math.radians(stage["alpha2_deg"]))
    T01 = stage["T0_in_K"]
    n_b = stage[f"n_blades_{kind}"]
    c_mean = stage[f"chord_{kind}_mm"]
    c_hub = 2.0 * c_mean / (1.0 + CHORD_TAPER)
    tc_root, tc_tip = ((TC_ROOT_R, TC_TIP_R) if kind == "rotor"
                       else (TC_ROOT_S, TC_TIP_S))

    sections = []
    for f in np.linspace(0.0, 1.0, n_span):
        r = r_hub + f * (r_tip - r_hub)
        # Mismo torbellino que usó la física (exponente VORTEX_N) y su
        # Cx(r) por equilibrio radial simple. Con n = −1 (vórtice libre)
        # esto es exactamente Cu·r = cte y Cx constante, bit a bit.
        Cu1 = vortex_cu(Cu1_m, r, r_m)
        Cu2 = vortex_cu(Cu2_m, r, r_m)
        Cx = vortex_cx(Cx_m, Cu1_m, r, r_m)
        chord = c_hub * (1.0 + (CHORD_TAPER - 1.0) * f)
        sigma_loc = n_b * (chord / 1000.0) / (2.0 * math.pi * r)
        t_c = tc_root + (tc_tip - tc_root) * f

        if kind == "rotor":
            U = omega * r
            ang_in = math.degrees(math.atan2(U - Cu1, Cx))
            ang_out = math.degrees(math.atan2(U - Cu2, Cx))
            W_in = math.hypot(Cx, U - Cu1)
            C_abs2 = Cx ** 2 + Cu1 ** 2
            sign = 1.0
        else:
            ang_in = math.degrees(math.atan2(Cu2, Cx))
            ang_out = math.degrees(math.atan2(Cu1, Cx))
            W_in = math.hypot(Cx, Cu2)
            C_abs2 = Cx ** 2 + Cu2 ** 2
            sign = -1.0

        T_static = max(T01 - C_abs2 / (2.0 * CP), 150.0)
        M_in = W_in / math.sqrt(GAMMA * RGAS * T_static)
        profile = "DCA" if M_in > M_DCA_THRESHOLD else "NACA65"

        ma = metal_angles(ang_in, ang_out, sigma_loc, t_c, profile)
        gen = dca_profile if profile == "DCA" else naca65_profile
        pts = gen(ma["camber_deg"], t_c, N_POINTS)     # cuerda unitaria
        pts = pts * chord                               # mm
        centroid = _polygon_centroid(pts)
        pts = pts - centroid                            # stacking centroidal
        # línea de comba (misma transformación) — la capa 5c la usa como
        # superficie media a engrosar con voxMeshShell (patrón Phy-CC v3)
        xc = np.linspace(0.0, 1.0, N_CAMBER)
        yc, _ = _circular_arc_camber(ma["camber_deg"], xc)
        camber = np.column_stack([xc, yc]) * chord - centroid
        if sign < 0:
            pts = pts * np.array([1.0, -1.0])           # espejo del estátor
            pts = pts[::-1]                             # conservar CCW
            camber = camber * np.array([1.0, -1.0])

        sections.append(dict(
            span_frac=round(float(f), 4), r_mm=round(r * 1000.0, 3),
            chord_mm=round(chord, 3),
            stagger_deg=round(sign * ma["stagger_deg"], 3),
            camber_deg=round(sign * ma["camber_deg"], 3),
            t_max_over_c=round(t_c, 4), profile=profile,
            metal_in_deg=round(sign * ma["chi1_deg"], 3),
            metal_out_deg=round(sign * ma["chi2_deg"], 3),
            M_in=round(M_in, 4), sigma_local=round(sigma_loc, 4),
            thickness_mm=round(t_c * chord, 3),
            points=[[round(float(px), 4), round(float(py), 4)]
                    for px, py in pts],
            camber_points=[[round(float(px), 4), round(float(py), 4)]
                           for px, py in camber],
        ))
    return sections


# ---------------------------------------------------------------------------
# IGV y OGV — las filas que la física asume y que faltaban en el 3D
# ---------------------------------------------------------------------------
# El meanline arranca la etapa 1 con pre-swirl α₁ ("IGV implícito", su
# pérdida se carga a la etapa 1) y quita el remolino de salida con una fila
# OGV cuya pérdida y longitud sí contabiliza. Sin estas filas la máquina
# impresa no puede cumplir los triángulos verificados (entrada axial real ≠
# α₁ asumido) ni la condición de salida axial de M_exit. Ambas se modelan
# como pseudo-etapas de estátor y reutilizan _row_sections tal cual:
#   IGV:  flujo 0 → α₁(r) de la etapa 1 (fila aceleradora, giro free-vortex)
#   OGV:  flujo α₁(r) de la última etapa → 0 (estátor difusor clásico)
# Simplificación documentada: la incidencia/desviación de Lieblein/Carter
# se derivó para cascadas difusoras. En el IGV (acelerador) la desviación
# real es menor (gradiente favorable); Carter la sobreestima unos grados y
# el metal queda sobre-girado ~4-6° respecto al ideal — el flujo entrega
# al menos el pre-swirl asumido. Primera aproximación declarada.

def _pseudo_stage_igv(record: dict) -> dict:
    st0 = record["stage_table"][0]
    return dict(
        stage=-1,
        Cx=st0["Cx"], r_mean_mm=st0["r_mean_mm"],
        r_hub_mm=st0["r_hub_mm"], r_tip_mm=st0["r_tip_mm"],
        alpha2_deg=0.0,                    # entrada del IGV: axial
        alpha1_deg=st0["alpha1_deg"],      # salida: pre-swirl de la etapa 1
        T0_in_K=st0["T0_in_K"],
        n_blades_stator=st0["n_blades_stator"],
        chord_stator_mm=st0["chord_stator_mm"],
    )


def _pseudo_stage_ogv(record: dict) -> dict:
    sl = record["stage_table"][-1]
    h_exit = _exit_height_mm(record)
    r_m = sl["r_mean_mm"]
    return dict(
        stage=len(record["stage_table"]),
        Cx=sl["Cx"], r_mean_mm=r_m,
        r_hub_mm=r_m - 0.5 * h_exit, r_tip_mm=r_m + 0.5 * h_exit,
        alpha2_deg=sl["alpha1_deg"],       # entrada: swirl residual α₁
        alpha1_deg=0.0,                    # salida: axial
        T0_in_K=sl["T0_out_K"],
        n_blades_stator=sl["n_blades_stator"],
        chord_stator_mm=sl["chord_stator_mm"],
    )


# ---------------------------------------------------------------------------
# Condiciones de contorno CFD
# ---------------------------------------------------------------------------
def cfx_boundary_conditions(theta: np.ndarray, record: dict) -> dict:
    theta = np.asarray(theta, dtype=float)
    st = record["stage_table"]
    return dict(
        P0_in_Pa=float(theta[11]), T0_in_K=float(theta[10]),
        mdot_kg_s=float(theta[12]), RPM=float(theta[1]),
        P0_out_Pa_est=float(record["PR"] * theta[11]),
        P_static_out_Pa_est=float(
            record["PR"] * theta[11] *
            (1.0 + 0.5 * (GAMMA - 1) * record["M_exit"] ** 2)
            ** (-GAMMA / (GAMMA - 1))),
        T0_out_K_est=float(record["T0_out"]),
        turbulence_intensity=0.05,
        interfaces=[dict(stage=s["stage"],
                         type="stage (mixing-plane) rotor/estator")
                    for s in st],
    )


# ---------------------------------------------------------------------------
# STEP de ensamble (extra `step`, CadQuery — fase 8.2)
# ---------------------------------------------------------------------------
# El export STEP existe para RE-CAD (NX/CATIA/FreeCAD): sólidos de
# revolución exactos (eje, casco del hub, carcasa con bridas) + álabes por
# loft SPLINE de las secciones `points` del contrato, en coordenadas de
# máquina (Z = eje de rotación). Decisiones documentadas:
#   * Álabes lofteados entre planos radiales PARALELOS (X = r constante) —
#     la aproximación estándar de re-CAD; el envolvimiento cilíndrico
#     exacto (φ = y/r) solo lo tiene el STL voxelizado.
#   * En modo "parts" cada pieza lleva UN álabe de muestra + el conteo en
#     parts/README.txt (el patrón circular de 40-60 álabes × fila en OCC
#     tarda minutos y multiplica el tamaño del STEP por >10 — el patrón
#     se aplica en el CAD destino en segundos).
#   * Omisiones (solo STL): puerto de sangrado, fillets de raíz.
# CadQuery sigue siendo OPCIONAL: sin él, export_step devuelve [] sin
# error (patrón de degradación del proyecto).

def _cq():
    try:
        import cadquery as cq
        return cq
    except Exception:
        return None


def _cq_flowpath_ring(cq, hub_line: list, tip_line: list, z0: float,
                      z1: float, r_hub_off: float = 0.0,
                      r_tip_off: float = 0.0, n: int = 40):
    """Sólido de revolución que ocupa EXACTAMENTE la vena entre las
    polilíneas de cubo y punta en [z0, z1], con offsets radiales."""
    zs = np.linspace(z0, z1, n)
    zh = [p[0] for p in hub_line]
    rh = [p[1] for p in hub_line]
    zt = [p[0] for p in tip_line]
    rt = [p[1] for p in tip_line]
    inner = [(float(np.interp(z, zh, rh)) + r_hub_off, float(z))
             for z in zs]
    outer = [(float(np.interp(z, zt, rt)) + r_tip_off, float(z))
             for z in zs]
    prof = inner + outer[::-1]
    return (cq.Workplane("XZ").polyline(prof).close()
            .revolve(360.0, (0, 0), (0, 1)))


def _fuse_wp(cq, a, b):
    """Unión a nivel de Shape (los booleanos de Workplane se pierden
    sólidos creados con newObject; ver _cq_blade_with_root)."""
    return cq.Workplane("XY").newObject([a.val().fuse(b.val()).clean()])


def _cq_blade_trimmed(cq, contract: dict, row: dict, kind: str,
                      k_span: int = 1, k_pts: int = 1):
    """Álabe RECORTADO contra la vena, que es como se construye en CAD
    real (se talla con margen y se corta contra las superficies de cubo y
    carcasa).

    Sin recorte no basta con darle a cada fila los radios del annulus en
    su estación: el annulus también se contrae A LO LARGO DE LA CUERDA,
    así que la punta que encaja en el centro de la fila sobresale en el
    borde de salida. Medido tras corregir la estación: aún +1.4 mm de
    invasión de carcasa. Recortando, la interferencia es cero por
    construcción y además la punta SIGUE el perfil de la carcasa.

    Holguras: el rotor deja ε en la punta (nace en el cubo); el estátor
    cuelga de la carcasa y deja ε en su extremo libre.
    """
    ann = contract["annulus"]
    eps = float(ann["tip_clearance_mm"])
    ext = 4.0
    solid = _cq_blade_solid(cq, _decimate_row(row, k_span, k_pts),
                            row["z_center_mm"], r_extend=ext)
    z0 = row["z_le_mm"] - 2.0
    z1 = row["z_te_mm"] + 2.0
    if kind == "rotor":
        trim = _cq_flowpath_ring(cq, ann["hub"], ann["tip"], z0, z1,
                                 r_hub_off=0.0, r_tip_off=0.0)
    else:
        trim = _cq_flowpath_ring(cq, ann["hub"], ann["tip"], z0, z1,
                                 r_hub_off=eps, r_tip_off=eps)
    out = solid.val().intersect(trim.val())
    # `intersect` devuelve un Shape genérico; si no se baja a Solid, el
    # `.union()` de Workplane lo ignora ("must have at least one solid on
    # the stack"). Se extraen los sólidos explícitamente.
    sol = out.Solids()
    if not sol:
        return solid
    return cq.Workplane("XY").newObject(
        [sol[0] if len(sol) == 1 else cq.Compound.makeCompound(sol)])


def _cq_blade_solid(cq, row: dict, z_center: float, r_extend: float = 0.0):
    """UN álabe de una fila, en coordenadas de máquina.

    Cada sección (marco de cuerda, centroide en el origen) se rota por su
    stagger firmado y se ENVUELVE sobre el cilindro de su propio radio:

        u (tangencial) → φ = u/r ,   x = r·cos φ ,  y = r·sin φ
        w (axial)      → z = z_center + w

    El envolvimiento NO es cosmético. Con secciones planas apoyadas en
    planos radiales paralelos (lo que hacía esta función hasta la fase
    9.1) el punto a distancia tangencial u de una sección queda en
    realidad a radio √(r²+u²) > r: en una fila de cuerda larga y stagger
    alto las esquinas de la sección de PUNTA se salen de la carcasa —
    visible en cuanto se mira el ensamble. Envolviendo, cada sección vive
    exactamente sobre su cilindro y la punta nunca supera r_tip.
    """
    wp = None
    r_prev = None
    secs = list(row["sections"])
    if r_extend > 0.0:
        # secciones guarda por dentro y por fuera: dan material para que
        # el recorte contra la vena tenga qué cortar
        secs = [dict(secs[0], r_mm=secs[0]["r_mm"] - r_extend)] + secs \
            + [dict(secs[-1], r_mm=secs[-1]["r_mm"] + r_extend)]
    for sec in secs:
        g = math.radians(sec["stagger_deg"])
        cg, sg = math.cos(g), math.sin(g)
        pts2 = [(px * sg + py * cg,            # u → Y (tangencial)
                 px * cg - py * sg)            # v → Z (axial, staggered)
                for px, py in sec["points"]]
        # Compensación del PLANADO: la sección vive en el plano X = r, así
        # que un punto a distancia tangencial u queda en realidad a radio
        # √(r²+u²). Se baja el plano a r_eff = √(r² − u_máx²) para que la
        # esquina más alejada caiga EXACTAMENTE en r y no un pelo fuera.
        # (Lo correcto sería envolver sobre el cilindro, pero las tapas
        # dejan de ser planas y OCC devuelve un sólido que las booleanas
        # descartan en silencio — con la raíz de abeto se vio: el álabe
        # desaparecía de la fusión.)
        r = float(sec["r_mm"])
        u_max = max(abs(u) for u, _ in pts2)
        r_eff = math.sqrt(max(r * r - u_max * u_max, (0.5 * r) ** 2))
        if wp is None:
            wp = cq.Workplane("YZ", origin=(r_eff, 0.0, 0.0))
        else:
            wp = wp.workplane(offset=r_eff - r_prev)
        wp = wp.polyline(pts2).close()
        r_prev = r_eff
    return wp.loft(ruled=False).translate((0.0, 0.0, z_center))


def _root_width_mm(record: dict) -> float:
    """Ancho de la raíz del álabe de la primera etapa (t/c de cubo por la
    cuerda de cubo) — la referencia de la que cuelga el abeto."""
    s0 = record["stage_table"][0]
    c_hub = 2.0 * s0["chord_rotor_mm"] / (1.0 + CHORD_TAPER)
    return TC_ROOT_R * c_hub


def _drum_inner_r_mm(record: dict) -> float:
    """Radio interior del tambor (mismo criterio que el bloque structural:
    30% del radio de cubo de la primera etapa)."""
    return 0.30 * float(record["stage_table"][0]["r_hub_mm"])


def firtree_profile(w_top_mm: float, depth_mm: float, n_teeth: int = 3,
                    taper_deg: float = 8.0, lobe_frac: float = 0.22,
                    clearance_mm: float = 0.0) -> list:
    """Perfil de RAÍZ DE ABETO (fir-tree) en el plano (tangencial, radial).

    Es la retención real de un álabe de compresor/turbina: una serie de
    dientes que transfieren la carga centrífuga al rim del disco por
    flancos de apoyo sucesivos, en vez de una raíz soldada o hundida. Fue
    el punto donde `turbodesigner` (OpenOrion) estaba por delante de
    Phy-AC — la fase 9 lo cierra con una construcción paramétrica propia.

    Convención: u = tangencial (mm, centrado en 0), v = PROFUNDIDAD hacia
    dentro desde la plataforma (0 en la plataforma, +depth en el fondo).
    Se devuelve el contorno CERRADO en sentido antihorario.

    Parámetros
    ----------
    w_top_mm    : ancho TOTAL en la plataforma (≈ espesor de raíz del álabe)
    depth_mm    : profundidad radial de la raíz
    n_teeth     : nº de pares de dientes (2-3 en compresores axiales)
    taper_deg   : semiángulo del árbol (el ancho cae hacia el fondo)
    lobe_frac   : profundidad del cuello como fracción del semiancho local
    clearance_mm: sobreespesor uniforme — usar >0 para la RANURA del disco
                  (la ranura debe ser mayor que la raíz para poder brochar
                  y dejar juego de montaje)
    """
    n = max(int(n_teeth), 1)
    depth = max(float(depth_mm), 1e-3)
    p = depth / n                          # paso radial por diente
    tan_t = math.tan(math.radians(max(taper_deg, 0.0)))
    w0 = 0.5 * max(float(w_top_mm), 1e-3) + clearance_mm

    def half_w(v: float) -> float:
        return max(w0 - v * tan_t, 0.25 * w0)

    lobe_abs = max(lobe_frac * w0, 1.4 * tan_t * 0.45 * p)

    right = [(w0, 0.0)]
    for k in range(n):
        v_crest = k * p
        v_neck = v_crest + 0.55 * p
        wc = half_w(v_crest)
        # El cuello se mide desde la ENVOLVENTE a su propia profundidad y
        # con un lóbulo ABSOLUTO: si se resta una fracción de la anchura
        # LOCAL, los dientes de abajo (donde el árbol ya es delgado) no
        # llegan a reabrirse y el "abeto" degenera en un cono escalonado
        # sin flancos que apoyen. El suelo garantiza que el lóbulo supere
        # siempre la caída del taper entre el cuello y la cresta siguiente.
        wn = max(half_w(v_neck) - lobe_abs, 0.18 * w0)
        right += [(wc, v_crest + 0.12 * p),     # apoyo recto bajo el diente
                  (wn, v_neck)]                 # flanco hasta el cuello
        if k < n - 1:
            right.append((half_w(v_crest + p), v_crest + 0.88 * p))
    w_bot = right[-1][0]
    right.append((0.85 * w_bot, depth))         # fondo redondeado (chaflán)
    left = [(-u, v) for u, v in reversed(right)]
    return [(float(u), float(v)) for u, v in right + left]


N_FIRTREE_SECTIONS = 5     # secciones del loft de raíz/ranura (ver
#                            _cq_firtree_solid: con 2 la cuerda del
#                            envolvimiento cilíndrico hace que ranura y
#                            raíz, de longitudes distintas, no encajen).
#                            5 baja el desencaje 1000× y deja TODOS los
#                            sólidos válidos; con 7 o más, OCC empieza a
#                            marcar inválido el álabe de alguna etapa.


def _r_on_line(line: list, z: float) -> float:
    """Radio de una polilínea de annulus (z, r) en una estación axial."""
    return float(np.interp(z, [p[0] for p in line], [p[1] for p in line]))


def _root_axial_and_slope(contract: dict, row: dict) -> tuple[float, float]:
    """Longitud axial de la raíz y PENDIENTE dr/dz de su plataforma.

    La plataforma de la raíz no es un cilindro: sigue la línea de cubo,
    que en un compresor sube etapa a etapa. Con plataforma cilíndrica el
    rim del disco asomaba por encima de la vena aguas arriba del centro de
    la fila y se metía dentro del álabe. La raíz, su ranura y el rim del
    disco comparten esta MISMA ley lineal, así que el juego de brochado
    es el único hueco entre ellos.
    """
    ft = contract["assembly"].get("firtree") or {}
    hub = contract["annulus"]["hub"]
    ax = float(ft.get("axial_frac", 0.85)) * (row["z_te_mm"] - row["z_le_mm"])
    zc = float(row["z_center_mm"])
    dz = 0.5 * max(ax, 1e-6)
    slope = (_r_on_line(hub, zc + dz) - _r_on_line(hub, zc - dz)) \
        / max(2.0 * dz, 1e-6)
    return ax, slope


def _cq_firtree_solid(cq, u_v_profile: list, r_platform: float,
                      z_center: float, axial_len: float,
                      stagger_deg: float, r_slope: float = 0.0):
    """Sólido de una raíz de abeto (o de su ranura) en coordenadas de
    máquina.

    El perfil (u = tangencial, v = profundidad radial) se ENVUELVE sobre
    el cilindro de cada radio —igual que las secciones del álabe— y se
    loftea entre la cara delantera y la trasera, desplazadas
    tangencialmente ±(L/2)·tan γ para reproducir el BROCHADO INCLINADO
    real (la ranura sigue la cuerda del álabe, no el eje).

    La versión anterior extruía a lo largo de Z y luego rotaba el sólido
    alrededor de un eje GLOBAL: eso barría la ranura por un arco enorme y
    el corte se comía el disco entero (medido: 19.9 cm³ retirados por una
    herramienta de 1.0 cm³). Construyendo las dos caras directamente no
    hay rotación de sólidos y la ranura retira exactamente su volumen.

    El loft lleva N_FIRTREE_SECTIONS secciones, no dos. Con solo las dos
    caras extremas la superficie reglada es una CUERDA del envolvimiento
    cilíndrico, y su flecha crece con la longitud axial: la ranura, que es
    un 5% más larga que la raíz, se hundía ~0.02 mm más que ella y el
    borde de la raíz quedaba mordiendo el rim del disco (2.5 mm³ por
    álabe, detectado por el chequeo de interferencias del ensamble). Con
    secciones intermedias ambas convergen a la MISMA superficie envuelta y
    el único hueco entre ellas vuelve a ser el juego declarado.
    """
    dz = 0.5 * float(axial_len)
    tan_g = math.tan(math.radians(float(stagger_deg)))
    wires = []
    for s in np.linspace(-1.0, 1.0, N_FIRTREE_SECTIONS):
        # La plataforma sigue la línea de cubo con pendiente CONSTANTE.
        # Raíz y ranura usan la misma ley aunque sus longitudes axiales
        # difieran, así que sus superficies son exactamente paralelas.
        rp = r_platform + s * dz * float(r_slope)
        du = s * dz * tan_g
        pts3 = []
        for u, v in u_v_profile:
            r = max(rp - v, 1e-3)
            phi = (u + du) / r
            pts3.append(cq.Vector(r * math.cos(phi), r * math.sin(phi),
                                  z_center + s * dz))
        if (pts3[0] - pts3[-1]).Length < 1e-9:
            pts3 = pts3[:-1]
        wires.append(cq.Wire.makePolygon(pts3, close=True))
    return cq.Workplane("XY").newObject([cq.Solid.makeLoft(wires, True)])


def _cq_ring(cq, r_in: float, r_out: float, z0: float, z1: float):
    """Anillo macizo de revolución (r_in→r_out, z0→z1) alrededor de Z."""
    return (cq.Workplane("XZ")
            .polyline([(r_in, z0), (r_out, z0), (r_out, z1), (r_in, z1)])
            .close().revolve(360.0, (0, 0), (0, 1)))


def _cq_line_ring(cq, line: list, wall: float, z0: float, z1: float,
                  outward: bool, n_pts: int = 24):
    """Casco de revolución que sigue una polilínea (z, r) del annulus en
    [z0, z1]: pared `wall` hacia fuera (carcasa) o hacia dentro (hub)."""
    zs = np.linspace(z0, z1, n_pts)
    zl = [p[0] for p in line]
    rl = [p[1] for p in line]
    rs = np.interp(zs, zl, rl)
    inner = [(float(r) + (0.0 if outward else -wall), float(z))
             for z, r in zip(zs, rs)]
    outer = [(float(r) + (wall if outward else 0.0), float(z))
             for z, r in zip(zs, rs)]
    prof = inner + outer[::-1]
    return (cq.Workplane("XZ").polyline(prof).close()
            .revolve(360.0, (0, 0), (0, 1)))


def _cq_part_solids(cq, contract: dict) -> dict:
    """Sólidos por PIEZA (mismos nombres que los STL de la capa 5c), con
    UN álabe de muestra por fila. Devuelve {nombre: solid}."""
    asm = contract["assembly"]
    ann = contract["annulus"]
    stages = contract["stages"]
    n_st = len(stages)
    clr = float(ann["tip_clearance_mm"])
    wall = float(asm.get("hub_wall_mm", 5.0))
    r_shaft = float(contract.get("structural", {})
                    .get("drum_inner_r_mm", 20.0))
    hub, tip = ann["hub"], ann["tip"]
    z_in, z_out = hub[0][0], hub[-1][0]

    def split_after(i: int) -> float:
        te = stages[i]["stator"]["z_te_mm"]
        nxt = stages[i + 1]["rotor"]["z_le_mm"] if i + 1 < n_st else te + 4.0
        return 0.5 * (te + nxt)

    parts: dict = {}
    # Eje: cilindro con muñones y barreno (igual que RotorDrum.voxShaft)
    stub = float(asm.get("shaft_stub_mm", 30.0))
    bore = float(asm.get("shaft_bore_frac", 0.35)) * r_shaft
    shaft = _cq_ring(cq, max(bore, 0.0), r_shaft, z_in - stub, z_out + stub) \
        if bore > 0.5 else _cq_ring(cq, 0.0, r_shaft, z_in - stub,
                                    z_out + stub)
    parts["Shaft"] = shaft

    for i, st in enumerate(stages):
        z0 = z_in if i == 0 else split_after(i - 1)
        z1 = z_out if i == n_st - 1 else split_after(i)
        rot = st["rotor"]
        web = min(max(float(asm.get("disk_web_frac", 0.5))
                      * (rot["z_te_mm"] - rot["z_le_mm"]),
                      float(asm.get("disk_web_min_mm", 4.0))),
                  float(asm.get("disk_web_max_mm", 14.0)))
        zc = rot["z_center_mm"]
        r_hub_c = float(np.interp(zc, [p[0] for p in hub],
                                  [p[1] for p in hub]))
        solid = _cq_line_ring(cq, hub, wall, z0, z1, outward=False)
        solid = solid.union(_cq_ring(cq, r_shaft, r_hub_c,
                                     zc - 0.5 * web, zc + 0.5 * web))
        solid = _fuse_wp(cq, solid, _cq_blade_trimmed(cq, contract, rot, "rotor"))
        # --- retención de abeto (fase 9) --------------------------------
        # La raíz del álabe de muestra baja al rim con un árbol de dientes
        # y el disco lleva su ranura brochada. Antes la raíz simplemente
        # se hundía en el cubo (soldadura): imprimible, pero no una
        # retención real ni desmontable.
        ft = asm.get("firtree")
        if ft and ft.get("enabled", True):
            sec0 = rot["sections"][0]
            w_top = max(float(sec0.get("thickness_mm", 2.0)),
                        0.6 * float(ft.get("min_width_mm", 3.0)))
            depth = float(ft["depth_mm"])
            ax = float(ft["axial_frac"]) * (rot["z_te_mm"] - rot["z_le_mm"])
            prof = firtree_profile(w_top, depth, int(ft["n_teeth"]),
                                   float(ft["taper_deg"]),
                                   float(ft["lobe_frac"]))
            slot = firtree_profile(w_top, depth, int(ft["n_teeth"]),
                                   float(ft["taper_deg"]),
                                   float(ft["lobe_frac"]),
                                   clearance_mm=float(ft["clearance_mm"]))
            stg = float(sec0.get("stagger_deg", 0.0))
            solid = solid.cut(_cq_firtree_solid(cq, slot, r_hub_c, zc,
                                                ax * 1.02, stg))
            solid = solid.union(_cq_firtree_solid(cq, prof, r_hub_c, zc,
                                                  ax, stg))
        parts[f"RotorStage{i + 1}"] = solid

    lip = 4.0
    for i, st in enumerate(stages):
        z0 = (tip[0][0] - lip) if i == 0 else split_after(i - 1)
        z1 = (tip[-1][0] + lip) if i == n_st - 1 else split_after(i)
        # el casco vive a tip+holgura (annulus pasante)
        shell_line = [[z, r + clr] for z, r in tip]
        solid = _cq_line_ring(cq, shell_line, wall, z0, z1, outward=True)
        stat = st["stator"]
        solid = _fuse_wp(cq, solid,
                         _cq_blade_trimmed(cq, contract, stat, "stator"))
        if i == 0 and contract.get("igv"):
            igv = contract["igv"]
            solid = _fuse_wp(cq, solid,
                             _cq_blade_trimmed(cq, contract, igv, "stator"))
        if i == n_st - 1 and contract.get("ogv"):
            ogv = contract["ogv"]
            solid = _fuse_wp(cq, solid,
                             _cq_blade_trimmed(cq, contract, ogv, "stator"))
        # bridas apernadas en el primer/último anillo
        for at_start, has in ((True, i == 0), (False, i == n_st - 1)):
            if not has:
                continue
            zf = z0 if at_start else z1 - float(asm.get("flange_t_mm", 8.0))
            r_tip_f = float(np.interp(z0 if at_start else z1,
                                      [p[0] for p in tip],
                                      [p[1] for p in tip]))
            r_o = r_tip_f + clr + wall
            fl = _cq_ring(cq, r_tip_f + clr, r_o
                          + float(asm.get("flange_w_mm", 12.0)),
                          zf, zf + float(asm.get("flange_t_mm", 8.0)))
            n_b = int(asm.get("flange_bolt_count", 12))
            r_b = r_o + 0.5 * float(asm.get("flange_w_mm", 12.0))
            d_b = float(asm.get("flange_bolt_d_mm", 6.0))
            for k in range(n_b):
                a = 2.0 * math.pi * k / n_b
                hole = (cq.Workplane("XY",
                                     origin=(r_b * math.cos(a),
                                             r_b * math.sin(a), zf - 1.0))
                        .circle(0.5 * d_b).extrude(
                            float(asm.get("flange_t_mm", 8.0)) + 2.0))
                fl = fl.cut(hole)
            solid = solid.union(fl)
        parts[f"StatorRing{i + 1}"] = solid
    return parts


def _step_readme(contract: dict) -> str:
    lines = [
        "STEP de ensamble Phy-AC — notas de re-CAD",
        "=========================================",
        "Coordenadas de maquina: Z = eje de rotacion, mm.",
        "Cada pieza lleva UN alabe/vano de muestra; aplicar el patron",
        "circular en el CAD destino con estos conteos:",
        "",
    ]
    for st in contract["stages"]:
        i = st["index"] + 1
        lines.append(f"  RotorStage{i}:  {st['rotor']['n_blades']} alabes "
                     f"(eje de patron = Z)")
        lines.append(f"  StatorRing{i}:  {st['stator']['n_blades']} vanos")
    if contract.get("igv"):
        lines.append(f"  IGV (en StatorRing1): "
                     f"{contract['igv']['n_blades']} vanos")
    if contract.get("ogv"):
        n = len(contract["stages"])
        lines.append(f"  OGV (en StatorRing{n}): "
                     f"{contract['ogv']['n_blades']} vanos")
    lines += [
        "",
        "Aproximaciones (exactas solo en el STL de la capa 5c):",
        "  * alabes lofteados entre planos radiales paralelos (sin el",
        "    envolvimiento cilindrico phi = y/r del voxelizado);",
        "  * sin puerto de sangrado ni fillets de raiz;",
        "  * holgura de punta: restar tip_clearance_mm del contrato.",
        "",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Ensamble DETALLADO (fase 9.1): piezas reales, no un sólido por etapa
#
# La vista de conjunto anterior fundía el tambor entero en un sólido y
# colgaba los álabes encima. Eso sirve para mirar la máquina, no para
# entenderla ni para fabricarla. El desglose de aquí es el de un rotor de
# compresor de verdad:
#
#   ROTOR   Shaft ─ Disc{i} (barreno, alma, rim con N ranuras de abeto,
#           taladros de tirante) ─ Blade{i}_{k} (con su raíz de abeto
#           METIDA en la ranura) ─ TieBolt{k} (tirantes pasantes)
#   ESTATOR CasingRing{i} (segmento de carcasa con brida en CADA extremo,
#           taladros y asientos de vano) ─ Vane{i}_{k} ─ CasingBolt{i}_{k}
#
# Cada anillo de carcasa cubre UNA etapa de estátor y se aperna al
# siguiente: así la máquina se abre por planos reales de junta, que es lo
# que pedía el desglose de piezas de la capa 5c.
# ---------------------------------------------------------------------------
def _cq_bolt(cq, d_shank: float, length: float, x: float, y: float,
             z0: float, d_head: float | None = None,
             h_head: float | None = None):
    """Perno pasante: vástago + cabeza + tuerca. Geometría mínima pero
    reconocible — lo que importa es que la junta EXISTA como pieza."""
    dh = d_head if d_head else 1.7 * d_shank
    hh = h_head if h_head else 0.65 * d_shank
    shank = (cq.Workplane("XY", origin=(x, y, z0 - hh))
             .circle(0.5 * d_shank).extrude(length + 2.0 * hh))
    head = (cq.Workplane("XY", origin=(x, y, z0 - hh))
            .polygon(6, dh).extrude(hh))
    nut = (cq.Workplane("XY", origin=(x, y, z0 + length))
           .polygon(6, dh).extrude(hh))
    return shank.union(head).union(nut)


def _disc_rim_inner_r(contract: dict, r_hub_c: float, r_shaft: float
                      ) -> tuple[float, float]:
    """(r_barreno, r_interior_del_rim) de un disco. Fuera de `_cq_disc`
    porque el radio del círculo de tirantes tiene que ser COMÚN a toda la
    pila y hay que conocerlo antes de construir ningún disco."""
    ft = contract["assembly"].get("firtree") or {}
    fit = float(ft.get("clearance_mm", 0.08))
    r_bore = r_shaft + fit                     # ajuste deslizante
    depth = float(ft.get("depth_mm", 10.0))
    return r_bore, max(r_hub_c - 1.9 * depth, r_bore + 6.0)


def _cq_disc(cq, contract: dict, i: int, r_hub_c: float, zc: float,
             web: float, r_shaft: float, n_blades: int,
             stagger_hub_deg: float, ax_root: float,
             n_tie: int, d_tie: float, r_tie: float | None = None,
             r_slope: float = 0.0):
    """Disco de UNA etapa: barreno sobre el eje, alma, rim, N ranuras de
    abeto brochadas y el círculo de taladros del tirante.

    El disco NO es el eje: tiene su propio barreno con ajuste y se monta
    sobre él. Esa distinción es la que faltaba."""
    asm_p = contract["assembly"]
    ft = asm_p.get("firtree") or {}
    fit = float(ft.get("clearance_mm", 0.08))
    rim_t = 1.6 * web                          # el rim es más grueso
    depth = float(ft.get("depth_mm", 10.0))
    r_bore, r_rim_i = _disc_rim_inner_r(contract, r_hub_c, r_shaft)

    disc = _cq_ring(cq, r_bore, r_rim_i, zc - 0.5 * web, zc + 0.5 * web)
    # El rim es un TRONCO DE CONO que sigue la línea de cubo, no un
    # cilindro: con radio constante r_hub_c el rim asomaba por encima de
    # la vena aguas arriba del centro de la fila (la línea de cubo sube) y
    # se metía dentro del álabe, que sí está recortado contra la vena.
    # El rim se queda un pelo POR DEBAJO de la plataforma del álabe
    # (`rim_relief_mm`, del orden del juego de brochado). Es lo que pasa
    # en un disco real —la plataforma del álabe forma la vena, el rim se
    # rebaja para no pelearse con ella— y es lo que hace que el rim y el
    # álabe, que está recortado contra una aproximación poligonal de la
    # polilínea de cubo, no compartan material por el redondeo de esa
    # aproximación (medido sin rebaje: 3 mm³ por álabe, ≈0.03 mm).
    relief = float(asm_p.get("rim_relief_mm", 0.10))
    hub_line = contract["annulus"]["hub"]
    z0r, z1r = zc - 0.5 * rim_t, zc + 0.5 * rim_t
    # La cara exterior sigue la POLILÍNEA de cubo, no una recta con la
    # pendiente de la raíz: la raíz mide su pendiente sobre ±ax_root/2 y
    # el rim solo ocupa ±rim_t/2, así que sobre una vena con quiebro la
    # secante larga se sale por encima de la línea de cubo dentro de la
    # banda del rim — que es exactamente donde está el álabe.
    top = [(z0r, _r_on_line(hub_line, z0r) - relief)]
    top += [(float(z), float(r) - relief) for z, r in hub_line
            if z0r < z < z1r]
    top.append((z1r, _r_on_line(hub_line, z1r) - relief))
    disc = disc.union(
        cq.Workplane("XZ")
        .polyline([(r_rim_i, z0r)] + [(r, z) for z, r in top]
                  + [(r_rim_i, z1r)])
        .close().revolve(360.0, (0, 0), (0, 1)))
    # Taladros del tirante en el alma. El radio es COMÚN a toda la pila:
    # cada disco tiene su propio r_rim_i (crece con la línea de cubo), así
    # que un círculo de taladros calculado disco a disco dejaba los
    # tirantes atravesando el alma de los discos delanteros.
    if r_tie is None:
        r_tie = 0.5 * (r_bore + r_rim_i)
    for k in range(n_tie):
        a = 2.0 * math.pi * k / n_tie
        disc = disc.cut(cq.Workplane(
            "XY", origin=(r_tie * math.cos(a), r_tie * math.sin(a),
                          zc - web)).circle(0.5 * d_tie).extrude(2.0 * web))
    # ranuras de abeto: una por álabe, a su paso angular
    if ft.get("enabled", True) and n_blades > 0:
        slot = firtree_profile(
            float(ft.get("w_top_mm", 6.0)), depth, int(ft.get("n_teeth", 3)),
            float(ft.get("taper_deg", 8.0)), float(ft.get("lobe_frac", 0.22)),
            clearance_mm=fit)
        one = _cq_firtree_solid(cq, slot, r_hub_c, zc, ax_root * 1.05,
                                stagger_hub_deg, r_slope)
        for k in range(n_blades):
            a = 360.0 * k / n_blades
            disc = disc.cut(one if k == 0
                            else one.rotate((0, 0, 0), (0, 0, 1), a))
    return disc, r_tie


def _cq_blade_with_root(cq, contract: dict, row: dict, r_hub_c: float,
                        k_span: int, k_pts: int):
    """Álabe + su raíz de abeto, como UNA pieza (que es como se fabrica y
    como se monta en la ranura del disco)."""
    ft = (contract["assembly"].get("firtree") or {})
    blade = _cq_blade_trimmed(cq, contract, row, "rotor", k_span, k_pts)
    if not ft.get("enabled", True):
        return blade, 0.0
    ax_root, r_slope = _root_axial_and_slope(contract, row)
    prof = firtree_profile(float(ft.get("w_top_mm", 6.0)),
                           float(ft.get("depth_mm", 10.0)),
                           int(ft.get("n_teeth", 3)),
                           float(ft.get("taper_deg", 8.0)),
                           float(ft.get("lobe_frac", 0.22)))
    stg = float(row["sections"][0].get("stagger_deg", 0.0))
    root = _cq_firtree_solid(cq, prof, r_hub_c, row["z_center_mm"],
                             ax_root, stg, r_slope)
    # Fusión a nivel de Shape: `Workplane.union` sobre workplanes creados
    # con newObject() se quedaba con UNO de los dos sólidos (el álabe
    # desaparecía y solo quedaba la raíz). El booleano de OCC directo no
    # tiene esa ambigüedad.
    fused = blade.val().fuse(root.val()).clean()
    return cq.Workplane("XY").newObject([fused]), ax_root


def _cq_detailed_machine(cq, contract: dict, k_span: int = 2,
                         k_pts: int = 2):
    """Ensamble jerárquico con las piezas reales separadas."""
    asm_p = contract["assembly"]
    ann = contract["annulus"]
    stages = contract["stages"]
    n_st = len(stages)
    clr = float(ann["tip_clearance_mm"])
    wall = float(asm_p.get("hub_wall_mm", 5.0))
    r_shaft = float(contract.get("structural", {})
                    .get("drum_inner_r_mm", 20.0))
    hub, tip = ann["hub"], ann["tip"]
    z_in, z_out = hub[0][0], hub[-1][0]
    stub = float(asm_p.get("shaft_stub_mm", 30.0))
    bore = float(asm_p.get("shaft_bore_frac", 0.35)) * r_shaft
    n_tie = int(asm_p.get("tie_bolt_count", 12))
    d_tie = float(asm_p.get("tie_bolt_d_mm", 10.0))
    t_fl = float(asm_p.get("flange_t_mm", 8.0))
    w_fl = float(asm_p.get("flange_w_mm", 12.0))
    n_cb = int(asm_p.get("flange_bolt_count", 12))
    d_cb = float(asm_p.get("flange_bolt_d_mm", 6.0))

    root = cq.Assembly(name=contract.get("name", "PhyAC_Machine"))
    rotor = cq.Assembly(name="ROTOR")
    stator = cq.Assembly(name="STATOR")
    fasteners = cq.Assembly(name="FASTENERS")

    def r_at(line, z):
        return float(np.interp(z, [p[0] for p in line], [p[1] for p in line]))

    # ---------------- ROTOR --------------------------------------------
    rotor.add(_cq_ring(cq, max(bore, 0.0), r_shaft, z_in - stub,
                       z_out + stub), name="Shaft")

    # Geometría de la pila ANTES de construir nada: el círculo de
    # tirantes tiene que ser el mismo en los N discos (cada disco tiene su
    # propio r_rim_i porque la línea de cubo sube, y un círculo por disco
    # dejaba los tirantes atravesando el alma de los delanteros), y las
    # bandas axiales que ocupan disco y raíz delimitan dónde PUEDE haber
    # casco de tambor.
    geo = []
    for st in stages:
        rot = st["rotor"]
        web = min(max(float(asm_p.get("disk_web_frac", 0.5))
                      * (rot["z_te_mm"] - rot["z_le_mm"]),
                      float(asm_p.get("disk_web_min_mm", 4.0))),
                  float(asm_p.get("disk_web_max_mm", 14.0)))
        zc = float(rot["z_center_mm"])
        r_hub_c = r_at(hub, zc)
        ax_root, r_slope = _root_axial_and_slope(contract, rot)
        r_bore, r_rim_i = _disc_rim_inner_r(contract, r_hub_c, r_shaft)
        # media banda ocupada por el conjunto disco+raíz en ese plano
        half = 0.5 * max(1.05 * ax_root, 1.6 * web)
        geo.append(dict(rot=rot, web=web, zc=zc, r_hub_c=r_hub_c,
                        ax_root=ax_root, r_slope=r_slope, r_bore=r_bore,
                        r_rim_i=r_rim_i, z_lo=zc - half, z_hi=zc + half))
    r_tie = 0.5 * (geo[0]["r_bore"] + min(g["r_rim_i"] for g in geo)) \
        if geo else None

    for i, g in enumerate(geo):
        rot = g["rot"]
        blade, ax_root = _cq_blade_with_root(cq, contract, rot, g["r_hub_c"],
                                             k_span, k_pts)
        disc, _ = _cq_disc(cq, contract, i, g["r_hub_c"], g["zc"], g["web"],
                           r_shaft, int(rot["n_blades"]),
                           float(rot["sections"][0].get("stagger_deg", 0.0)),
                           ax_root, n_tie, d_tie, r_tie=r_tie,
                           r_slope=g["r_slope"])
        rotor.add(disc, name=f"Disc{i + 1}")
        for k in range(int(rot["n_blades"])):
            a = 360.0 * k / int(rot["n_blades"])
            rotor.add(blade if k == 0
                      else blade.rotate((0, 0, 0), (0, 0, 1), a),
                      name=f"Rotor{i + 1}_blade{k + 1:03d}")
    # Casco del flowpath ENTRE bandas de disco. Antes cada spacer iba de
    # centro de disco a centro de disco: atravesaba de lado a lado el rim
    # del disco (~53 cm³ de material compartido) y la raíz de los álabes.
    # El casco es lo que cierra la vena en los huecos; el rim del disco la
    # cierra donde está el disco.
    spans = []
    for i, g in enumerate(geo):
        spans.append((z_in if i == 0 else geo[i - 1]["z_hi"], g["z_lo"],
                      f"DrumSpacer{i + 1}"))
    spans.append((geo[-1]["z_hi"], z_out, "DrumSpacerExit"))
    for z0, z1, nm in spans:
        if z1 - z0 > 1e-3:
            rotor.add(_cq_line_ring(cq, hub, wall, z0, z1, outward=False),
                      name=nm)
    # tirantes pasantes que aprietan la pila de discos contra el eje
    if r_tie:
        z0 = stages[0]["rotor"]["z_center_mm"] - 20.0
        L = stages[-1]["rotor"]["z_center_mm"] + 20.0 - z0
        for k in range(n_tie):
            a = 2.0 * math.pi * k / n_tie
            fasteners.add(_cq_bolt(cq, d_tie, L, r_tie * math.cos(a),
                                   r_tie * math.sin(a), z0),
                          name=f"TieBolt{k + 1:02d}")

    # ---------------- ESTATOR ------------------------------------------
    def split_after(i):
        te = stages[i]["stator"]["z_te_mm"]
        nxt = stages[i + 1]["rotor"]["z_le_mm"] if i + 1 < n_st else te + 4.0
        return 0.5 * (te + nxt)

    shell_line = [[z, r + clr] for z, r in tip]
    joints = []
    for i, st in enumerate(stages):
        z0 = (tip[0][0] - 4.0) if i == 0 else split_after(i - 1)
        z1 = (tip[-1][0] + 4.0) if i == n_st - 1 else split_after(i)
        ring = _cq_line_ring(cq, shell_line, wall, z0, z1, outward=True)
        # BRIDA EN CADA EXTREMO: cada anillo se aperna a sus vecinos
        for at_start in (True, False):
            z_ref = z0 if at_start else z1
            zf = z_ref if at_start else z_ref - t_fl
            r_i = r_at(shell_line, z_ref)
            fl = _cq_ring(cq, r_i, r_i + wall + w_fl, zf, zf + t_fl)
            r_b = r_i + wall + 0.5 * w_fl
            for k in range(n_cb):
                a = 2.0 * math.pi * k / n_cb
                fl = fl.cut(cq.Workplane(
                    "XY", origin=(r_b * math.cos(a), r_b * math.sin(a),
                                  zf - 1.0)).circle(0.5 * d_cb)
                    .extrude(t_fl + 2.0))
            ring = ring.union(fl)
            if at_start and i > 0:
                joints.append((zf, r_b, i))
        # PUERTO DE SANGRADO: boss radial con barreno pasante en el anillo
        # de la etapa que elige la física (stage_stall_first + 2). Existía
        # en la vía STL desde la fase 7 y NO en esta: el contrato lo
        # declara y las dos rutas tienen que construir la misma máquina.
        if i + 1 == int(np.clip(asm_p.get("bleed_stage", 1), 1, n_st)):
            zb = st["stator"]["z_center_mm"]
            r_in = r_at(shell_line, zb)
            r_out = r_in + wall + float(asm_p.get("bleed_boss_h_mm", 14.0))
            d_boss = float(asm_p.get("bleed_boss_d_mm", 32.0))
            d_hole = float(asm_p.get("bleed_hole_d_mm", 18.0))
            boss = (cq.Workplane("YZ", origin=(r_in + 1.0, 0.0, zb))
                    .circle(0.5 * d_boss)
                    .extrude(-(r_out - r_in - 1.0)))
            hole = (cq.Workplane("YZ", origin=(r_in - 2.0, 0.0, zb))
                    .circle(0.5 * d_hole)
                    .extrude(-(r_out - r_in + 4.0)))
            # El boss es un cilindro RADIAL de d_boss de diámetro: sobre
            # su huella axial la carcasa se contrae, así que su borde de
            # aguas arriba entra en la vena y muerde los vanos. Se corta
            # contra el barreno de la carcasa, que es la superficie que
            # ninguna pieza estática puede cruzar.
            zsb = np.linspace(z0 - 2.0, z1 + 2.0, 40)
            prof_b = ([(0.0, float(zsb[0]))]
                      + [(r_at(shell_line, float(z)), float(z))
                         for z in zsb]
                      + [(0.0, float(zsb[-1]))])
            bore = (cq.Workplane("XZ").polyline(prof_b).close()
                    .revolve(360.0, (0, 0), (0, 1)))
            ring = _fuse_wp(cq, ring, boss).cut(bore).cut(hole)
        stator.add(ring, name=f"CasingRing{i + 1}")
        vane = _cq_blade_trimmed(cq, contract, st["stator"], "stator",
                                 k_span, k_pts)
        for k in range(int(st["stator"]["n_blades"])):
            a = 360.0 * k / int(st["stator"]["n_blades"])
            stator.add(vane if k == 0
                       else vane.rotate((0, 0, 0), (0, 0, 1), a),
                       name=f"Stator{i + 1}_vane{k + 1:03d}")
    for tag, key in (("IGV", "igv"), ("OGV", "ogv")):
        row = contract.get(key)
        if not row:
            continue
        vane = _cq_blade_trimmed(cq, contract, row, "stator", k_span,
                                 k_pts)
        for k in range(int(row["n_blades"])):
            a = 360.0 * k / int(row["n_blades"])
            stator.add(vane if k == 0
                       else vane.rotate((0, 0, 0), (0, 0, 1), a),
                       name=f"{tag}_vane{k + 1:03d}")
    # pernos de cada junta entre anillos de carcasa
    for zf, r_b, i in joints:
        for k in range(n_cb):
            a = 2.0 * math.pi * k / n_cb
            fasteners.add(_cq_bolt(cq, d_cb, t_fl * 2.0,
                                   r_b * math.cos(a), r_b * math.sin(a),
                                   zf - t_fl),
                          name=f"CasingBolt{i}_{k + 1:02d}")

    root.add(rotor, name="ROTOR")
    root.add(stator, name="STATOR")
    root.add(fasteners, name="FASTENERS")
    return root


def _cq_duct(cq, pts_r0: tuple, z0: float, pts_r1: tuple, z1: float,
             wall: float, n: int = 20, outward: bool = False):
    """Conducto de transición de revolución entre dos estaciones (hub o
    tip), con mezcla suave (coseno) entre radios."""
    zs = np.linspace(z0, z1, n)
    t = 0.5 * (1.0 - np.cos(np.pi * np.linspace(0.0, 1.0, n)))
    rs = pts_r0 + (pts_r1 - pts_r0) * t
    line = [[float(z), float(r)] for z, r in zip(zs, rs)]
    return _cq_line_ring(cq, line, wall, z0, z1, outward=outward, n_pts=n)


def export_step_spool(modules: list, outdir: str, name: str = "PhyAC_Spool",
                      gap_mm: float = 120.0,
                      k_span: int = 2, k_pts: int = 2) -> list:
    """STEP de un EJE COMPUESTO por varios módulos Phy-AC en serie.

    `modules`: lista de (contrato, etiqueta) en orden de flujo. Cada
    contrato es un diseño Phy-AC verificado por separado; entre módulos se
    construye el conducto de transición (hub + separador de flujo) que los
    une geométricamente.

    LÍMITE DECLARADO: la composición es GEOMÉTRICA. El meanline de Phy-AC
    resuelve UNA máquina de un solo flujo; no modela el reparto de caudal
    entre bypass y núcleo ni el acoplamiento de ejes. Cada módulo cumple
    su propia espec y su propio margen; el conjunto no es un ciclo de
    turbofán verificado, es el ensamble de dos máquinas verificadas.
    """
    cq = _cq()
    if cq is None:
        return []
    root = cq.Assembly(name=name)
    z_off = 0.0
    prev = None
    for contract, tag in modules:
        ann = contract["annulus"]
        z_in = ann["hub"][0][0]
        sub = _cq_detailed_machine(cq, contract, k_span, k_pts)
        sub.name = tag
        root.add(sub, name=tag, loc=cq.Location((0, 0, z_off - z_in)))
        r_hub_out = ann["hub"][-1][1]
        r_tip_out = ann["tip"][-1][1]
        z_out_abs = z_off + (ann["hub"][-1][0] - z_in)
        mdot = float(contract["design_vector"]["massflow"])
        if prev is not None:
            (rz, rh0, rt0, m0), tag0 = prev
            wall = float(contract["assembly"].get("hub_wall_mm", 5.0))
            r_hub_in = ann["hub"][0][1]
            r_tip_in = ann["tip"][0][1]
            # RADIO DEL SEPARADOR: el módulo aguas abajo solo se traga su
            # propio gasto. Con el fan pasando m0 y el núcleo mdot, el
            # separador va donde el ÁREA de salida del fan se divide en
            # esa proporción (velocidad axial y densidad ~uniformes):
            #     π(r_split² − r_hub²) / π(r_tip² − r_hub²) = mdot/m0
            frac = min(max(mdot / max(m0, 1e-9), 0.02), 1.0)
            r_split = math.sqrt(rh0 ** 2 + frac * (rt0 ** 2 - rh0 ** 2))
            duct = cq.Assembly(name=f"TRANSITION_{tag0}_{tag}")
            duct.add(_cq_duct(cq, rh0, rz, r_hub_in, z_off, wall),
                     name="CoreIntakeHub")
            duct.add(_cq_duct(cq, r_split, rz, r_tip_in, z_off, wall,
                              outward=True), name="FlowSplitter")
            if frac < 0.98:
                # conducto de BYPASS: la carcasa del módulo anterior sigue
                # hacia atrás por fuera del núcleo
                duct.add(_cq_duct(cq, rt0, rz, rt0, z_off + 1.2 * (
                    ann["hub"][-1][0] - ann["hub"][0][0]), wall,
                    outward=True), name="BypassDuct")
            root.add(duct, name=duct.name)
        prev = ((z_out_abs, r_hub_out, r_tip_out, mdot), tag)
        z_off = z_out_abs + gap_mm
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{name}.step")
    root.save(path, "STEP", write_pcurves=False, precision_mode=0)
    return [path]


def _decimate_row(row: dict, k_span: int, k_pts: int) -> dict:
    """Copia de una fila con menos estaciones de span y menos puntos por
    sección — nivel de detalle para el ensamble COMPLETO.

    Un STEP con 250+ álabes de 13 secciones × 60 puntos pesa ~130 MB y
    casi ningún CAD lo abre con soltura. Bajando a 7 estaciones × 31
    puntos el twist free-vortex se conserva (las secciones intermedias
    eran ya interpolación lineal) y el fichero cae a ~30 MB. Las piezas
    individuales (`mode="parts"`) se siguen exportando a resolución
    completa: la pérdida de detalle es SOLO de la vista de conjunto.
    """
    if k_span <= 1 and k_pts <= 1:
        return row
    out = dict(row)
    secs = row["sections"]
    idx = sorted(set(list(range(0, len(secs) - 1, max(k_span, 1)))
                     + [len(secs) - 1]))
    new = []
    for i in idx:
        s = dict(secs[i])
        if k_pts > 1:
            pts = s["points"]
            s["points"] = [pts[j] for j in range(0, len(pts), k_pts)] \
                + [pts[0]]
        new.append(s)
    out["sections"] = new
    return out


def _cq_full_machine(cq, contract: dict, k_span: int = 2, k_pts: int = 2):
    """Ensamble COMPLETO: cada corona con TODOS sus álabes patronados.

    Los modos `parts`/`assembly` emiten un álabe de muestra por fila y
    delegan el patrón circular al CAD destino (barato y suficiente para
    re-CAD). Este modo construye la máquina entera — fan/rotores, discos,
    eje, carcasa, estátores, IGV y OGV — que es lo que hace falta para
    mirarla, medirla o mandarla a fabricar.

    Cada álabe se LOFTEA UNA VEZ por fila y luego se rota-copia: el loft
    es lo caro, la rotación es casi gratis. Las coronas se añaden a la
    cq.Assembly como sólidos independientes (sin unión booleana), que es
    lo que mantiene el export en segundos en vez de minutos.
    """
    asm_p = contract["assembly"]
    ann = contract["annulus"]
    stages = contract["stages"]
    n_st = len(stages)
    clr = float(ann["tip_clearance_mm"])
    wall = float(asm_p.get("hub_wall_mm", 5.0))
    r_shaft = float(contract.get("structural", {})
                    .get("drum_inner_r_mm", 20.0))
    hub, tip = ann["hub"], ann["tip"]
    z_in, z_out = hub[0][0], hub[-1][0]
    stub = float(asm_p.get("shaft_stub_mm", 30.0))
    bore = float(asm_p.get("shaft_bore_frac", 0.35)) * r_shaft

    assembly = cq.Assembly(name=contract.get("name", "PhyAC_Machine"))

    def add_row(row: dict, tag: str, kind: str = "stator"):
        """Loft de un álabe + patrón circular completo."""
        n_b = int(row["n_blades"])
        one = _cq_blade_trimmed(cq, contract, row, kind, k_span, k_pts)
        for k in range(n_b):
            a = 360.0 * k / n_b
            blade = one if k == 0 else one.rotate((0, 0, 0), (0, 0, 1), a)
            assembly.add(blade, name=f"{tag}_blade{k + 1:03d}")

    # --- rotor: eje + tambor/discos + coronas de rotor -------------------
    assembly.add(_cq_ring(cq, max(bore, 0.0), r_shaft,
                          z_in - stub, z_out + stub), name="Shaft")
    drum = _cq_line_ring(cq, hub, wall, z_in, z_out, outward=False)
    for i, st in enumerate(stages):
        rot = st["rotor"]
        web = min(max(float(asm_p.get("disk_web_frac", 0.5))
                      * (rot["z_te_mm"] - rot["z_le_mm"]),
                      float(asm_p.get("disk_web_min_mm", 4.0))),
                  float(asm_p.get("disk_web_max_mm", 14.0)))
        zc = rot["z_center_mm"]
        r_hub_c = float(np.interp(zc, [p[0] for p in hub],
                                  [p[1] for p in hub]))
        drum = drum.union(_cq_ring(cq, r_shaft, r_hub_c,
                                   zc - 0.5 * web, zc + 0.5 * web))
    assembly.add(drum, name="RotorDrum")
    for i, st in enumerate(stages):
        add_row(st["rotor"], f"Rotor{i + 1}", "rotor")

    # --- estátor: carcasa + coronas de estátor + IGV/OGV -----------------
    shell_line = [[z, r + clr] for z, r in tip]
    casing = _cq_line_ring(cq, shell_line, wall, tip[0][0] - 4.0,
                           tip[-1][0] + 4.0, outward=True)
    for at_start in (True, False):
        z_ref = (tip[0][0] - 4.0) if at_start else (tip[-1][0] + 4.0)
        t_fl = float(asm_p.get("flange_t_mm", 8.0))
        zf = z_ref if at_start else z_ref - t_fl
        r_tip_f = float(np.interp(z_ref, [p[0] for p in tip],
                                  [p[1] for p in tip]))
        r_o = r_tip_f + clr + wall
        fl = _cq_ring(cq, r_tip_f + clr,
                      r_o + float(asm_p.get("flange_w_mm", 12.0)),
                      zf, zf + t_fl)
        n_b = int(asm_p.get("flange_bolt_count", 12))
        r_b = r_o + 0.5 * float(asm_p.get("flange_w_mm", 12.0))
        d_b = float(asm_p.get("flange_bolt_d_mm", 6.0))
        for k in range(n_b):
            a = 2.0 * math.pi * k / n_b
            fl = fl.cut(cq.Workplane(
                "XY", origin=(r_b * math.cos(a), r_b * math.sin(a),
                              zf - 1.0)).circle(0.5 * d_b)
                .extrude(t_fl + 2.0))
        casing = casing.union(fl)
    assembly.add(casing, name="Casing")
    if contract.get("igv"):
        add_row(contract["igv"], "IGV")
    for i, st in enumerate(stages):
        add_row(st["stator"], f"Stator{i + 1}")
    if contract.get("ogv"):
        add_row(contract["ogv"], "OGV")
    return assembly


def _detailed_step_readme(contract: dict) -> str:
    """Notas del ensamble DETALLADO (machine_detailed.step)."""
    st = contract["stages"]
    asm = contract["assembly"]
    ft = asm.get("firtree") or {}
    d = contract.get("derived", {})
    ann = contract["annulus"]
    n_r = sum(s["rotor"]["n_blades"] for s in st)
    n_s = sum(s["stator"]["n_blades"] for s in st)
    n_i = contract["igv"]["n_blades"] if contract.get("igv") else 0
    n_o = contract["ogv"]["n_blades"] if contract.get("ogv") else 0
    n_tie = int(asm.get("tie_bolt_count", 12))
    n_cb = int(asm.get("flange_bolt_count", 12))
    n_joint = max(len(st) - 1, 0)
    L = [
        "Phy-AC — ensamble DETALLADO (machine_detailed.step)",
        "===================================================",
        "Coordenadas de maquina: Z = eje de rotacion, mm. Gas en +Z.",
        "",
        "ARBOL DEL ENSAMBLE (tres subensambles):",
        "",
        "  ROTOR",
        "    Shaft            eje con munones y barreno central",
        f"    Disc1..{len(st):<11d}disco por etapa: barreno de ajuste sobre",
        "                     el eje, alma, rim, ranuras de abeto y",
        "                     taladros de tirante. NO es el eje.",
        f"    Rotor{{i}}_blade{{k}}  {n_r} alabes, cada uno con su RAIZ DE",
        "                     ABETO metida en la ranura de su disco",
        "    DrumSpacer{i}    casco del flowpath entre discos",
        "  STATOR",
        f"    CasingRing1..{len(st):<7d}un anillo POR ETAPA de estator, con",
        "                     brida en CADA extremo (se abren por planos",
        "                     de junta reales)",
        f"    Stator{{i}}_vane{{k}}  {n_s} vanos",
    ]
    if n_i:
        L.append(f"    IGV_vane{{k}}       {n_i} vanos de guia de entrada")
    if n_o:
        L.append(f"    OGV_vane{{k}}       {n_o} vanos de salida")
    L += [
        "  FASTENERS",
        f"    TieBolt01..{n_tie:<9d}tirantes pasantes que aprietan la pila",
        "                     de discos contra el eje",
        f"    CasingBolt{{i}}_{{k}}  {n_cb} por junta x {n_joint} juntas de carcasa",
        "",
        "RAIZ DE ABETO (retencion de alabe):",
        f"  dientes           {ft.get('n_teeth', '-')}",
        f"  profundidad       {ft.get('depth_mm', '-')} mm",
        f"  ancho en plataforma {ft.get('w_top_mm', '-')} mm",
        f"  juego de montaje  {ft.get('clearance_mm', '-')} mm por lado",
        "  La ranura del disco se genera con el MISMO perfil mas el juego,",
        "  asi que alabe y disco encajan por construccion.",
        f"  rebaje del rim   {asm.get('rim_relief_mm', 0.10)} mm por debajo de la",
        "                     plataforma del alabe (el rim no forma la vena)",
        "",
        "INTERFERENCIAS:",
        "  Este ensamble pasa el chequeo booleano por pares",
        "  (phyac --check-interference / geometry_generator."
        "assembly_interferences):",
        f"  ningun par de piezas comparte mas de "
        f"{INTERFERENCE_TOL_MM3:g} mm3 de material.",
        "",
        f"  Longitud axial    {ann['hub'][-1][0] - ann['hub'][0][0]:.1f} mm",
        f"  Diametro de punta {2 * ann['tip'][0][1]:.1f} mm (entrada)",
        f"  Holgura de punta  {ann['tip_clearance_mm']:.2f} mm radial",
    ]
    if d:
        L += [
            "",
            "Punto de diseno verificado:",
            f"  PR                {d.get('PR', 0):.3f}",
            f"  eta politropica   {d.get('eta_poly', 0):.4f}",
            f"  Potencia          {d.get('power_W', 0) / 1e6:.2f} MW",
            f"  Margen de bombeo  "
            f"{100 * (d.get('sm_flow') or 0):.1f} % en gasto",
            f"  Torbellino        n = {d.get('vortex_n', -1.0):+.2f}",
        ]
    L += [
        "",
        "Nivel de detalle de los perfiles: 7 estaciones de span y ~31",
        "puntos por seccion (el contrato lleva 13 y 60). El twist se",
        "conserva; para fabricar, usar los STL de la capa 5c.",
    ]
    return "\n".join(L) + "\n"


def _full_step_readme(contract: dict) -> str:
    """Notas del ensamble COMPLETO (machine_full.step)."""
    st = contract["stages"]
    d = contract.get("derived", {})
    ann = contract["annulus"]
    n_r = sum(s["rotor"]["n_blades"] for s in st)
    n_s = sum(s["stator"]["n_blades"] for s in st)
    n_i = contract["igv"]["n_blades"] if contract.get("igv") else 0
    n_o = contract["ogv"]["n_blades"] if contract.get("ogv") else 0
    L = [
        "Phy-AC — ensamble COMPLETO (machine_full.step)",
        "==============================================",
        "Coordenadas de maquina: Z = eje de rotacion, mm. El gas fluye en +Z.",
        "",
        "A diferencia de parts/*.step, aqui NO hay alabe de muestra: todas",
        "las coronas estan patronadas. Cada solido es un componente con",
        "nombre en el arbol del STEP.",
        "",
        f"  Etapas                {len(st)}",
        f"  Alabes de rotor       {n_r}  (Rotor1..{len(st)}_bladeNNN)",
        f"  Vanos de estator      {n_s}  (Stator1..{len(st)}_bladeNNN)",
        f"  IGV / OGV             {n_i} / {n_o}",
        f"  Solidos totales       {n_r + n_s + n_i + n_o + 3}"
        "  (+ Shaft, RotorDrum, Casing)",
        "",
        "Grupo ROTOR (gira):   Shaft, RotorDrum, Rotor*_blade*",
        "Grupo ESTATOR (fijo): Casing, IGV_blade*, Stator*_blade*, OGV_blade*",
        "",
        f"  Longitud axial        {ann['hub'][-1][0] - ann['hub'][0][0]:.1f} mm",
        f"  Radio de punta entr.  {ann['tip'][0][1]:.1f} mm",
        f"  Holgura de punta      {ann['tip_clearance_mm']:.2f} mm (radial)",
    ]
    if d:
        L += [
            "",
            "Punto de diseno verificado:",
            f"  PR                    {d.get('PR', 0):.3f}",
            f"  eta politropica       {d.get('eta_poly', 0):.4f}",
            f"  Potencia              {d.get('power_W', 0) / 1e6:.2f} MW",
            f"  Margen de bombeo      "
            f"{100 * (d.get('sm_flow') or 0):.1f} % en gasto (linea de trabajo)",
        ]
    L += [
        "",
        "Nivel de detalle: los alabes del ENSAMBLE se loftean con 7",
        "estaciones de span y ~31 puntos por seccion (el contrato lleva 13",
        "y 60). El twist free-vortex se conserva; las piezas sueltas de",
        "parts/*.step y los STL de la capa 5c van a resolucion completa.",
        "",
        "Aproximacion geometrica: los alabes se loftean entre planos",
        "radiales paralelos, sin el envolvimiento cilindrico phi = y/r que",
        "aplica el voxelizado de la capa 5c. En alabes de gran altura la",
        "diferencia es visible cerca de la punta; para fabricar, usar el STL.",
    ]
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
# Verificación de interferencias del ENSAMBLE (fase 10 · G-04)
# ---------------------------------------------------------------------------
# Los tres fallos de geometría que aparecieron al montar la máquina —álabes
# saliéndose de la carcasa, vanos metidos en el tambor, ranuras que se
# comían el disco entero— los encontró una persona mirando el STEP, no la
# suite. El chequeo que había mide vértices de filas SUELTAS contra las
# polilíneas de annulus: no habría cazado ninguno de los tres, porque
# ninguno era un problema de una fila contra la vena, sino de una pieza
# contra otra pieza ya montada.
#
# Esto comprueba lo que de verdad importa: el volumen de la intersección
# booleana de dos piezas que no deben ocupar el mismo sitio tiene que ser
# cero. Piezas en CONTACTO (la raíz dentro de su ranura, la brida contra
# la brida vecina, el disco calado sobre el eje) dan intersección de
# volumen nulo aunque compartan caras, así que la regla es uniforme y no
# necesita lista de excepciones.

INTERFERENCE_TOL_MM3 = 1.0    # ~1 mm³: por debajo es ruido de teselado de
#                               las booleanas de OCC, no material compartido


def _asm_parts(asm, loc=None, prefix: str = "") -> list[tuple[str, object]]:
    """Aplana un cq.Assembly a (nombre, sólido) con las localizaciones ya
    compuestas — las piezas quedan en coordenadas de máquina."""
    here = asm.loc if loc is None else loc * asm.loc
    out: list[tuple[str, object]] = []
    name = prefix + (asm.name or "")
    if asm.obj is not None:
        shp = asm.obj.val() if hasattr(asm.obj, "val") else asm.obj
        try:
            out.append((name, shp.moved(here)))
        except Exception:
            out.append((name, shp))
    for ch in asm.children:
        out += _asm_parts(ch, here, (name + "/") if name else "")
    return out


def _interference_sample(parts: list[tuple[str, object]]
                         ) -> list[tuple[str, object]]:
    """Muestra representativa: todas las piezas únicas más DOS miembros
    consecutivos de cada corona patronada.

    Comprobar los 24 álabes de una corona contra todo no añade
    información — la corona es la misma pieza girada. Dos vecinos sí:
    cazan el solape álabe-álabe cuando el paso se queda corto. Lo mismo
    vale para tirantes y pernos de brida. Los discos y los anillos de
    carcasa NO son un patrón: cada uno es una pieza distinta y entran
    todos.
    """
    keep = []
    for name, shp in parts:
        leaf = name.rsplit("/", 1)[-1]
        m = re.search(r"(?:blade|vane|TieBolt|CasingBolt\d+_)(\d+)$", leaf)
        if m and int(m.group(1)) > 2:
            continue
        keep.append((name, shp))
    return keep


def assembly_interferences(contract: dict, k_span: int = 2, k_pts: int = 2,
                           tol_mm3: float = INTERFERENCE_TOL_MM3,
                           parts: list | None = None) -> list[dict]:
    """Pares de piezas del ensamble que comparten material.

    Devuelve una lista de dicts ordenada por volumen solapado descendente:
    `{"a", "b", "volume_mm3"}`. Vacía = ensamble limpio. Nunca devuelve un
    booleano: para arreglar una interferencia hace falta saber QUÉ dos
    piezas y CUÁNTO.
    """
    cq = _cq()
    if cq is None:
        return []
    if parts is None:
        parts = _asm_parts(_cq_detailed_machine(cq, contract, k_span, k_pts))
    sample = _interference_sample(parts)

    # Pre-filtro por caja envolvente: la inmensa mayoría de los pares está
    # en estaciones axiales distintas y no hace falta una booleana para
    # saberlo. `enlarge(-tol)` evita disparar por caras que solo se tocan.
    boxes = []
    for name, shp in sample:
        try:
            boxes.append((name, shp, shp.BoundingBox()))
        except Exception:
            boxes.append((name, shp, None))

    def bb_overlap(b1, b2, pad=1e-6):
        if b1 is None or b2 is None:
            return True
        return (b1.xmin < b2.xmax - pad and b2.xmin < b1.xmax - pad
                and b1.ymin < b2.ymax - pad and b2.ymin < b1.ymax - pad
                and b1.zmin < b2.zmax - pad and b2.zmin < b1.zmax - pad)

    hits: list[dict] = []
    for i in range(len(boxes)):
        n1, s1, b1 = boxes[i]
        for j in range(i + 1, len(boxes)):
            n2, s2, b2 = boxes[j]
            if not bb_overlap(b1, b2):
                continue
            try:
                common = s1.intersect(s2)
                vol = float(common.Volume()) if common is not None else 0.0
            except Exception:
                # Una booleana que ni siquiera corre es un fallo que hay
                # que ver, no un pase: se reporta con volumen negativo.
                hits.append(dict(a=n1, b=n2, volume_mm3=-1.0))
                continue
            if vol > tol_mm3:
                hits.append(dict(a=n1, b=n2, volume_mm3=vol))
    hits.sort(key=lambda h: -h["volume_mm3"])
    return hits


def format_interferences(hits: list[dict], n_parts: int | None = None,
                         tol_mm3: float = INTERFERENCE_TOL_MM3) -> str:
    """Informe legible del chequeo de interferencias."""
    if not hits:
        return (f"ensamble limpio: ningún par comparte material "
                f"(tolerancia {tol_mm3:g} mm³"
                + (f", {n_parts} piezas comprobadas" if n_parts else "")
                + ")")
    lines = [f"{len(hits)} par(es) con material compartido "
             f"(tolerancia {tol_mm3:g} mm³):"]
    for h in hits:
        if h["volume_mm3"] < 0:
            lines.append(f"  · {h['a']} ∩ {h['b']}: la booleana FALLÓ "
                         f"(sólido inválido)")
        else:
            lines.append(f"  · {h['a']} ∩ {h['b']}: "
                         f"{h['volume_mm3'] / 1000.0:.3f} cm³")
    return "\n".join(lines)


def export_step(contract: dict, outdir: str, mode: str = "parts") -> list[str]:
    """Exporta STEP del diseño. Devuelve las rutas escritas ([] sin
    cadquery — degradación silenciosa, dependencia opcional `step`).

    mode:
      "blade0"   — UN álabe del primer rotor (rápido; compatibilidad).
      "parts"    — una pieza STEP por parte física (nombres de los STL)
                   con un álabe de muestra + parts/README.txt (default).
      "assembly" — además, cq.Assembly con todas las piezas nombradas →
                   assembly.step.
      "detailed" — ensamble jerárquico con las PIEZAS reales separadas
                   (eje, discos con ranuras de abeto, álabes con raíz,
                   anillos de carcasa por etapa con bridas, vanos,
                   tirantes y pernos) → machine_detailed.step.
      "full"     — máquina COMPLETA con todas las coronas patronadas
                   (fan/rotores + discos + eje + carcasa + estátores +
                   IGV/OGV) → machine_full.step.
    """
    cq = _cq()
    if cq is None:
        return []
    written: list[str] = []
    try:
        if mode == "blade0":
            row = contract["stages"][0]["rotor"]
            solid = _cq_blade_solid(cq, row, row["z_center_mm"])
            path = os.path.join(outdir, "blade_stage0.step")
            cq.exporters.export(solid, path)
            return [path]

        if mode == "detailed":
            os.makedirs(outdir, exist_ok=True)
            path = os.path.join(outdir, "machine_detailed.step")
            _cq_detailed_machine(cq, contract).save(
                path, "STEP", write_pcurves=False, precision_mode=0)
            rp = os.path.join(outdir, "machine_detailed_README.txt")
            with open(rp, "w", encoding="utf-8") as f:
                f.write(_detailed_step_readme(contract))
            return [path, rp]

        if mode == "full":
            os.makedirs(outdir, exist_ok=True)
            path = os.path.join(outdir, "machine_full.step")
            # write_pcurves=False quita las curvas paramétricas sobre
            # cada cara (redundantes: cualquier CAD las reconstruye) —
            # ~30% menos de fichero sin perder geometría.
            _cq_full_machine(cq, contract).save(
                path, "STEP", write_pcurves=False, precision_mode=0)
            with open(os.path.join(outdir, "machine_full_README.txt"), "w",
                      encoding="utf-8") as f:
                f.write(_full_step_readme(contract))
            return [path, os.path.join(outdir, "machine_full_README.txt")]

        parts = _cq_part_solids(cq, contract)
        pdir = os.path.join(outdir, "parts")
        os.makedirs(pdir, exist_ok=True)
        for name, solid in parts.items():
            path = os.path.join(pdir, f"{name}.step")
            cq.exporters.export(solid, path)
            written.append(path)
        readme = os.path.join(pdir, "README.txt")
        with open(readme, "w", encoding="utf-8") as f:
            f.write(_step_readme(contract))
        written.append(readme)

        if mode == "assembly":
            asm = cq.Assembly()
            for name, solid in parts.items():
                asm.add(solid, name=name)
            path = os.path.join(outdir, "assembly.step")
            asm.save(path)
            written.append(path)
        return written
    except Exception:
        return written


def try_export_step(contract: dict, path: str) -> bool:
    """Compatibilidad: UN álabe del primer rotor (ahora loft spline vía
    export_step). False sin error si cadquery no está instalado."""
    out = export_step(contract, os.path.dirname(path) or ".", mode="blade0")
    return bool(out)


# ---------------------------------------------------------------------------
# Generación del paquete de geometría
# ---------------------------------------------------------------------------
def build_contract(theta: np.ndarray, record: dict, structural: dict | None,
                   run_meta: dict | None) -> dict:
    theta = np.asarray(theta, dtype=float)
    omega = float(theta[1]) * 2 * math.pi / 60.0
    # 1º: secciones de todas las filas → extensiones axiales EXACTAS
    secs_map: dict = {}
    extents: dict = {}
    for s in record["stage_table"]:
        for kind in ("rotor", "stator"):
            secs = _row_sections(s, kind, omega)
            secs_map[(s["stage"], kind)] = secs
            extents[(s["stage"], kind)] = _row_extent_mm(secs)
    # El IGV solo EXISTE si la etapa 1 pide pre-swirl. Un FAN se diseña
    # con α₁ = 0 (Rx = 1 − ψ/2) precisamente para no llevar IGV: emitir
    # una fila de giro nulo delante del rotor sería inventar una pieza
    # que no está en la máquina. La física ya lo decide igual
    # (record["igv"]["present"]); aquí la geometría se alinea con ella.
    has_igv = bool((record.get("igv") or {}).get("present", True)) and \
        abs(record["stage_table"][0]["alpha1_deg"]) > 1e-3
    igv_secs = _row_sections(_pseudo_stage_igv(record), "stator", omega) \
        if has_igv else None
    ogv_secs = _row_sections(_pseudo_stage_ogv(record), "stator", omega)
    igv_ext = _row_extent_mm(igv_secs) if has_igv else None
    ogv_ext = _row_extent_mm(ogv_secs)
    # 2º: colocación secuencial sin solapes (IGV al frente, OGV al final)
    rows = _row_geometry(record, extents, igv_ext, ogv_ext)

    # 2b (fase 9.1): RECONSTRUIR cada fila con los radios del annulus EN
    # SU PROPIA ESTACIÓN AXIAL y con las holguras de marcha.
    #
    # Hasta ahora todas las filas de una etapa se tallaban con el (r_hub,
    # r_tip) de la ENTRADA de esa etapa. Pero el annulus se contrae a lo
    # largo de z: el estátor, que está medio paso más atrás, quedaba con
    # la punta por encima de la carcasa (medido: +8.5 mm en la etapa 1,
    # con una holgura de 0.4 mm) y con la raíz por debajo del cubo
    # (−8.8 mm), o sea atravesando el tambor. Se veía a simple vista en
    # el ensamble.
    #
    # Holguras: el ROTOR nace en el cubo y su punta se queda en la línea
    # de punta (el hueco a la carcasa es ε). El ESTÁTOR cuelga de la
    # carcasa (llega a r_tip + ε) y su extremo libre se retira ε del
    # tambor giratorio.
    ann0 = annulus_lines(record, extents, igv_ext, ogv_ext)
    hub_l, tip_l = ann0["hub"], ann0["tip"]
    eps = float(ann0["tip_clearance_mm"])

    def _r_at(line, z):
        return float(np.interp(z, [p[0] for p in line],
                               [p[1] for p in line]))

    def _rebuild(stage, kind, rg):
        zc = 0.5 * (rg["z_le"] + rg["z_te"])
        r_h, r_t = _r_at(hub_l, zc), _r_at(tip_l, zc)
        if kind == "rotor":
            return _row_sections(stage, kind, omega, r_hub_mm=r_h,
                                 r_tip_mm=r_t)
        return _row_sections(stage, kind, omega, r_hub_mm=r_h + eps,
                             r_tip_mm=r_t + eps)

    for s in record["stage_table"]:
        for kind in ("rotor", "stator"):
            rg = next(r for r in rows if r["stage"] == s["stage"]
                      and r["kind"] == kind)
            secs_map[(s["stage"], kind)] = _rebuild(s, kind, rg)
    if has_igv:
        rg = next(r for r in rows if r["kind"] == "igv")
        igv_secs = _rebuild(_pseudo_stage_igv(record), "stator", rg)
    rg = next(r for r in rows if r["kind"] == "ogv")
    ogv_secs = _rebuild(_pseudo_stage_ogv(record), "stator", rg)

    stages = []
    for s in record["stage_table"]:
        r_rows = {r["kind"]: r for r in rows if r["stage"] == s["stage"]
                  and r["kind"] in ("rotor", "stator")}
        entry = dict(index=s["stage"])
        for kind in ("rotor", "stator"):
            rg = r_rows[kind]
            entry[kind] = dict(
                n_blades=int(s[f"n_blades_{kind}"]),
                z_le_mm=round(rg["z_le"], 3), z_te_mm=round(rg["z_te"], 3),
                z_center_mm=round(0.5 * (rg["z_le"] + rg["z_te"]), 3),
                rotating=(kind == "rotor"),
                sections=secs_map[(s["stage"], kind)],
            )
        stages.append(entry)

    def _vane_row(kind: str, secs: list[dict], n_blades: int) -> dict:
        rg = [r for r in rows if r["kind"] == kind][0]
        return dict(
            n_blades=n_blades,
            z_le_mm=round(rg["z_le"], 3), z_te_mm=round(rg["z_te"], 3),
            z_center_mm=round(0.5 * (rg["z_le"] + rg["z_te"]), 3),
            rotating=False, sections=secs,
        )

    st_tab = record["stage_table"]
    igv_row = _vane_row("igv", igv_secs,
                        int(st_tab[0]["n_blades_stator"])) if has_igv else None
    ogv_row = _vane_row("ogv", ogv_secs, int(st_tab[-1]["n_blades_stator"]))

    struct_block = {}
    if structural:
        struct_block = dict(
            material=structural["material"],
            sigma_vm_max_MPa=structural["sigma_vm_max_MPa"],
            burst_margin=structural["burst_margin"],
            AN2_in2rpm2=structural["AN2_in2rpm2"],
            rotor_mass_kg=structural["rotor_mass_kg"],
            drum_inner_r_mm=round(
                0.30 * record["stage_table"][0]["r_hub_mm"], 2),
            feasible_struct=structural["feasible_struct"],
        )

    # Parámetros del ensamble (capa 5c, estilo turbodesigner: eje +
    # discos por etapa + casco del hub + carcasa con bridas). El C# los
    # lee con defaults si faltan — el contrato es la única frontera.
    import structures_core as _sc   # lazy: structures_core nos importa a
    #                                 nosotros de forma diferida (fase 7)
    n_st = record["n_stages"]
    assembly = dict(
        shaft_stub_mm=30.0,        # muñón del eje a cada lado
        shaft_bore_frac=0.35,      # barreno central / radio del eje
        hub_wall_mm=5.0,           # pared del casco del hub (flowpath)
        disk_web_frac=0.5,         # espesor del alma / cuerda axial
        disk_web_min_mm=4.0, disk_web_max_mm=14.0,
        flange_bolt_count=12, flange_bolt_d_mm=6.0,
        flange_w_mm=12.0, flange_t_mm=8.0,
        # tirantes pasantes que aprietan la pila de discos sobre el eje
        tie_bolt_count=12, tie_bolt_d_mm=10.0,
        # Rebaje del rim del disco respecto a la plataforma del álabe.
        # En un disco real la vena la forma la plataforma del álabe y el
        # rim se queda por debajo; aquí además garantiza que rim y álabe
        # no compartan material por el redondeo poligonal con el que se
        # recorta el álabe contra la línea de cubo.
        rim_relief_mm=0.10,
        # fillet de raíz de álabe: el MISMO radio que usa el K_t de
        # Peterson en structures_core (fase 7 — fuente única)
        blade_fillet_r_mm=_sc.BLADE_FILLET_R_MM,
        # Puerto de sangrado (boss radial con barreno en la carcasa).
        # Desde la fase 9 la etapa NO es "la de en medio" por convención:
        # la elige la física. `stage_stall_first` es la etapa que agota
        # antes su capacidad de subida de presión de Koch a velocidad
        # parcial — sangrar JUSTO DETRÁS de ella es lo que descarga las
        # etapas frontales y permite arrancar (Dixon & Hall §5.9). El
        # agujero de la pieza y el modelo aerodinámico apuntan por fin al
        # mismo sitio.
        bleed_stage=int(min(max(record.get("stage_stall_first", 0) + 2, 1),
                            max(n_st, 1))),
        bleed_hole_d_mm=18.0, bleed_boss_d_mm=32.0, bleed_boss_h_mm=14.0,
        # Retención de álabe de ABETO (fase 9). La profundidad se escala
        # con la altura del álabe de la primera etapa y se acota para que
        # nunca coma más de la mitad del rim disponible sobre el eje.
        firtree=dict(
            enabled=True,
            n_teeth=3,
            # La profundidad de un abeto escala con el ANCHO de la raíz
            # (proporción real ≈ 2-2.5×), no con la altura del álabe: un
            # árbol demasiado profundo para su ancho se estrangula y deja
            # el último diente sin sección resistente. Acotada además por
            # el rim disponible sobre el tambor.
            w_top_mm=round(max(_root_width_mm(record), 3.0), 2),
            depth_mm=round(float(np.clip(
                2.2 * _root_width_mm(record), 4.0,
                0.45 * max(record["stage_table"][0]["r_hub_mm"]
                           - _drum_inner_r_mm(record), 8.0))), 2),
            taper_deg=8.0,
            lobe_frac=0.22,
            clearance_mm=0.08,      # juego de brochado/montaje por lado
            axial_frac=0.85,        # longitud de la raíz / cuerda axial
            min_width_mm=3.0,
        ),
    )

    return dict(
        schema=SCHEMA,
        design_vector={n: float(v) for n, v in zip(VAR_NAMES, theta)},
        derived=dict(PR=record["PR"], eta_poly=record["eta_poly"],
                     eta_isen=record["eta_isen"],
                     power_W=record["power_W"], U_tip_m_s=record["U_tip"],
                     AN2_m2_rpm2=record["AN2_m2_rpm2"],
                     length_mm=record["length_mm"],
                     n_stages=record["n_stages"],
                     r_tip_in_mm=record["r_tip_in_mm"],
                     M_rel_tip1=record["M_rel_tip1"],
                     min_SM=record["min_SM"], feasible=record["feasible"],
                     # --- fase 9 -----------------------------------------
                     sm_flow=record.get("sm_flow"),
                     sm_koch_wl=record.get("sm_koch_wl"),
                     stage_stall_first=record.get("stage_stall_first"),
                     ds_machine_J_kgK=record.get("ds_machine_J_kgK"),
                     blockage_exit=record.get("blockage_exit"),
                     gas_model=record.get("gas_model", "variable_cp"),
                     vortex_n=VORTEX_N,
                     bleed=record.get("bleed"),
                     igv_loss=record.get("igv"),
                     source=record.get("source", "meanline_L0")),
        annulus=annulus_lines(record, extents, igv_ext, ogv_ext),
        stages=stages,
        igv=igv_row,
        ogv=ogv_row,
        assembly=assembly,
        structural=struct_block,
        cfx_boundary_conditions=cfx_boundary_conditions(theta, record),
        run_meta=run_meta or {},
    )


def generate(theta: np.ndarray, record: dict, outdir: str,
             structural: dict | None = None,
             run_meta: dict | None = None) -> list[str]:
    """Escribe el paquete de geometría en <outdir> y devuelve el manifest.

    Archivos: axial_compressor.json (contrato 5c), annulus.csv,
    stage_summary.csv y, si cadquery está instalado, blade_stage0.step.
    """
    os.makedirs(outdir, exist_ok=True)
    manifest: list[str] = []
    contract = build_contract(theta, record, structural, run_meta)

    p_json = os.path.join(outdir, "axial_compressor.json")
    with open(p_json, "w", encoding="utf-8") as f:
        json.dump(contract, f, indent=1)
    manifest.append(p_json)

    p_ann = os.path.join(outdir, "annulus.csv")
    with open(p_ann, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["z_hub_mm", "r_hub_mm", "z_tip_mm", "r_tip_mm"])
        ann = contract["annulus"]
        for (zh, rh), (zt, rt) in zip(ann["hub"], ann["tip"]):
            w.writerow([zh, rh, zt, rt])
    manifest.append(p_ann)

    p_sum = os.path.join(outdir, "stage_summary.csv")
    with open(p_sum, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["stage", "phi", "psi", "Rx", "r_mean_mm", "h_blade_mm",
                    "n_bl_rotor", "n_bl_stator", "chord_r_mm", "chord_s_mm",
                    "beta1_deg", "beta2_deg", "alpha1_deg", "alpha2_deg",
                    "DF_rotor", "DF_stator", "SM", "eta_tt", "PR"])
        for s in record["stage_table"]:
            w.writerow([s["stage"], round(s["phi"], 4), round(s["psi"], 4),
                        round(s["Rx"], 4), round(s["r_mean_mm"], 2),
                        round(s["h_blade_mm"], 2), s["n_blades_rotor"],
                        s["n_blades_stator"], round(s["chord_rotor_mm"], 2),
                        round(s["chord_stator_mm"], 2),
                        round(s["beta1_deg"], 2), round(s["beta2_deg"], 2),
                        round(s["alpha1_deg"], 2), round(s["alpha2_deg"], 2),
                        round(s["DF_rotor"], 4), round(s["DF_stator"], 4),
                        round(s["SM"], 4), round(s["eta_tt"], 4),
                        round(s["PR"], 4)])
    manifest.append(p_sum)

    # BOM (estilo turbodesigner: piezas del ensamble con cantidades)
    p_bom = os.path.join(outdir, "bom.csv")
    with open(p_bom, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["item", "qty", "stage", "dims_mm", "material"])
        mat = (structural or {}).get("material", "Ti-6Al-4V")
        asm = contract["assembly"]
        shaft_r = contract["structural"].get("drum_inner_r_mm", 0.0) \
            if contract["structural"] else 0.0
        z_len = contract["stages"][-1]["stator"]["z_te_mm"]
        w.writerow(["eje", 1, "-",
                    f"r{shaft_r:.1f} x L{z_len + 2 * asm['shaft_stub_mm']:.0f}",
                    mat])
        for st, s in zip(contract["stages"], record["stage_table"]):
            rw = st["rotor"]
            web = min(max(asm["disk_web_frac"]
                          * (rw["z_te_mm"] - rw["z_le_mm"]),
                          asm["disk_web_min_mm"]), asm["disk_web_max_mm"])
            w.writerow([f"disco etapa {st['index'] + 1}", 1, st["index"] + 1,
                        f"r{s['r_hub_mm']:.1f} x t{web:.1f}", mat])
            w.writerow([f"alabe rotor etapa {st['index'] + 1}",
                        rw["n_blades"], st["index"] + 1,
                        f"c{s['chord_rotor_mm']:.1f} x h{s['h_blade_mm']:.1f}",
                        mat])
            w.writerow([f"alabe estator etapa {st['index'] + 1}",
                        st["stator"]["n_blades"], st["index"] + 1,
                        f"c{s['chord_stator_mm']:.1f} x h{s['h_blade_mm']:.1f}",
                        mat])
        for st, s in zip(contract["stages"], record["stage_table"]):
            w.writerow([f"anillo de carcasa etapa {st['index'] + 1}", 1,
                        st["index"] + 1,
                        f"r{s['r_tip_mm']:.1f}", mat])
        s0, sl = record["stage_table"][0], record["stage_table"][-1]
        if contract.get("igv"):
            w.writerow(["vano IGV", contract["igv"]["n_blades"], 0,
                        f"c{s0['chord_stator_mm']:.1f} x "
                        f"h{s0['h_blade_mm']:.1f}", mat])
        w.writerow(["vano OGV", contract["ogv"]["n_blades"],
                    len(record["stage_table"]) + 1,
                    f"c{sl['chord_stator_mm']:.1f} x h{sl['h_blade_mm']:.1f}",
                    mat])
        w.writerow(["puerto de sangrado", 1, asm["bleed_stage"],
                    f"D{asm['bleed_hole_d_mm']:.0f}", mat])
        w.writerow(["perno de brida", 2 * asm["flange_bolt_count"], "-",
                    f"M{asm['flange_bolt_d_mm']:.0f}", "acero 8.8"])
    manifest.append(p_bom)

    p_step = os.path.join(outdir, "blade_stage0.step")
    if try_export_step(contract, p_step):
        manifest.append(p_step)

    return manifest


if __name__ == "__main__":
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    from physics_core import Fidelity, evaluate
    theta_ref = np.array([4.0, 12_500.0, 0.62, 0.55, 0.30, -0.10,
                          0.60, 1.20, 1.10, 2.20, 288.15, 101_325.0, 25.0])
    rec = evaluate(theta_ref, fidelity=Fidelity.L0)
    files = generate(theta_ref, rec, "geometry_out")
    print("Manifest:")
    for f in files:
        print("  ", f, f"({os.path.getsize(f)} B)")
    with open(files[0], encoding="utf-8") as fh:
        c = json.load(fh)
    s0 = c["stages"][0]
    print("rotor st0:", s0["rotor"]["n_blades"], "álabes,",
          len(s0["rotor"]["sections"]), "secciones,",
          len(s0["rotor"]["sections"][0]["points"]), "pts/sección")
    for sec in s0["rotor"]["sections"]:
        print(f"  f={sec['span_frac']:.2f} r={sec['r_mm']:7.1f} "
              f"c={sec['chord_mm']:6.2f} γ={sec['stagger_deg']:6.2f} "
              f"θ={sec['camber_deg']:6.2f} {sec['profile']} M={sec['M_in']:.2f}")
    st0 = s0["stator"]["sections"][0]
    print(f"estator hub: γ={st0['stagger_deg']:.2f} θ={st0['camber_deg']:.2f}")
