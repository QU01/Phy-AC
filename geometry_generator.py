"""
QUASAR Phy-AC · geometry_generator.py
=====================================
Capa 5a del CEM: del record meanline verificado a la GEOMETRÍA de filas de
álabes — el contrato `axial_compressor.json` (schema phyac-axial-1) que
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

import numpy as np

from physics_core import (CP, GAMMA, RGAS, ROW_GAP_FRACTION,
                          TIP_CLEARANCE_MM, VAR_NAMES)
from blade_profiles import (N_POINTS, _circular_arc_camber, dca_profile,
                            metal_angles, naca65_profile)

N_CAMBER = 41              # puntos de la línea de comba (capa 5c)

N_SPAN = 13                # secciones por fila (twist free-vortex EXACTO en
#                            cada estación — con pocas secciones la
#                            interpolación lineal de la capa 5c deja
#                            crestas visibles en las estaciones originales)
CHORD_TAPER = 0.85         # c_tip/c_hub
TC_ROOT_R, TC_TIP_R = 0.10, 0.05   # t/c del rotor hub→tip
TC_ROOT_S, TC_TIP_S = 0.09, 0.06   # t/c del estátor
M_DCA_THRESHOLD = 0.80     # M de entrada sobre el cual la sección es DCA

SCHEMA = "phyac-axial-1"


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
                  n_span: int = N_SPAN) -> list[dict]:
    """Secciones de una fila con triángulos free-vortex y ángulos metálicos.

    Rotor: ángulos relativos β(r); estátor: absolutos α(r) con signo
    NEGATIVO (convención del contrato).
    """
    Cx = stage["Cx"]
    r_m = stage["r_mean_mm"] / 1000.0
    r_hub = stage["r_hub_mm"] / 1000.0
    r_tip = stage["r_tip_mm"] / 1000.0
    Cu1_m = Cx * math.tan(math.radians(stage["alpha1_deg"]))
    Cu2_m = Cx * math.tan(math.radians(stage["alpha2_deg"]))
    T01 = stage["T0_in_K"]
    n_b = stage[f"n_blades_{kind}"]
    c_mean = stage[f"chord_{kind}_mm"]
    c_hub = 2.0 * c_mean / (1.0 + CHORD_TAPER)
    tc_root, tc_tip = ((TC_ROOT_R, TC_TIP_R) if kind == "rotor"
                       else (TC_ROOT_S, TC_TIP_S))

    sections = []
    for f in np.linspace(0.0, 1.0, n_span):
        r = r_hub + f * (r_tip - r_hub)
        Cu1 = Cu1_m * r_m / r
        Cu2 = Cu2_m * r_m / r
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


def _cq_blade_solid(cq, row: dict, z_center: float):
    """UN álabe de una fila, en coordenadas de máquina.

    Cada sección (marco de cuerda, centroide en el origen) se rota por su
    stagger firmado y se posa en el plano radial X = r (workplane "YZ":
    u→Y tangencial, v→Z axial); loft spline entre planos paralelos y
    traslación al z del eje de stacking de la fila."""
    secs = row["sections"]
    wp = None
    r_prev = None
    for sec in secs:
        g = math.radians(sec["stagger_deg"])
        cg, sg = math.cos(g), math.sin(g)
        pts2 = [(px * sg + py * cg,            # u → Y (tangencial)
                 px * cg - py * sg)            # v → Z (axial, staggered)
                for px, py in sec["points"]]
        r = float(sec["r_mm"])
        if wp is None:
            wp = cq.Workplane("YZ", origin=(r, 0.0, 0.0))
        else:
            wp = wp.workplane(offset=r - r_prev)
        wp = wp.polyline(pts2).close()
        r_prev = r
    return wp.loft(ruled=False).translate((0.0, 0.0, z_center))


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
        solid = solid.union(_cq_blade_solid(cq, rot, zc))
        parts[f"RotorStage{i + 1}"] = solid

    lip = 4.0
    for i, st in enumerate(stages):
        z0 = (tip[0][0] - lip) if i == 0 else split_after(i - 1)
        z1 = (tip[-1][0] + lip) if i == n_st - 1 else split_after(i)
        # el casco vive a tip+holgura (annulus pasante)
        shell_line = [[z, r + clr] for z, r in tip]
        solid = _cq_line_ring(cq, shell_line, wall, z0, z1, outward=True)
        stat = st["stator"]
        solid = solid.union(_cq_blade_solid(cq, stat, stat["z_center_mm"]))
        if i == 0 and contract.get("igv"):
            igv = contract["igv"]
            solid = solid.union(_cq_blade_solid(cq, igv, igv["z_center_mm"]))
        if i == n_st - 1 and contract.get("ogv"):
            ogv = contract["ogv"]
            solid = solid.union(_cq_blade_solid(cq, ogv, ogv["z_center_mm"]))
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


def export_step(contract: dict, outdir: str, mode: str = "parts") -> list[str]:
    """Exporta STEP del diseño. Devuelve las rutas escritas ([] sin
    cadquery — degradación silenciosa, dependencia opcional `step`).

    mode:
      "blade0"   — UN álabe del primer rotor (rápido; compatibilidad).
      "parts"    — una pieza STEP por parte física (nombres de los STL)
                   con un álabe de muestra + parts/README.txt (default).
      "assembly" — además, cq.Assembly con todas las piezas nombradas →
                   assembly.step.
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
    igv_secs = _row_sections(_pseudo_stage_igv(record), "stator", omega)
    ogv_secs = _row_sections(_pseudo_stage_ogv(record), "stator", omega)
    igv_ext = _row_extent_mm(igv_secs)
    ogv_ext = _row_extent_mm(ogv_secs)
    # 2º: colocación secuencial sin solapes (IGV al frente, OGV al final)
    rows = _row_geometry(record, extents, igv_ext, ogv_ext)
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
    igv_row = _vane_row("igv", igv_secs, int(st_tab[0]["n_blades_stator"]))
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
        # fillet de raíz de álabe: el MISMO radio que usa el K_t de
        # Peterson en structures_core (fase 7 — fuente única)
        blade_fillet_r_mm=_sc.BLADE_FILLET_R_MM,
        # puerto de sangrado (boss radial con barreno en la carcasa)
        bleed_stage=max(n_st // 2, 1),   # etapa (1-based) del puerto
        bleed_hole_d_mm=18.0, bleed_boss_d_mm=32.0, bleed_boss_h_mm=14.0,
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
        w.writerow(["vano IGV", contract["igv"]["n_blades"], 0,
                    f"c{s0['chord_stator_mm']:.1f} x h{s0['h_blade_mm']:.1f}",
                    mat])
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
