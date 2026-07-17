"""
QUASAR Phy-AC · report_generator.py
===================================
Capa 5b del CEM: REPORTE HTML autocontenido del diseño (cero dependencias —
los gráficos son SVG generados en Python puro, sin matplotlib ni JS externo).

Secciones:
  1. Especificación de entrada y veredicto.
  2. Tabla de performance verificada + vector de diseño.
  3. Frente de Pareto verificado (SVG con hover por <title>).
  4. Tabla por etapa + carta de Smith (posición ψ–φ de cada etapa).
  5. Desglose de pérdidas por mecanismo (barras SVG).
  6. Triángulos de velocidad de la etapa 1 (SVG a escala) + annulus.
  7. Visualización matplotlib (visualization.py, PNG base64; opcional).
  8. Estructura del rotor (structures_core L0s, restricción dura).
  9. Condiciones de frontera CFD listas para copiar.
 10. Archivos de geometría y trazabilidad del proceso.
"""

from __future__ import annotations

import html
import json
import math
from datetime import datetime

import numpy as np

from physics_core import VAR_NAMES

ACCENT = "#0f6fff"
INK = "#1a1d27"
MUT = "#6b7280"
OK = "#0a8f5b"
BAD = "#c43d3d"


# ==========================================================================
# Helpers SVG
# ==========================================================================
def _svg_pareto(front: list[dict], spec, w=560, h=380) -> str:
    pts = [(p["PR"], p["eta_poly"], p["feasible"]) for p in front]
    if not pts:
        return "<p>(sin frente)</p>"
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    x0, x1 = min(xs) - 0.2, max(xs) + 0.2
    y0, y1 = min(ys) - 0.02, max(ys) + 0.02
    ml, mb, mt, mr = 56, 42, 16, 16

    def X(v):
        return ml + (v - x0) / (x1 - x0) * (w - ml - mr)

    def Y(v):
        return h - mb - (v - y0) / (y1 - y0) * (h - mb - mt)

    s = [f'<svg viewBox="0 0 {w} {h}" style="max-width:{w}px;width:100%">']
    for i in range(5):
        gx = x0 + i * (x1 - x0) / 4
        gy = y0 + i * (y1 - y0) / 4
        s.append(f'<line x1="{X(gx):.0f}" y1="{mt}" x2="{X(gx):.0f}" '
                 f'y2="{h - mb}" stroke="#e5e7eb"/>')
        s.append(f'<line x1="{ml}" y1="{Y(gy):.0f}" x2="{w - mr}" '
                 f'y2="{Y(gy):.0f}" stroke="#e5e7eb"/>')
        s.append(f'<text x="{X(gx):.0f}" y="{h - mb + 16}" font-size="11" '
                 f'text-anchor="middle" fill="{MUT}">{gx:.2f}</text>')
        s.append(f'<text x="{ml - 6}" y="{Y(gy) + 4:.0f}" font-size="11" '
                 f'text-anchor="end" fill="{MUT}">{gy:.3f}</text>')
    if x0 < spec.PR_target < x1:
        s.append(f'<line x1="{X(spec.PR_target):.0f}" y1="{mt}" '
                 f'x2="{X(spec.PR_target):.0f}" y2="{h - mb}" '
                 f'stroke="{ACCENT}" stroke-dasharray="5 4"/>')
        s.append(f'<text x="{X(spec.PR_target) + 4:.0f}" y="{mt + 12}" '
                 f'font-size="11" fill="{ACCENT}">PR* = '
                 f'{spec.PR_target}</text>')
    for pr, eta, feas in sorted(pts):
        c = OK if feas else BAD
        s.append(f'<circle cx="{X(pr):.1f}" cy="{Y(eta):.1f}" r="5" '
                 f'fill="{c}" fill-opacity="0.85" stroke="white">'
                 f'<title>PR={pr:.3f}  η={eta:.3f}  '
                 f'{"factible" if feas else "infactible"}</title></circle>')
    s.append(f'<text x="{(w + ml - mr) / 2:.0f}" y="{h - 6}" font-size="12" '
             f'text-anchor="middle" fill="{INK}">Relación de presiones PR '
             f'[-]</text>')
    s.append(f'<text x="14" y="{(h - mb + mt) / 2:.0f}" font-size="12" '
             f'text-anchor="middle" fill="{INK}" transform="rotate(-90 14 '
             f'{(h - mb + mt) / 2:.0f})">η politrópica [-]</text>')
    s.append("</svg>")
    return "".join(s)


def _machine_losses(rec: dict) -> dict:
    """Agrega el desglose por etapa a mecanismos de máquina [J/kg]."""
    agg = dict(rotor_perfil_endwall_choque=0.0, estator=0.0,
               holgura_de_punta=0.0)
    for s in rec.get("stage_table", []):
        L = s["losses"]
        agg["rotor_perfil_endwall_choque"] += L["rotor"]
        agg["estator"] += L["stator"]
        agg["holgura_de_punta"] += L["clearance"]
    return {k: v for k, v in agg.items() if v > 0}


def _svg_losses(losses: dict, w=560) -> str:
    if not losses:
        return "<p>(sin desglose)</p>"
    items = sorted(losses.items(), key=lambda kv: -kv[1])
    total = sum(v for _, v in items) or 1.0
    bh, gap, ml = 26, 10, 190
    h = len(items) * (bh + gap) + 14
    maxv = items[0][1] or 1.0
    names = dict(rotor_perfil_endwall_choque="Rotor (perfil+endwall+choque)",
                 estator="Estátor (perfil+endwall)",
                 holgura_de_punta="Holgura de punta")
    s = [f'<svg viewBox="0 0 {w} {h}" style="max-width:{w}px;width:100%">']
    for i, (k, v) in enumerate(items):
        y = i * (bh + gap) + 6
        bw = max((v / maxv) * (w - ml - 90), 2)
        s.append(f'<text x="{ml - 8}" y="{y + bh - 8}" font-size="12" '
                 f'text-anchor="end" fill="{INK}">{names.get(k, k)}</text>')
        s.append(f'<rect x="{ml}" y="{y}" width="{bw:.0f}" height="{bh}" '
                 f'rx="4" fill="{ACCENT}" fill-opacity="0.8"/>')
        s.append(f'<text x="{ml + bw + 6:.0f}" y="{y + bh - 8}" '
                 f'font-size="12" fill="{MUT}">{v / 1000:.2f} kJ/kg '
                 f'({100 * v / total:.0f}%)</text>')
    s.append("</svg>")
    return "".join(s)


def _svg_triangles(rec: dict, w=560, h=300) -> str:
    """Triángulos de velocidad de la etapa 1 (entrada y salida del rotor)."""
    st = rec.get("stage_table")
    if not st:
        return ""
    s1 = st[0]
    U = s1["U_m"]
    Cx = s1["Cx"]
    a1 = math.radians(s1["alpha1_deg"])
    a2 = math.radians(s1["alpha2_deg"])
    Cu1, Cu2 = Cx * math.tan(a1), Cx * math.tan(a2)
    if U <= 0:
        return ""
    sc = (w - 140) / max(U, 1.0)
    oy1, oy2 = h / 2 - 30, h - 40
    ox = 70

    def arrow(x1, y1, x2, y2, color, label, dy=-6):
        return (f'<line x1="{x1:.0f}" y1="{y1:.0f}" x2="{x2:.0f}" '
                f'y2="{y2:.0f}" stroke="{color}" stroke-width="2.5" '
                f'marker-end="url(#ar)"/>'
                f'<text x="{(x1 + x2) / 2 + 6:.0f}" '
                f'y="{(y1 + y2) / 2 + dy:.0f}"'
                f' font-size="12" fill="{color}">{label}</text>')

    s = [f'<svg viewBox="0 0 {w} {h}" style="max-width:{w}px;width:100%">',
         '<defs><marker id="ar" markerWidth="8" markerHeight="8" refX="6" '
         'refY="3" orient="auto"><path d="M0,0 L6,3 L0,6 z" '
         f'fill="{INK}"/></marker></defs>',
         f'<text x="{ox}" y="20" font-size="13" fill="{INK}" '
         f'font-weight="600">Etapa 1 — entrada del rotor (1) y salida (2), '
         'línea media</text>']
    # Estación 1: C1 (Cu1, Cx), W1 = C1 − U
    s.append(arrow(ox, oy1, ox + U * sc, oy1, MUT, f"U = {U:.0f} m/s"))
    s.append(arrow(ox, oy1, ox + Cu1 * sc, oy1 - Cx * sc, OK,
                   f"C₁ ({s1['C1']:.0f})"))
    s.append(arrow(ox + U * sc, oy1, ox + Cu1 * sc, oy1 - Cx * sc, ACCENT,
                   f"W₁ ({s1['W1']:.0f})"))
    # Estación 2
    s.append(arrow(ox, oy2, ox + U * sc, oy2, MUT, "U"))
    s.append(arrow(ox, oy2, ox + Cu2 * sc, oy2 - Cx * sc, OK,
                   f"C₂ ({s1['C2']:.0f})"))
    s.append(arrow(ox + U * sc, oy2, ox + Cu2 * sc, oy2 - Cx * sc, ACCENT,
                   f"W₂ ({s1['W2']:.0f})"))
    s.append(f'<text x="{ox}" y="{h - 8}" font-size="12" fill="{MUT}">'
             f'φ = {s1["phi"]:.3f}, ψ = {s1["psi"]:.3f}, '
             f'Rx = {s1["Rx"]:.2f}, β₁ = {s1["beta1_deg"]:.1f}°, '
             f'β₂ = {s1["beta2_deg"]:.1f}°</text>')
    s.append("</svg>")
    return "".join(s)


def _svg_smith(rec: dict, w=560, h=380) -> str:
    """Posición de cada etapa en el plano ψ–φ (región de Smith)."""
    st = rec.get("stage_table")
    if not st:
        return ""
    ml, mb, mt, mr = 56, 42, 16, 16
    x0, x1, y0, y1 = 0.2, 1.1, 0.1, 0.6

    def X(v):
        return ml + (v - x0) / (x1 - x0) * (w - ml - mr)

    def Y(v):
        return h - mb - (v - y0) / (y1 - y0) * (h - mb - mt)

    s = [f'<svg viewBox="0 0 {w} {h}" style="max-width:{w}px;width:100%">']
    # región de alta eficiencia de la carta de Smith (guía clásica)
    s.append(f'<rect x="{X(0.45):.0f}" y="{Y(0.42):.0f}" '
             f'width="{X(0.75) - X(0.45):.0f}" '
             f'height="{Y(0.20) - Y(0.42):.0f}" fill="{OK}" '
             f'fill-opacity="0.08" stroke="{OK}" stroke-dasharray="4 4"/>')
    s.append(f'<text x="{X(0.6):.0f}" y="{Y(0.42) - 6:.0f}" font-size="11" '
             f'text-anchor="middle" fill="{OK}">región η alta '
             '(Smith)</text>')
    for i in range(5):
        gx = x0 + i * (x1 - x0) / 4
        gy = y0 + i * (y1 - y0) / 4
        s.append(f'<line x1="{X(gx):.0f}" y1="{mt}" x2="{X(gx):.0f}" '
                 f'y2="{h - mb}" stroke="#e5e7eb"/>')
        s.append(f'<line x1="{ml}" y1="{Y(gy):.0f}" x2="{w - mr}" '
                 f'y2="{Y(gy):.0f}" stroke="#e5e7eb"/>')
        s.append(f'<text x="{X(gx):.0f}" y="{h - mb + 16}" font-size="11" '
                 f'text-anchor="middle" fill="{MUT}">{gx:.2f}</text>')
        s.append(f'<text x="{ml - 6}" y="{Y(gy) + 4:.0f}" font-size="11" '
                 f'text-anchor="end" fill="{MUT}">{gy:.2f}</text>')
    for q in st:
        fx = min(max(q["phi"], x0), x1)
        fy = min(max(q["psi"], y0), y1)
        s.append(f'<circle cx="{X(fx):.1f}" cy="{Y(fy):.1f}" r="7" '
                 f'fill="{ACCENT}" fill-opacity="0.85" stroke="white">'
                 f'<title>Etapa {q["stage"] + 1}: φ={q["phi"]:.3f} '
                 f'ψ={q["psi"]:.3f} η_tt={q["eta_tt"]:.3f}</title></circle>')
        s.append(f'<text x="{X(fx):.0f}" y="{Y(fy) + 4:.0f}" font-size="9" '
                 f'text-anchor="middle" fill="#fff">{q["stage"] + 1}</text>')
    s.append(f'<text x="{(w + ml - mr) / 2:.0f}" y="{h - 6}" font-size="12" '
             f'text-anchor="middle" fill="{INK}">Coeficiente de flujo φ '
             '[-]</text>')
    s.append(f'<text x="14" y="{(h - mb + mt) / 2:.0f}" font-size="12" '
             f'text-anchor="middle" fill="{INK}" transform="rotate(-90 14 '
             f'{(h - mb + mt) / 2:.0f})">Coef. de carga ψ [-]</text>')
    s.append("</svg>")
    return "".join(s)


def _svg_annulus(rec: dict, w=560, h=260) -> str:
    """Plano meridional: hub/tip con cajas de filas."""
    try:
        from geometry_generator import annulus_lines, _row_geometry
        ann = annulus_lines(rec)
        rows = _row_geometry(rec)
    except Exception:
        return ""
    hub, tip = ann["hub"], ann["tip"]
    z_max = max(p[0] for p in tip) or 1.0
    r_max = max(p[1] for p in tip) * 1.06
    ml, mb, mt, mr = 56, 36, 14, 14

    def X(z):
        return ml + z / z_max * (w - ml - mr)

    def Y(r):
        return h - mb - r / r_max * (h - mb - mt)

    s = [f'<svg viewBox="0 0 {w} {h}" style="max-width:{w}px;width:100%">']
    for name, line, color in (("tip", tip, INK), ("hub", hub, MUT)):
        pts = " ".join(f"{X(z):.1f},{Y(r):.1f}" for z, r in line)
        s.append(f'<polyline points="{pts}" fill="none" stroke="{color}" '
                 'stroke-width="2"/>')
    st = {q["stage"]: q for q in rec["stage_table"]}
    for rw in rows:
        q = st[rw["stage"]]
        color = ACCENT if rw["kind"] == "rotor" else OK
        s.append(f'<rect x="{X(rw["z_le"]):.1f}" y="{Y(q["r_tip_mm"]):.1f}" '
                 f'width="{max(X(rw["z_te"]) - X(rw["z_le"]), 2):.1f}" '
                 f'height="{Y(q["r_hub_mm"]) - Y(q["r_tip_mm"]):.1f}" '
                 f'fill="{color}" fill-opacity="0.35">'
                 f'<title>{rw["kind"]} etapa {rw["stage"] + 1}</title>'
                 '</rect>')
    s.append(f'<text x="{(w + ml) / 2:.0f}" y="{h - 4}" font-size="12" '
             f'text-anchor="middle" fill="{INK}">z [mm] — rotor azul, '
             'estátor verde</text>')
    s.append(f'<text x="14" y="{(h - mb + mt) / 2:.0f}" font-size="12" '
             f'text-anchor="middle" fill="{INK}" transform="rotate(-90 14 '
             f'{(h - mb + mt) / 2:.0f})">r [mm]</text>')
    s.append("</svg>")
    return "".join(s)


# ==========================================================================
# Reporte principal
# ==========================================================================
def _row(k, v):
    return f"<tr><td>{html.escape(str(k))}</td><td><b>{v}</b></td></tr>"


def _stage_table_html(rec: dict) -> str:
    head = ("<tr><th>Etapa</th><th>φ</th><th>ψ</th><th>PR</th><th>η_tt</th>"
            "<th>DF_r</th><th>DF_s</th><th>SM</th><th>M_rel,tip</th>"
            "<th>h [mm]</th><th>Z_r/Z_s</th></tr>")
    rows = "".join(
        f"<tr><td>{s['stage'] + 1}</td><td>{s['phi']:.3f}</td>"
        f"<td>{s['psi']:.3f}</td><td>{s['PR']:.3f}</td>"
        f"<td>{s['eta_tt']:.3f}</td><td>{s['DF_rotor']:.3f}</td>"
        f"<td>{s['DF_stator']:.3f}</td>"
        f"<td style='color:{OK if s['SM'] > 0.1 else BAD}'>{s['SM']:.3f}</td>"
        f"<td>{s['M_rel_tip']:.3f}</td><td>{s['h_blade_mm']:.1f}</td>"
        f"<td>{s['n_blades_rotor']}/{s['n_blades_stator']}</td></tr>"
        for s in rec.get("stage_table", []))
    return f"<table>{head}{rows}</table>"


def _figures_section(theta: np.ndarray, rec: dict,
                     figures_dir: str | None, log,
                     material: str | None = None) -> str:
    try:
        import visualization
        figs = visualization.generate_figures(theta, rec,
                                              outdir=figures_dir, log=log,
                                              material=material)
    except Exception as e:
        log(f"      ⚠ visualización omitida: {type(e).__name__}: {e}")
        figs = {}
    if not figs:
        return ("<h2>7 · Visualización matplotlib</h2>"
                "<p class='mut'>matplotlib no disponible — sección omitida "
                "(<code>pip install matplotlib</code> para habilitarla).</p>")
    cards = "".join(
        f'<div class="card" style="margin-bottom:20px">'
        f'<h3 style="font-size:15px;margin-top:2px">'
        f'{html.escape(f["title"])}</h3>'
        f'<img src="data:image/png;base64,{f["b64"]}" '
        f'alt="{html.escape(name)}" style="width:100%;height:auto"/>'
        f'<p class="mut" style="margin-bottom:2px">'
        f'{html.escape(f["caption"])}</p></div>'
        for name, f in figs.items())
    return f"<h2>7 · Visualización matplotlib</h2>{cards}"


def _structural_section(structural: dict | None, log) -> str:
    if not structural:
        return ""
    ok_chip = ('<span style="color:#fff;background:%s;padding:3px 10px;'
               'border-radius:6px;font-weight:600">%s</span>')
    chip = (ok_chip % (OK, "MÁRGENES OK") if structural["feasible_struct"]
            else ok_chip % (BAD, "MARGEN VIOLADO"))
    labels = ["Fluencia a sobrevelocidad", "Margen de estallido",
              "Límite AN²", "Raíz de álabe (K_t)", "Campbell (paso de álabes)"]
    dyn_rows = ""
    if "campbell" in structural:
        f1s = " · ".join(f"{p['f1_Hz']:.0f}"
                         for p in structural["campbell"]["per_row"])
        dyn_rows = (
            _row("f₁ flap por etapa (Southwell)", f"{f1s} Hz")
            + _row("Margen de Campbell mín (paso de álabes)",
                   f"{structural['campbell_margin_min']:.2f} "
                   f"(≥ {structural['params']['campbell_band']:.2f})")
            + _row("K_t raíz (Peterson, fillet "
                   f"{structural['blade_fillet_r_mm']:.1f} mm)",
                   f"{structural['k_t_root']:.2f}")
            + _row("σ_alt admisible (Goodman, peor etapa)",
                   f"{structural['sigma_alt_allow_MPa']:.0f} MPa")
            + _row("Cribado de flutter (peor V*/1.4)",
                   f"margen {structural['flutter_margin_min']:.2f} "
                   "(métrica, no g)"))
    rows = "".join(
        _row("Material", structural["material"])
        + _row("σ_vm máx (disco, peor etapa)",
               f"{structural['sigma_vm_max_MPa']:.0f} MPa "
               f"(adm. {structural['sigma_allow_MPa']:.0f} MPa, etapa "
               f"{structural['stage_worst_yield'] + 1})")
        + _row(f"Fluencia a {structural['params']['overspeed']:.0%} N",
               f"{structural['yield_ratio_overspeed']:.2f} de σ_y (≤ 1)")
        + _row("Margen de estallido",
               f"{structural['burst_margin']:.2f} "
               f"(mín {structural['burst_margin_min']:.2f})")
        + _row("AN²", f"{structural['AN2_in2rpm2']:.2e} in²·rpm² "
               f"(máx {structural['AN2_max_m2_rpm2'] / 6.4516e-4:.2e})")
        + _row("Raíz de álabe (con K_t)",
               f"{structural['blade_ratio']:.2f} de σ_y (≤ 1)")
        + dyn_rows
        + _row("Masa del rotor (discos)",
               f"{structural['rotor_mass_kg']:.1f} kg"))
    g_rows = "".join(
        f"<tr><td>{html.escape(lbl)}</td>"
        f"<td><b style='color:{OK if g <= 0 else BAD}'>{g:+.3f}</b></td></tr>"
        for lbl, g in zip(labels, structural["g_struct"]))
    warns = "".join(f"<li>{html.escape(w)}</li>"
                    for w in structural["structural_warnings"])
    warns_html = (f"<p class='mut'>Diagnósticos:</p><ul class='mut'>{warns}"
                  "</ul>") if warns else ""
    fig_html = ""
    try:
        import visualization
        b64 = visualization.structural_figure_b64(structural, log=log)
        if b64:
            fig_html = (f'<div class="card"><img src="data:image/png;'
                        f'base64,{b64}" alt="disk_stress" '
                        f'style="width:100%;height:auto"/></div>')
    except Exception as e:
        log(f"      ⚠ figura estructural omitida: {type(e).__name__}: {e}")
    return (f"<h2>8 · Estructura del rotor "
            f"<span class='mut'>(discos por etapa L0s, restricción del "
            f"optimizador)</span></h2>"
            f"<p>{chip}</p><div class='grid'><div><table>{rows}</table>"
            f"</div><div><h3 style='font-size:15px'>g_struct ≤ 0</h3>"
            f"<table>{g_rows}</table></div></div>{warns_html}{fig_html}"
            "<p class='mut'>El vector g_struct entra a la dominancia "
            "restringida del optimizador junto con las restricciones "
            "aerodinámicas: el diseño recomendado ya lo satisface por "
            "construcción. Solver de disco validado contra Timoshenko; "
            "materiales con valores típicos de handbook — sustituir por "
            "specs certificadas para diseño real.</p>")


def generate_report(spec, result: dict, path: str,
                    geometry_manifest: list | None = None,
                    cfx_bcs: dict | None = None,
                    figures_dir: str | None = None,
                    include_figures: bool = True,
                    structural: dict | None = None,
                    log=print) -> str:
    rec = result["record"]
    theta = np.asarray(result["theta"], dtype=float)
    front = result.get("pareto_front", [])
    hist = result.get("history", [])
    feas = rec.get("feasible", False)
    verdict = ("DISEÑO FACTIBLE VERIFICADO" if feas
               else "MEJOR CANDIDATO (revisar restricciones)")
    vcolor = OK if feas else BAD

    perf_rows = "".join(
        _row(lbl, f"{rec[k]:.3f}{u}") for k, lbl, u in [
            ("PR", "Relación de presiones PR_tt", ""),
            ("eta_poly", "Eficiencia politrópica η_poly", ""),
            ("eta_isen", "Eficiencia isentrópica η_tt", ""),
            ("U_tip", "Velocidad de punta U_tip", " m/s"),
            ("Mu", "Mach periférico M_u", ""),
            ("M_rel_tip1", "Mach relativo punta R1", ""),
            ("min_SM", "Margen de bombeo mín (Koch)", ""),
            ("length_mm", "Longitud axial", " mm"),
            ("r_tip_in_mm", "Radio de punta entrada", " mm"),
        ] if k in rec) \
        + _row("Etapas", rec.get("n_stages", "?")) \
        + _row("Potencia de eje", f"{rec.get('power_W', 0) / 1e6:.2f} MW")

    geo_rows = "".join(_row(n, f"{v:.3f}") for n, v in zip(VAR_NAMES, theta))

    hist_rows = "".join(
        f"<tr><td>{h['round']}</td><td>{h['PR']:.3f}</td>"
        f"<td>{h['eta']:.3f}</td><td>{h['Mu']:.3f}</td>"
        f"<td>{h['n_evals']}</td>"
        f"<td>{h['surrogate'].get('r2', {}).get('PR', 0):.3f} / "
        f"{h['surrogate'].get('r2', {}).get('eta_poly', 0):.3f}</td>"
        f"<td>{h['surrogate'].get('coverage_2sigma', 0):.2f}</td></tr>"
        for h in hist)

    figures_html = (_figures_section(
        theta, rec, figures_dir, log,
        material=(structural or {}).get("material"))
        if include_figures else "")
    structural_html = _structural_section(structural, log)

    cfx_html = ""
    if cfx_bcs:
        cfx_html = ("<h2>9 · Condiciones de frontera CFD</h2>"
                    "<pre>" + html.escape(json.dumps(cfx_bcs, indent=2))
                    + "</pre>")

    files_html = ""
    if geometry_manifest:
        items = "".join(f"<li><code>{html.escape(str(v))}</code></li>"
                        for v in geometry_manifest)
        files_html = f"<h2>10 · Archivos de geometría</h2><ul>{items}</ul>"

    n_feas = sum(1 for p in front if p["feasible"])
    doc = f"""<!DOCTYPE html><html lang="es"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Quasar Phy-AC — Reporte de diseño</title>
<style>
 body{{font-family:'Segoe UI',system-ui,sans-serif;color:{INK};margin:0;
      background:#f7f8fa}}
 .wrap{{max-width:920px;margin:0 auto;padding:32px 20px 60px}}
 h1{{font-size:26px;margin:.2em 0}} h2{{font-size:18px;margin-top:36px;
   border-bottom:2px solid #e5e7eb;padding-bottom:6px}}
 .verdict{{display:inline-block;padding:6px 14px;border-radius:8px;
   color:#fff;background:{vcolor};font-weight:600}}
 table{{border-collapse:collapse;width:100%;
   background:#fff;border-radius:8px;overflow:hidden;
   box-shadow:0 1px 3px rgba(0,0,0,.07)}}
 td,th{{padding:7px 12px;border-bottom:1px solid #eef0f3;font-size:14px;
   text-align:left}}
 .mut{{color:{MUT};font-size:13px}} pre{{background:#fff;padding:14px;
   border-radius:8px;overflow:auto;font-size:13px;
   box-shadow:0 1px 3px rgba(0,0,0,.07)}}
 .grid{{display:grid;grid-template-columns:1fr 1fr;gap:24px}}
 @media(max-width:760px){{.grid{{grid-template-columns:1fr}}}}
 .card{{background:#fff;border-radius:10px;padding:16px;
   box-shadow:0 1px 3px rgba(0,0,0,.07)}}
</style></head><body><div class="wrap">
<p class="mut">Quasar Solutions Lab · Phy-AC v0.1 ·
 {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
<h1>Reporte de diseño autónomo — Compresor axial multietapa</h1>
<p><span class="verdict">{verdict}</span></p>

<h2>1 · Especificación</h2>
<table>{_row("PR objetivo", spec.PR_target)}
{_row("Flujo másico", f"{spec.massflow} kg/s")}
{_row("RPM máx", f"{spec.RPM_max:,.0f}")}
{_row("U_tip máx", f"{spec.U_tip_max} m/s")}
{_row("r_tip máx", f"{spec.r_tip_max_mm} mm")}
{_row("Etapas máx", spec.n_stages_max)}
{_row("Entrada", f"{spec.T0_in} K · {spec.P0_in:,.0f} Pa")}</table>

<h2>2 · Resultado verificado</h2>
<div class="grid"><div>
<h3 style="font-size:15px">Performance (evaluación física, fuente:
 <code>{html.escape(str(rec.get('source', '')))}</code>)</h3>
<table>{perf_rows}</table></div><div>
<h3 style="font-size:15px">Vector de diseño θ</h3>
<table>{geo_rows}</table></div></div>
<p class="mut">Los valores de esta sección provienen de evaluaciones
físicas reales (no del surrogate); la incertidumbre del modelo aplica a
las predicciones de búsqueda, no a este punto verificado.</p>

<h2>3 · Frente de Pareto verificado
 <span class="mut">({len(front)} puntos, {n_feas} factibles)</span></h2>
<div class="card">{_svg_pareto(front, spec)}</div>

<h2>4 · Etapas y carta de Smith</h2>
{_stage_table_html(rec)}
<div class="card" style="margin-top:16px">{_svg_smith(rec)}</div>

<h2>5 · Desglose de pérdidas — por qué η es la que es</h2>
<div class="card">{_svg_losses(_machine_losses(rec))}</div>

<h2>6 · Triángulos de velocidad y annulus</h2>
<div class="card">{_svg_triangles(rec)}</div>
<div class="card" style="margin-top:16px">{_svg_annulus(rec)}</div>

{figures_html}
{structural_html}
{cfx_html}
{files_html}

<h2>11 · Trazabilidad del proceso</h2>
<table><tr><th>Ronda</th><th>PR</th><th>η</th><th>Mu</th>
<th>Evals acum.</th><th>R² (PR/η)</th><th>Cobertura ±2σ</th></tr>
{hist_rows}</table>
<p class="mut">Semilla 71 · reproducible bit a bit ·
 dataset completo en el checkpoint JSON.</p>
</div></body></html>"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(doc)
    return path
