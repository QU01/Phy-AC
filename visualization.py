"""
QUASAR Phy-AC · visualization.py
================================
Capa 5b (apoyo): figuras matplotlib del meanline axial, embebidas en el
reporte como PNG base64 y opcionalmente escritas a <outdir>/*.png.

Dependencia OPCIONAL (extra `viz`): si matplotlib no está instalado,
`generate_figures` devuelve {} y el reporte omite la sección con una nota.

Figuras:
  map            speedlines PR(ṁ) con línea de bombeo y choke (L0)
  stage_loading  ψ / DF / SM por etapa (barras)
  annulus        plano meridional hub/tip con cajas de filas
  sections       secciones del rotor 1 (hub/media/punta) staggered
  campbell       f₁(N) por etapa vs engine orders (fase 7; requiere
                 material — structures_core.blade_modes)
"""

from __future__ import annotations

import base64
import io
import math
import os

import numpy as np


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _fig_to_b64(fig, outdir: str | None, name: str) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        fig.savefig(os.path.join(outdir, f"{name}.png"), dpi=110,
                    bbox_inches="tight")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def generate_figures(theta: np.ndarray, rec: dict,
                     outdir: str | None = None, log=print,
                     material: str | None = None) -> dict:
    """Devuelve {nombre: {title, b64, caption}}. Nunca lanza por figura.

    `material`: aleación del rotor (activa la figura de Campbell de la
    fase 7; None = omitirla)."""
    try:
        plt = _mpl()
    except Exception:
        return {}
    theta = np.asarray(theta, dtype=float)
    figs: dict = {}

    # ---- mapa del compresor ----------------------------------------------
    try:
        from physics_core import compressor_map
        mp = compressor_map(theta)
        fig, ax = plt.subplots(figsize=(7, 4.6))
        surge_pts = []
        for sl in mp["speedlines"]:
            pts = [p for p in sl["points"] if not p["unstable"]]
            m = [p["mdot"] for p in pts]
            pr = [p["PR"] for p in pts]
            ax.plot(m, pr, "-o", ms=3,
                    label=f"N={sl['speed_frac']:.2f}")
            if sl["mdot_surge"] is not None:
                p0 = next(p for p in pts if p["mdot"] == sl["mdot_surge"])
                surge_pts.append((p0["mdot"], p0["PR"]))
        if len(surge_pts) > 1:
            ax.plot(*zip(*surge_pts), "r--", lw=2, label="línea de bombeo")
        ax.plot(mp["design"]["mdot"], mp["design"]["PR"], "k*", ms=14,
                label="diseño")
        ax.set_xlabel("ṁ [kg/s]")
        ax.set_ylabel("PR_tt [-]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        figs["map"] = dict(
            title="Mapa de operación (L0)",
            b64=_fig_to_b64(fig, outdir, "map"),
            caption="Speedlines con geometría congelada; bombeo por SM de "
                    "Koch ≤ 0 y rama de pendiente positiva. Tendencia L0 — "
                    "recalibrar contra CFD/banco.")
    except Exception as e:
        log(f"      ⚠ figura map: {type(e).__name__}: {e}")

    # ---- carga por etapa ---------------------------------------------------
    try:
        st = rec["stage_table"]
        idx = np.arange(len(st)) + 1
        fig, ax = plt.subplots(figsize=(7, 3.8))
        w = 0.27
        ax.bar(idx - w, [s["psi"] for s in st], w, label="ψ")
        ax.bar(idx, [s["DF_rotor"] for s in st], w, label="DF rotor")
        ax.bar(idx + w, [s["SM"] for s in st], w, label="SM (Koch)")
        ax.axhline(0.10, color="r", ls="--", lw=1, alpha=0.6)
        ax.set_xticks(idx)
        ax.set_xlabel("etapa")
        ax.grid(alpha=0.3, axis="y")
        ax.legend(fontsize=9)
        figs["stage_loading"] = dict(
            title="Carga, difusión y margen por etapa",
            b64=_fig_to_b64(fig, outdir, "stage_loading"),
            caption="ψ distribuido por psi_mid/psi_slope; línea roja = "
                    "margen mínimo de bombeo exigido (g).")
    except Exception as e:
        log(f"      ⚠ figura stage_loading: {type(e).__name__}: {e}")

    # ---- annulus meridional -----------------------------------------------
    try:
        from geometry_generator import annulus_lines, _row_geometry
        ann = annulus_lines(rec)
        rows = _row_geometry(rec)
        st = {q["stage"]: q for q in rec["stage_table"]}
        fig, ax = plt.subplots(figsize=(7, 3.6))
        for line, c in ((ann["hub"], "k"), (ann["tip"], "k")):
            z = [p[0] for p in line]
            r = [p[1] for p in line]
            ax.plot(z, r, c, lw=2)
        for rw in rows:
            q = st[rw["stage"]]
            color = "#0f6fff" if rw["kind"] == "rotor" else "#0a8f5b"
            ax.fill_between([rw["z_le"], rw["z_te"]],
                            q["r_hub_mm"], q["r_tip_mm"],
                            color=color, alpha=0.35)
        ax.set_xlabel("z [mm]")
        ax.set_ylabel("r [mm]")
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.3)
        figs["annulus"] = dict(
            title="Plano meridional",
            b64=_fig_to_b64(fig, outdir, "annulus"),
            caption="Annulus a r_medio constante; rotor azul, estátor "
                    "verde. La altura de álabe cae con la compresión.")
    except Exception as e:
        log(f"      ⚠ figura annulus: {type(e).__name__}: {e}")

    # ---- secciones del rotor 1 ---------------------------------------------
    try:
        from geometry_generator import _row_sections
        omega = float(theta[1]) * 2 * math.pi / 60.0
        secs = _row_sections(rec["stage_table"][0], "rotor", omega)
        fig, ax = plt.subplots(figsize=(6, 6))
        for sec in (secs[0], secs[len(secs) // 2], secs[-1]):
            g = math.radians(sec["stagger_deg"])
            cg, sg = math.cos(g), math.sin(g)
            p = np.array(sec["points"])
            x = p[:, 0] * cg - p[:, 1] * sg
            y = p[:, 0] * sg + p[:, 1] * cg
            ax.plot(np.append(x, x[0]), np.append(y, y[0]),
                    label=f"r={sec['r_mm']:.0f} mm ({sec['profile']})")
        ax.set_aspect("equal")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
        ax.set_xlabel("axial [mm]")
        ax.set_ylabel("tangencial [mm]")
        figs["sections"] = dict(
            title="Secciones del rotor 1 (hub/media/punta, staggered)",
            b64=_fig_to_b64(fig, outdir, "sections"),
            caption="Free vortex: el stagger crece y la comba cae hacia la "
                    "punta; DCA donde M_entrada > 0.8.")
    except Exception as e:
        log(f"      ⚠ figura sections: {type(e).__name__}: {e}")

    # ---- diagrama de Campbell (fase 7) -------------------------------------
    if material is not None:
        try:
            from structures_core import blade_modes, CAMPBELL_BAND
            st = rec["stage_table"]
            RPM_d = float(theta[1])
            n_pts = 40
            Ns = np.linspace(0.65 * RPM_d, 1.10 * RPM_d, n_pts)
            fig, ax = plt.subplots(figsize=(7, 4.6))
            eo_max_hz = 0.0
            for s in st:
                f1s = [blade_modes(s, material, N, kind="rotor")["f1_Hz"]
                       for N in Ns]
                ax.plot(Ns, f1s, lw=2,
                        label=f"f₁ rotor {int(s['stage']) + 1}")
                eo_bp = int(s["n_blades_stator"]) * RPM_d / 60.0
                eo_max_hz = max(eo_max_hz, max(f1s))
            # abanico de engine orders k=1..6 + paso de álabes (min y máx)
            for k in range(1, 7):
                ax.plot(Ns, k * Ns / 60.0, color="0.65", lw=0.8)
                ax.annotate(f"{k}E", (Ns[-1], k * Ns[-1] / 60.0),
                            fontsize=7, color="0.4")
            n_bp = sorted({int(s["n_blades_stator"]) for s in st})
            for nb in (n_bp[0], n_bp[-1]) if n_bp else ():
                ax.plot(Ns, nb * Ns / 60.0, "m-", lw=1.0, alpha=0.7)
                ax.annotate(f"{nb}E (paso)", (Ns[0], nb * Ns[0] / 60.0),
                            fontsize=7, color="m")
            ax.axvspan(RPM_d * (1 - CAMPBELL_BAND),
                       RPM_d * (1 + CAMPBELL_BAND),
                       color="r", alpha=0.08,
                       label=f"±{CAMPBELL_BAND:.0%} N diseño")
            ax.axvline(RPM_d, color="r", ls="--", lw=1)
            ax.set_xlabel("N [rpm]")
            ax.set_ylabel("f [Hz]")
            ax.set_ylim(0, 1.4 * eo_max_hz)
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8, loc="upper left")
            figs["campbell"] = dict(
                title="Diagrama de Campbell (flap 1º, Southwell)",
                b64=_fig_to_b64(fig, outdir, "campbell"),
                caption="f₁(N) por etapa de rotor con rigidización "
                        "centrífuga vs engine orders k=1..6 (gris) y paso "
                        "de álabes (magenta, líneas fuera de escala hacia "
                        "arriba). El margen duro de g_struct es frente al "
                        "paso de álabes de los estátores vecinos a N de "
                        "diseño; los EO bajos son métrica reportada.")
        except Exception as e:
            log(f"      ⚠ figura campbell: {type(e).__name__}: {e}")

    return figs


def structural_figure_b64(structural: dict, log=print) -> str | None:
    """σ_vm de disco y raíz de álabe por etapa vs admisible."""
    try:
        plt = _mpl()
    except Exception:
        return None
    try:
        ps = structural["per_stage"]
        idx = np.arange(len(ps)) + 1
        fig, ax = plt.subplots(figsize=(7, 3.6))
        w = 0.35
        ax.bar(idx - w / 2, [p["sigma_vm_MPa"] for p in ps], w,
               label="σ_vm disco")
        ax.bar(idx + w / 2, [p["sigma_root_MPa"] for p in ps], w,
               label="σ raíz de álabe")
        ax.axhline(structural["sigma_allow_MPa"], color="r", ls="--",
                   label="σ_y derateada (peor etapa)")
        ax.set_xticks(idx)
        ax.set_xlabel("etapa")
        ax.set_ylabel("MPa")
        ax.grid(alpha=0.3, axis="y")
        ax.legend(fontsize=9)
        return _fig_to_b64(fig, None, "structural")
    except Exception as e:
        log(f"      ⚠ figura estructural: {type(e).__name__}: {e}")
        return None
