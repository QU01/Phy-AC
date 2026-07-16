#!/usr/bin/env python3
"""
QUASAR Phy-AC · phyac_cli.py
============================
Punto de entrada de PRODUCTO: el ingeniero escribe la especificación y recibe
el paquete completo — diseño verificado, frente de Pareto, geometría de filas
de álabes, reporte HTML y dataset acumulado.

Uso:
    python phyac_cli.py --pr 4.0 --mdot 25 --rpm-max 18000 \\
        --utip-max 460 --rtip-max 400 --rounds 5 --n-init 320 --outdir runs/x

    # asistente interactivo (pregunta los números de la espec):
    python phyac_cli.py --interactive
    python phyac_cli.py                  # sin argumentos = asistente

    # corrida rápida de humo:
    python phyac_cli.py --pr 4.0 --mdot 25 --quick

    # con espec en JSON:
    python phyac_cli.py --spec-file espec.json --outdir runs/cliente_x

Salidas en <outdir>/:
    phyac_run.json   checkpoint completo (espec, historial, dataset)
    dataset.csv      dataset físico acumulado (semilla del flywheel)
    report.html      reporte autocontenido (Capa 5)
    geometry/        axial_compressor.json, annulus.csv, stage_summary.csv
                     [, blade_stage0.step]
    *.stl            con --stl/--voxel (capa 5c, C#/PicoGK)
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys

import numpy as np

from physics_core import (Fidelity, register_hifi_pair, VAR_NAMES,
                          l1_available, l1_unavailable_reason,
                          set_cache_path, compressor_map)
from neural_optimizer import DesignSpec, design, OUT_KEYS
import geometry_generator
import report_generator
import structures_core
from geometry_generator import cfx_boundary_conditions
import cli_style as ui


# ----------------------------------------------------------- argumentos
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="phyac", add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Quasar Phy-AC — autonomous axial compressor "
        "design (spec → verified blade rows → STL).",
        epilog="examples:\n"
        "  phyac --interactive                guided wizard\n"
        "  phyac --pr 4 --mdot 25 --quick     quick smoke run (~2 min)\n"
        "  phyac --pr 4 --mdot 25 --rounds 5 --outdir runs/x\n"
        "  phyac --spec-file spec.json --map  spec from JSON + map\n")
    s = p.add_argument_group("specification")
    s.add_argument("--pr", type=float, default=4.0, help="Target PR (default 4.0)")
    s.add_argument("--mdot", type=float, default=25.0, help="kg/s (default 25)")
    s.add_argument("--rpm-max", type=float, default=20_000)
    s.add_argument("--utip-max", type=float, default=480.0, help="m/s")
    s.add_argument("--rtip-max", type=float, default=400.0, help="mm")
    s.add_argument("--nstages-max", type=int, default=8)
    s.add_argument("--power-max", type=float, default=None, help="W")
    s.add_argument("--t0", type=float, default=288.15, help="K")
    s.add_argument("--p0", type=float, default=101_325.0, help="Pa")
    s.add_argument("--material", type=str, default="Ti-6Al-4V",
                   choices=structures_core.list_materials(),
                   help="rotor material for the structural margins "
                        "(default Ti-6Al-4V)")
    s.add_argument("--spec-file", type=str, default=None,
                   help="JSON with the spec (overwrites flags)")
    b = p.add_argument_group("budget")
    b.add_argument("--rounds", type=int, default=5)
    b.add_argument("--n-init", type=int, default=320)
    b.add_argument("--batch-size", type=int, default=14)
    b.add_argument("--workers", type=int, default=1,
                   help="parallel processes for L1 SCM")
    b.add_argument("--fidelity", choices=["L0", "L1"], default="L1")
    b.add_argument("--quick", action="store_true",
                   help="quick smoke run (n_init=100, rounds=2, batch=8)")
    o = p.add_argument_group("outputs")
    o.add_argument("--outdir", type=str, default="phyac_out")
    o.add_argument("--map", action="store_true",
                   help="generate off-design map of the final design")
    o.add_argument("--hifi-pairs", type=str, default=None,
                   help="JSON with CFD calibration pairs")
    o.add_argument("--stl", action="store_true",
                   help="generate rotor/casing STL files (default 0.5mm voxel)")
    o.add_argument("--voxel", type=float, default=None,
                   help="voxel resolution in mm for STL generation "
                        "(e.g. 0.3, 0.5, 0.8); implies --stl")
    o.add_argument("--no-figures", action="store_true",
                   help="skip the matplotlib visualization section")
    u = p.add_argument_group("interface")
    u.add_argument("--interactive", "-i", action="store_true",
                   help="wizard that asks for specification interactively")
    u.add_argument("--no-color", action="store_true",
                   help="disable color and ASCII art")
    u.add_argument("--quiet", "-q", action="store_true",
                   help="milestones only")
    return p.parse_args(argv)


def build_spec(args) -> DesignSpec:
    kw = dict(PR_target=args.pr, massflow=args.mdot, RPM_max=args.rpm_max,
              U_tip_max=args.utip_max, r_tip_max_mm=args.rtip_max,
              n_stages_max=args.nstages_max,
              power_max_W=args.power_max, T0_in=args.t0, P0_in=args.p0,
              material=args.material)
    if args.spec_file:
        with open(args.spec_file) as f:
            kw.update(json.load(f))
    return DesignSpec(**kw)


# ----------------------------------------------------------- interactive
def _ask(prompt: str, default, cast=float):
    raw = input(ui.c("  ? ", "brand", bold=True)
                + prompt + ui.c(f" [{default}]: ", "dim"))
    raw = raw.strip()
    if not raw:
        return default
    try:
        return cast(raw)
    except ValueError:
        print(ui.warn(f"invalid value, using {default}"))
        return default


def run_wizard(args):
    print(ui.rule("Specification Wizard", "accent"))
    print(ui.c("  Press Enter to accept default values.\n", "dim"))
    args.pr = _ask("Target PR (pressure ratio)", args.pr)
    args.mdot = _ask("Mass flow [kg/s]", args.mdot)
    args.rpm_max = _ask("Max RPM", args.rpm_max)
    args.utip_max = _ask("Max tip speed [m/s]", args.utip_max)
    args.rtip_max = _ask("Max tip radius [mm]", args.rtip_max)
    args.nstages_max = _ask("Max stages", args.nstages_max, int)
    args.t0 = _ask("Inlet T0 [K]", args.t0)
    args.p0 = _ask("Inlet P0 [Pa]", args.p0)
    mats = structures_core.list_materials()
    raw = input(ui.c("  ? ", "brand", bold=True)
                + f"Rotor material ({', '.join(mats)})"
                + ui.c(f" [{args.material}]: ", "dim")).strip()
    if raw:
        if raw in mats:
            args.material = raw
        else:
            print(ui.warn(f"unknown material, using {args.material}"))
    mode = input(ui.c("  ? ", "brand", bold=True)
                 + "mode: [q]uick / [f]ull"
                 + ui.c(" [q]: ", "dim")).strip().lower()
    args.quick = not mode.startswith("f")
    print()
    return args


# ----------------------------------------------------------- live logger
class StyledLog:
    """Logger that colorizes and formats optimizer output."""

    def __init__(self, quiet=False):
        self.quiet = quiet

    def __call__(self, msg="") -> None:
        text = str(msg)
        stripped = text.strip()
        m = re.search(r"Round\s+(\d+)/(\d+)", stripped)
        if m:
            n, tot = int(m.group(1)), int(m.group(2))
            print()
            print(ui.rule(f"Round {n}/{tot}", "brand2"))
            print("  " + ui.progress_bar(n / tot, ui.c("exploring…", "dim")))
            return
        if "⚠" in text or "FAILED" in text:
            print("  " + ui.warn(stripped.lstrip("⚠ ").strip()))
            return
        if "Best verified" in stripped:
            print("  " + ui.ok(stripped))
            return
        if stripped and not self.quiet:
            print("  " + ui.c(stripped, "dim"))


# ----------------------------------------------------------- export
def export_dataset_csv(run_json: str, csv_path: str) -> int:
    if not os.path.exists(run_json):
        return 0
    with open(run_json) as f:
        data = json.load(f)
    rows = data.get("dataset", [])
    if not rows:
        return 0
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(VAR_NAMES + OUT_KEYS + ["g", "source"])
        for r in rows:
            w.writerow(list(r["theta"]) + [r[k] for k in OUT_KEYS]
                       + [json.dumps(r.get("g")), r.get("source")])
    return len(rows)


def generate_stls(contract_json: str, outdir: str, voxel: float) -> bool:
    import subprocess
    import shutil

    dotnet_bin = shutil.which("dotnet")
    if not dotnet_bin:
        standard_path = r"C:\Program Files\dotnet\dotnet.exe"
        if os.path.exists(standard_path):
            dotnet_bin = standard_path
        else:
            print(ui.warn("No se encontró 'dotnet' — instala .NET SDK 9 "
                          "para la capa 5c."))
            return False

    print(ui.c("\nCompilando AxialCompressorDesigner...", "dim"))
    try:
        res_build = subprocess.run(
            [dotnet_bin, "build", "AxialCompressorDesigner.sln", "-c",
             "Debug"], capture_output=True, text=True)
        if res_build.returncode != 0:
            print(ui.warn("Error compilando AxialCompressorDesigner.sln:"))
            print(res_build.stderr or res_build.stdout)
            return False
    except Exception as e:
        print(ui.warn(f"Error ejecutando compilación dotnet: {e}"))
        return False

    dll_path = os.path.join("AxialCompressorDesigner.Example", "bin",
                            "Debug", "net9.0-windows",
                            "AxialCompressorDesigner.Example.dll")
    if not os.path.exists(dll_path):
        print(ui.warn(f"No se encontró el ejecutable en: {dll_path}"))
        return False

    print(ui.c(f"Generando STLs (voxel={voxel} mm)... puede tardar unos "
               "segundos/minutos...", "accent"))
    try:
        res_run = subprocess.run(
            [dotnet_bin, dll_path, contract_json, outdir, str(voxel)],
            capture_output=True, text=True)
        if res_run.returncode != 0:
            print(ui.warn("Error ejecutando AxialCompressorDesigner:"))
            print(res_run.stderr or res_run.stdout)
            return False
        for line in res_run.stdout.splitlines():
            line_s = line.strip()
            if line_s and ("→" in line_s or "PhyAC" in line_s
                           or "Exported" in line_s or "[PhyAC]" in line_s):
                print("  " + ui.ok(line_s))
        return True
    except Exception as e:
        print(ui.warn(f"Error generando STLs: {e}"))
        return False


# ----------------------------------------------------------- main
def main(argv=None) -> int:
    # consolas/pipes Windows reportan cp1252 y truenan con η/₂/✓
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    args = parse_args(argv)
    ui.init(use_color=False if args.no_color else None)

    raw_argv = sys.argv[1:] if argv is None else list(argv)
    if not raw_argv:
        args.interactive = True

    print()
    print(ui.banner())
    print()

    if args.interactive:
        args = run_wizard(args)
    if args.quick:
        args.n_init, args.rounds, args.batch_size = 100, 2, 8
    os.makedirs(args.outdir, exist_ok=True)
    spec = build_spec(args)
    fidelity = Fidelity.L1 if args.fidelity == "L1" else Fidelity.L0

    if fidelity == Fidelity.L1 and not l1_available():
        print(ui.warn(f"turbo-design not available ({l1_unavailable_reason()})"))
        print(ui.c("    → entire run will be L0 (meanline), tagged as "
                   "'source'.", "dim"))
        fidelity = Fidelity.L0
        args.fidelity = "L0"

    set_cache_path(os.environ.get("PHYAC_CACHE")
                   or os.path.join(args.outdir, "phys_cache.jsonl"))

    if args.hifi_pairs:
        with open(args.hifi_pairs) as f:
            for pair in json.load(f):
                register_hifi_pair(np.array(pair["theta"]),
                                   {k: v for k, v in pair.items()
                                    if k != "theta"})
        print(ui.ok(f"calibration loaded: {args.hifi_pairs}"))

    budget = (f"{args.n_init} init + {args.rounds}×{args.batch_size} "
              "acquisition" + ("  (quick)" if args.quick else ""))
    print(ui.panel("Specification", [
        ("Target PR", ui.c(f"{spec.PR_target:.2f}", "val", bold=True)),
        ("Mass flow", f"{spec.massflow:.2f} kg/s"),
        ("Max RPM", f"{spec.RPM_max:,.0f}"),
        ("Max U_tip", f"{spec.U_tip_max:.0f} m/s"),
        ("Max r_tip", f"{spec.r_tip_max_mm:.0f} mm"),
        ("Max stages", f"{spec.n_stages_max}"),
        ("Inlet", f"T0={spec.T0_in:.1f} K   P0={spec.P0_in:,.0f} Pa"),
        ("Budget", ui.c(budget, "accent")),
        ("Fidelity", ui.c(args.fidelity, "brand")),
    ]))
    print()
    print(ui.rule("Autonomous optimization", "brand"))

    ckpt = os.path.join(args.outdir, "phyac_run.json")
    logger = StyledLog(quiet=args.quiet)
    result = design(spec, rounds=args.rounds, n_init=args.n_init,
                    batch_size=args.batch_size, fidelity=fidelity,
                    n_workers=args.workers, checkpoint_path=ckpt, log=logger)

    rec, theta = result["record"], np.asarray(result["theta"])

    # ---- Structural margins (structures_core L0s — hard constraint) ----
    structural = None
    if spec.material is not None:
        try:
            structural = structures_core.evaluate_structural(
                theta, rec, material=spec.material)
        except Exception as e:
            print(ui.warn(f"structural evaluation failed: {e}"))

    # ---- Layer 5: geometry + report + dataset ----
    print()
    print(ui.rule("Deliverables (Layer 5)", "brand"))
    geo_dir = os.path.join(args.outdir, "geometry")
    manifest = geometry_generator.generate(
        theta, rec, geo_dir, structural=structural,
        run_meta=dict(checkpoint=ckpt, rounds=args.rounds,
                      n_init=args.n_init, fidelity=args.fidelity))
    print("  " + ui.ok(f"geometry   → {geo_dir}/ ({len(manifest)} files)"))
    bcs = cfx_boundary_conditions(theta, rec)
    figs_dir = os.path.join(args.outdir, "figures")
    report = report_generator.generate_report(
        spec, result, os.path.join(args.outdir, "report.html"),
        geometry_manifest=manifest, cfx_bcs=bcs,
        figures_dir=figs_dir, include_figures=not args.no_figures,
        structural=structural, log=logger)
    print("  " + ui.ok(f"report     → {report}"))
    if not args.no_figures and os.path.isdir(figs_dir):
        n_figs = len([f for f in os.listdir(figs_dir) if f.endswith(".png")])
        if n_figs:
            print("  " + ui.ok(f"figures    → {figs_dir}/ ({n_figs} PNG)"))
    n_rows = export_dataset_csv(ckpt, os.path.join(args.outdir, "dataset.csv"))
    print("  " + ui.ok(f"dataset    → {n_rows} rows (flywheel)"))

    # ---- Generación de STL con C# (opcional) ----
    if args.stl or args.voxel is not None:
        voxel_val = args.voxel if args.voxel is not None else 0.5
        contract_path = os.path.join(geo_dir, "axial_compressor.json")
        if os.path.exists(contract_path):
            generate_stls(contract_path, args.outdir, voxel_val)
        else:
            print(ui.warn(f"Could not find {contract_path} to generate STLs."))

    # ---- Mapa off-design opcional ----
    if args.map:
        mp = compressor_map(theta)
        map_path = os.path.join(args.outdir, "map.csv")
        with open(map_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["speed_frac", "rpm", "mdot_kgs", "PR", "eta_poly",
                        "min_SM", "M_rel_tip1", "unstable", "stall",
                        "choke"])
            for sl in mp["speedlines"]:
                for q in sl["points"]:
                    w.writerow([sl["speed_frac"], sl["rpm"], q["mdot"],
                                q["PR"], q["eta_poly"], q["min_SM"],
                                q["M_rel_tip1"], int(q["unstable"]),
                                int(q["stall"]), int(q["choke"])])
        print("  " + ui.ok(f"map        → {map_path}"))
        print(ui.c("\n  Off-design map (L0 proxy):", "dim"))
        for sl in mp["speedlines"]:
            su, ch = sl["mdot_surge"], sl["mdot_choke"]
            rng_txt = (f"{su:.2f} – {ch:.2f} kg/s" if su is not None
                       else ui.c("no stable range", "warn"))
            print(f"    {ui.c('N=%.2f' % sl['speed_frac'], 'key')} "
                  f"({sl['rpm']:,.0f} RPM): {rng_txt}")

    # ---- Final result ----
    feasible = rec.get("feasible")
    status = (ui.c(" FEASIBLE ", "ok", bold=True) if feasible
              else ui.c(" INFEASIBLE ", "err", bold=True))
    feas = sum(1 for q in result["pareto_front"] if q["feasible"])
    print()
    rows = [
        ("PR", ui.c(f"{rec['PR']:.3f}", "val", bold=True)
         + ui.c(f"  (target {spec.PR_target:.2f})", "dim")),
        ("polytropic η", ui.c(f"{rec['eta_poly']:.3f}", "val", bold=True)),
        ("Stages", f"{rec.get('n_stages', '?')}"),
        ("U_tip", f"{rec['U_tip']:.0f} m/s"),
        ("Status", status),
        ("Pareto", f"{len(result['pareto_front'])} points ({feas} feasible)"),
        ("Checkpoint", ui.c(ckpt, "dim")),
    ]
    if structural:
        s_ok = structural["feasible_struct"]
        rows.insert(5, ("Structure",
                        (ui.c(" MARGINS OK ", "ok", bold=True) if s_ok
                         else ui.c(" MARGIN VIOLATED ", "err", bold=True))
                        + ui.c(f"  {structural['material']} · "
                               f"burst {structural['burst_margin']:.2f} · "
                               f"AN² {structural['AN2_in2rpm2']:.1e} "
                               "in²rpm²", "dim")))
    print(ui.panel("Recommended design", rows,
                   style="ok" if feasible else "warn"))
    print()
    print("  " + ui.c("Open report:", "dim") + " "
          + ui.c(report, "brand", bold=True))
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
