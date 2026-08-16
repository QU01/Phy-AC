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

from physics_core import (Fidelity, register_hifi_pair, VAR_NAMES, NDIM,
                          BOUNDS_LO, BOUNDS_HI,
                          l1_available, l1_unavailable_reason,
                          set_cache_path, compressor_map,
                          SM_FLOW_MIN as pc_SM_MIN)
from neural_optimizer import (DesignSpec, design, OUT_KEYS,
                              AutonomousAxialDesigner, evaluate_design,
                              pareto_from_checkpoint)
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
    s.add_argument("--vortex-n", type=float, default=None, metavar="N",
                   help="vortex exponent Cu(r) = Cu_m*(r/r_m)^N. -1 = free "
                        "vortex (default, classic); -0.5..0 = controlled "
                        "vortex, which is what a low hub-to-tip FAN needs — "
                        "free vortex collapses its hub reaction")
    s.add_argument("--m-exit-max", type=float, default=None, metavar="M",
                   help="max absolute exit Mach (default 0.55, diffuser/"
                        "combustor-friendly). A fan discharging into a "
                        "bypass duct runs 0.6-0.7")
    s.add_argument("--sm-flow-min", type=float, default=None, metavar="F",
                   help="minimum working-line surge margin in mass flow "
                        "(default 0.15 = 15%%)")
    s.add_argument("--spec-file", type=str, default=None,
                   help="JSON with the spec (overwrites flags)")
    s.add_argument("--fix", action="append", default=None,
                   metavar="VAR=VALUE",
                   help="pin a design variable and optimize the rest "
                        "(repeatable), e.g. --fix n_stages=5 --fix "
                        "phi1=0.60. Vars: " + ", ".join(VAR_NAMES[:10]))
    b = p.add_argument_group("budget")
    b.add_argument("--rounds", type=int, default=5)
    b.add_argument("--n-init", type=int, default=320)
    b.add_argument("--batch-size", type=int, default=14)
    b.add_argument("--workers", type=int, default=1,
                   help="parallel processes for L1 SCM")
    b.add_argument("--fidelity", choices=["L0", "L1"], default="L1")
    b.add_argument("--quick", action="store_true",
                   help="quick smoke run (n_init=100, rounds=2, batch=8)")
    b.add_argument("--seed", type=int, default=None,
                   help="random seed of the whole run (default 71; a "
                        "resumed run keeps its checkpoint's seed)")
    b.add_argument("--resume", nargs="?", const="__auto__", default=None,
                   metavar="CHECKPOINT",
                   help="warm-start from a phyac_run.json checkpoint "
                        "(without a value: <outdir>/phyac_run.json); the "
                        "spec is taken from the checkpoint")
    d = p.add_argument_group("direct evaluation (no optimization)")
    d.add_argument("--eval-theta", type=str, default=None,
                   metavar="V1,V2,...",
                   help="evaluate ONE given design and emit the full "
                        "deliverables. Accepts 10 design values "
                        f"({', '.join(VAR_NAMES[:10])}), 12 (those 10 + "
                        "phi_slope,Rx_slope), 13 (legacy full theta incl. "
                        f"T0,P0,mdot) or {NDIM} (full theta, op point at "
                        "positions 11-13)")
    d.add_argument("--per-stage", type=str, default=None, metavar="FILE",
                   help="expert mode (with --eval-theta): JSON with "
                        "per-stage overrides {\"phi\": [...], \"psi\": "
                        "[...], \"Rx\": [...]} (lists of n_stages values) "
                        "replacing the linear-slope distributions; "
                        "meanline L0 only")
    d.add_argument("--from-checkpoint", type=str, default=None,
                   help="checkpoint for --list-pareto / --pareto-pick "
                        "(default <outdir>/phyac_run.json)")
    d.add_argument("--list-pareto", action="store_true",
                   help="print the verified Pareto front of a checkpoint "
                        "and exit")
    d.add_argument("--pareto-pick", type=int, default=None, metavar="N",
                   help="regenerate the deliverables for point #N of the "
                        "checkpoint's Pareto front (see --list-pareto); "
                        "written to <outdir>/pareto_NN/")
    o = p.add_argument_group("outputs")
    o.add_argument("--outdir", type=str, default="phyac_out")
    o.add_argument("--map", action="store_true",
                   help="generate off-design map of the final design")
    o.add_argument("--map-vsv", action="store_true",
                   help="close the front variable stators along each "
                        "part-speed line (auto schedule) — without it, "
                        "low-speed lines of PR>4 machines are unphysical "
                        "fixed-geometry extrapolation; implies --map")
    o.add_argument("--hifi-pairs", type=str, default=None,
                   help="JSON with CFD calibration pairs")
    o.add_argument("--stl", action="store_true",
                   help="generate rotor/casing STL files (default 0.5mm voxel)")
    o.add_argument("--step", action="store_true",
                   help="export STEP for re-CAD (optional `step` extra / "
                        "CadQuery): exact revolved solids + one sample "
                        "blade per row in machine coordinates")
    o.add_argument("--step-mode",
                   choices=["parts", "assembly", "blade0", "detailed",
                            "full"],
                   default="parts",
                   help="STEP scope: one file per physical part (+ "
                        "parts/README.txt with blade counts, default), "
                        "also a named cq.Assembly, a single rotor-1 blade, "
                        "the detailed hierarchical assembly (shaft, discs "
                        "with broached fir-tree slots, bladed rotors, "
                        "casing rings with flanges, tie bolts), or the "
                        "full machine with every row patterned")
    o.add_argument("--check-interference", action="store_true",
                   help="pairwise boolean check over the assembled machine: "
                        "no two parts may share material. Reports the "
                        "offending pair and the overlapped volume, and "
                        "fails the run if anything overlaps")
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


def parse_fixed(items) -> dict:
    """['phi1=0.6', 'n_stages=5'] → {'phi1': 0.6, 'n_stages': 5.0}.
    La validación de nombres la hace DesignSpec."""
    fixed = {}
    for it in items or []:
        if "=" not in it:
            raise ValueError(f"--fix espera VAR=VALOR, llegó '{it}'")
        name, _, raw = it.partition("=")
        try:
            fixed[name.strip()] = float(raw)
        except ValueError:
            raise ValueError(f"--fix {name}: valor no numérico '{raw}'")
    return fixed


def apply_machine_class(args) -> dict:
    """Parámetros de CLASE DE MÁQUINA: no son variables de diseño (el
    optimizador no los busca) sino decisiones del ingeniero que cambian
    la física del meanline. Se escriben en el módulo y se registran en
    run_meta para que la corrida sea reproducible.

    * `vortex_n`: un compresor clásico va a vórtice libre (−1); un FAN de
      HTR bajo NO puede — el vórtice libre le hunde la reacción de raíz.
    * `m_exit_max`: 0.55 es el límite amable para un difusor/cámara; un
      fan que descarga a un conducto de bypass corre 0.6-0.7.
    """
    import physics_core as _pc
    import geometry_generator as _gg
    meta = {}
    if args.vortex_n is not None:
        _pc.VORTEX_N = _gg.VORTEX_N = float(args.vortex_n)
        meta["vortex_n"] = float(args.vortex_n)
    if args.m_exit_max is not None:
        _pc.M_EXIT_MAX = float(args.m_exit_max)
        meta["m_exit_max"] = float(args.m_exit_max)
    meta.setdefault("vortex_n", _pc.VORTEX_N)
    meta.setdefault("m_exit_max", _pc.M_EXIT_MAX)
    return meta


def build_spec(args) -> DesignSpec:
    kw = dict(PR_target=args.pr, massflow=args.mdot, RPM_max=args.rpm_max,
              U_tip_max=args.utip_max, r_tip_max_mm=args.rtip_max,
              n_stages_max=args.nstages_max,
              power_max_W=args.power_max, T0_in=args.t0, P0_in=args.p0,
              material=args.material, fixed_vars=parse_fixed(args.fix))
    if args.sm_flow_min is not None:
        kw["sm_flow_min"] = float(args.sm_flow_min)
    if args.spec_file:
        with open(args.spec_file) as f:
            kw.update(json.load(f))
    return DesignSpec(**kw)


def warn_fixed_out_of_bounds(spec: DesignSpec):
    """Fijar fuera de los bounds del optimizador es legal (la física los
    tolera) pero merece aviso: el surrogate extrapola en esa dimensión."""
    for name, val in spec.fixed_vars.items():
        i = VAR_NAMES.index(name)
        if not (BOUNDS_LO[i] <= val <= BOUNDS_HI[i]):
            print(ui.warn(f"--fix {name}={val:g} fuera de los bounds "
                          f"[{BOUNDS_LO[i]:g}, {BOUNDS_HI[i]:g}] del "
                          "optimizador — el surrogate extrapolará"))


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


# ----------------------------------------------------------- deliverables
def emit_deliverables(args, spec, result, outdir, ckpt=None, logger=None,
                      title="Recommended design"):
    """Capa 5 compartida por los tres caminos del CLI (optimización,
    --eval-theta, --pareto-pick): estructural + geometría + reporte +
    dataset (si hay checkpoint) + STL/map opcionales + panel final."""
    logger = logger or StyledLog(quiet=args.quiet)
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
    geo_dir = os.path.join(outdir, "geometry")
    manifest = geometry_generator.generate(
        theta, rec, geo_dir, structural=structural,
        run_meta=dict(checkpoint=ckpt or "", rounds=args.rounds,
                      n_init=args.n_init, fidelity=args.fidelity))
    print("  " + ui.ok(f"geometry   → {geo_dir}/ ({len(manifest)} files)"))
    bcs = cfx_boundary_conditions(theta, rec)
    figs_dir = os.path.join(outdir, "figures")
    report = report_generator.generate_report(
        spec, result, os.path.join(outdir, "report.html"),
        geometry_manifest=manifest, cfx_bcs=bcs,
        figures_dir=figs_dir, include_figures=not args.no_figures,
        structural=structural, log=logger)
    print("  " + ui.ok(f"report     → {report}"))
    if not args.no_figures and os.path.isdir(figs_dir):
        n_figs = len([f for f in os.listdir(figs_dir) if f.endswith(".png")])
        if n_figs:
            print("  " + ui.ok(f"figures    → {figs_dir}/ ({n_figs} PNG)"))
    if ckpt and os.path.exists(ckpt):
        n_rows = export_dataset_csv(ckpt, os.path.join(outdir, "dataset.csv"))
        print("  " + ui.ok(f"dataset    → {n_rows} rows (flywheel)"))

    # ---- STEP de ensamble para re-CAD (opcional, CadQuery) ----
    if getattr(args, "step", False):
        contract_path = os.path.join(geo_dir, "axial_compressor.json")
        with open(contract_path, encoding="utf-8") as f:
            contract = json.load(f)
        paths = geometry_generator.export_step(contract, geo_dir,
                                               mode=args.step_mode)
        if paths:
            print("  " + ui.ok(f"step       → {len(paths)} files "
                               f"({args.step_mode}) in {geo_dir}/"))
        else:
            print(ui.warn("cadquery not installed — STEP export skipped "
                          "(pip install cadquery)"))

    # ---- Interferencias del ensamble (G-04) ----
    # Los fallos de geometría que han aparecido hasta ahora los encontró
    # una persona abriendo el STEP. Esto los busca solo: el volumen de la
    # intersección de dos piezas que no deben ocupar el mismo sitio tiene
    # que ser cero.
    if getattr(args, "check_interference", False):
        contract_path = os.path.join(geo_dir, "axial_compressor.json")
        with open(contract_path, encoding="utf-8") as f:
            contract = json.load(f)
        cq_mod = geometry_generator._cq()
        if cq_mod is None:
            print(ui.warn("cadquery not installed — interference check "
                          "skipped (pip install cadquery)"))
        else:
            parts = geometry_generator._asm_parts(
                geometry_generator._cq_detailed_machine(cq_mod, contract))
            n_sample = len(geometry_generator._interference_sample(parts))
            hits = geometry_generator.assembly_interferences(
                contract, parts=parts)
            txt = geometry_generator.format_interferences(hits, n_sample)
            if hits:
                args._interference_failed = True
                print(ui.err("interference  → " + txt.splitlines()[0]))
                for line in txt.splitlines()[1:]:
                    print("  " + ui.c(line, "err"))
            else:
                print("  " + ui.ok("interference → " + txt))

    # ---- Generación de STL con C# (opcional) ----
    if args.stl or args.voxel is not None:
        voxel_val = args.voxel if args.voxel is not None else 0.5
        contract_path = os.path.join(geo_dir, "axial_compressor.json")
        if os.path.exists(contract_path):
            generate_stls(contract_path, outdir, voxel_val)
        else:
            print(ui.warn(f"Could not find {contract_path} to generate STLs."))

    # ---- Mapa off-design opcional ----
    if args.map or args.map_vsv:
        vsv_mode = "auto" if args.map_vsv else "none"
        mp = compressor_map(theta, vsv=vsv_mode)
        map_path = os.path.join(outdir, "map.csv")
        with open(map_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["speed_frac", "rpm", "vsv_deg", "mdot_kgs", "PR",
                        "eta_poly", "min_SM", "M_rel_tip1",
                        "stage_stall_first", "stall", "choke",
                        "mdot_surge", "mdot_choke", "mdot_wl", "PR_wl",
                        "sm_flow"])
            for sl in mp["speedlines"]:
                wl = sl.get("working_point") or {}
                for q in sl["points"]:
                    w.writerow([sl["speed_frac"], sl["rpm"],
                                sl.get("vsv_deg", 0.0), q["mdot"],
                                q["PR"], q["eta_poly"], q["min_SM"],
                                q["M_rel_tip1"], q["stage_stall_first"],
                                int(q["stall"]), int(q["choke"]),
                                sl["mdot_surge"], sl["mdot_choke"],
                                wl.get("mdot"), wl.get("PR"),
                                wl.get("sm_flow")])
        print("  " + ui.ok(f"map        → {map_path}"))
        print(ui.c(f"\n  Off-design map (VSV {vsv_mode}) — surge…choke, "
                   "working point and margin:", "dim"))
        for sl in mp["speedlines"]:
            su, ch = sl["mdot_surge"], sl["mdot_choke"]
            wl = sl.get("working_point")
            rng_txt = (f"{su:.2f} – {ch:.2f} kg/s" if su is not None
                       else ui.c("no stable range", "warn"))
            if wl:
                pct = 100.0 * wl["sm_flow"]
                sty = ("ok" if wl["sm_flow"] >= pc_SM_MIN
                       else ("warn" if pct > 5.0 else "err"))
                wl_txt = ("  " + ui.c(f"SM {pct:.0f}%", sty)
                          + ui.c(f" @ {wl['mdot']:.2f} kg/s · PR "
                                 f"{wl['PR']:.2f} · 1st stall st."
                                 f"{wl['stage_stall_first'] + 1}", "dim"))
            else:
                wl_txt = ""
            vsv_txt = (ui.c(f"  VSV {sl['vsv_deg']:.0f}°", "accent")
                       if sl.get("vsv_deg") else "")
            print(f"    {ui.c('N=%.2f' % sl['speed_frac'], 'key')} "
                  f"({sl['rpm']:,.0f} RPM): {rng_txt}{vsv_txt}{wl_txt}")

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
        ("Surge margin", _sm_row(rec)),
        ("U_tip", f"{rec['U_tip']:.0f} m/s"),
        ("Status", status),
        ("Pareto", f"{len(result['pareto_front'])} points ({feas} feasible)"),
    ]
    if ckpt:
        rows.append(("Checkpoint", ui.c(ckpt, "dim")))
    if structural:
        s_ok = structural["feasible_struct"]
        rows.insert(5, ("Structure",
                        (ui.c(" MARGINS OK ", "ok", bold=True) if s_ok
                         else ui.c(" MARGIN VIOLATED ", "err", bold=True))
                        + ui.c(f"  {structural['material']} · "
                               f"burst {structural['burst_margin']:.2f} · "
                               f"AN² {structural['AN2_in2rpm2']:.1e} "
                               "in²rpm²", "dim")))
    print(ui.panel(title, rows, style="ok" if feasible else "warn"))
    print()
    print("  " + ui.c("Open report:", "dim") + " "
          + ui.c(report, "brand", bold=True))
    print()
    return report


def _exit_code(args) -> int:
    """3 si el chequeo de interferencias encontró piezas solapadas.

    Un entregable con interferencias no es un entregable: el CLI tiene que
    decirlo con el código de salida, no solo por pantalla, para que un CI
    o un script lo note.
    """
    return 3 if getattr(args, "_interference_failed", False) else 0


# ----------------------------------------------------------- pareto mode
def run_pareto_mode(args) -> int:
    """--list-pareto / --pareto-pick: opera sobre un checkpoint existente
    sin re-optimizar. El frente se re-verifica con la física (L0)."""
    ck = args.from_checkpoint or os.path.join(args.outdir, "phyac_run.json")
    if not os.path.exists(ck):
        print(ui.warn(f"checkpoint not found: {ck} (use --from-checkpoint)"))
        return 2
    spec, front = pareto_from_checkpoint(ck)
    front.sort(key=lambda q: (not q["feasible"], -q["eta_poly"]))
    print(ui.rule("Verified Pareto front", "brand"))
    print(ui.c(f"  {ck} — {len(front)} points\n", "dim"))
    print(ui.c(f"  {'#':>3}  {'feas':^4} {'stg':>3} {'PR':>7} "
               f"{'eta_poly':>8} {'Mu':>6} {'P [kW]':>9}  source", "key"))
    for i, q in enumerate(front):
        feas_txt = (ui.c("  OK ", "ok") if q["feasible"]
                    else ui.c("  -- ", "warn"))
        print(f"  {i:>3} {feas_txt} {q.get('n_stages') or 0:>3} "
              f"{q['PR']:7.3f} {q['eta_poly']:8.3f} {q['Mu']:6.3f} "
              f"{q['power_W'] / 1e3:9.1f}  {q.get('source', '')}")
    if args.pareto_pick is None:
        print(ui.c("\n  deliverables for a point:  --pareto-pick N "
                   "[--stl] [--map]", "dim"))
        return 0
    n = args.pareto_pick
    if not 0 <= n < len(front):
        print(ui.warn(f"--pareto-pick {n} out of range 0..{len(front) - 1}"))
        return 2
    fidelity = Fidelity.L1 if args.fidelity == "L1" else Fidelity.L0
    result = evaluate_design(spec, np.array(front[n]["theta"]),
                             fidelity=fidelity)
    outdir = os.path.join(args.outdir, f"pareto_{n:02d}")
    os.makedirs(outdir, exist_ok=True)
    print()
    print("  " + ui.ok(f"Pareto point #{n} → {outdir}/"))
    emit_deliverables(args, spec, result, outdir, ckpt=None,
                      title=f"Pareto point #{n}")
    return _exit_code(args)


# ----------------------------------------------------------- main
def _sm_row(rec: dict) -> str:
    """Margen de bombeo en la LÍNEA DE TRABAJO — el número que un ingeniero
    de compresores mira antes que ningún otro. Sin línea de trabajo la
    palabra "margen" no tiene denominador (Dixon & Hall §5.9)."""
    sm = rec.get("sm_flow")
    if sm is None:
        return ui.c("not evaluated", "dim")
    pct = 100.0 * float(sm)
    sty = "ok" if sm >= pc_SM_MIN else ("warn" if pct > 5.0 else "err")
    txt = ui.c(f"{pct:.1f}%", sty, bold=True) \
        + ui.c(f"  (min {100 * pc_SM_MIN:.0f}%)", "dim")
    srg = rec.get("surge") or {}
    if srg.get("ok"):
        txt += ui.c(f" · {srg['mdot_surge']:.1f} → {srg['mdot_wl']:.1f} → "
                    f"{srg['mdot_choke']:.1f} kg/s · 1st stall st."
                    f"{srg['stage_stall_first'] + 1}", "dim")
    return txt


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

    # ---- modos sobre checkpoint existente (sin optimizar) ----
    if args.list_pareto or args.pareto_pick is not None:
        return run_pareto_mode(args)

    if args.interactive:
        args = run_wizard(args)
    if args.quick:
        args.n_init, args.rounds, args.batch_size = 100, 2, 8
    os.makedirs(args.outdir, exist_ok=True)
    machine_class = apply_machine_class(args)
    try:
        spec = build_spec(args)
    except ValueError as e:
        print(ui.warn(str(e)))
        return 2
    warn_fixed_out_of_bounds(spec)
    if abs(machine_class["vortex_n"] + 1.0) > 1e-9:
        print(ui.c(f"  vortex exponent n = {machine_class['vortex_n']:+.2f} "
                   "(controlled vortex, not free)", "accent"))
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

    logger = StyledLog(quiet=args.quiet)

    # ---- evaluación directa de un θ dado (sin optimizar) ----
    if args.eval_theta:
        try:
            vals = [float(x) for x in re.split(r"[,\s]+",
                                               args.eval_theta.strip()) if x]
            # θ completo (13 legacy o 15): el punto de operación viaja en
            # las posiciones 10-12 y sobreescribe los flags del spec
            if len(vals) in (13, len(VAR_NAMES)):
                spec.T0_in, spec.P0_in, spec.massflow = (vals[10], vals[11],
                                                         vals[12])
            per_stage = None
            if args.per_stage:
                with open(args.per_stage, encoding="utf-8") as f:
                    per_stage = json.load(f)
            result = evaluate_design(spec, np.array(vals), fidelity=fidelity,
                                     per_stage=per_stage)
        except (ValueError, OSError) as e:
            print(ui.warn(str(e)))
            return 2
        print(ui.panel("Direct evaluation (no optimization)", [
            (name, f"{val:g}")
            for name, val in zip(VAR_NAMES, result["theta"])]))
        emit_deliverables(args, spec, result, args.outdir, ckpt=None,
                          logger=logger, title="Evaluated design")
        return _exit_code(args)

    # ---- optimización (nueva o reanudada) ----
    ckpt = os.path.join(args.outdir, "phyac_run.json")
    designer = None
    if args.resume is not None:
        resume_path = (ckpt if args.resume == "__auto__" else args.resume)
        if not os.path.exists(resume_path):
            print(ui.warn(f"checkpoint not found: {resume_path}"))
            return 2
        designer = AutonomousAxialDesigner.load(
            resume_path, fidelity=fidelity, n_workers=args.workers,
            log=logger, seed=args.seed)
        spec = designer.spec       # el dataset se evaluó bajo este spec
        print(ui.ok(f"resumed: {resume_path} · {len(designer.recs)} previous "
                    f"evaluations · seed {designer.seed}"))

    budget = (f"{args.n_init} init + {args.rounds}×{args.batch_size} "
              "acquisition" + ("  (quick)" if args.quick else ""))
    if designer is not None:
        budget = (f"{len(designer.recs)} resumed + {args.rounds}×"
                  f"{args.batch_size} acquisition")
    panel_rows = [
        ("Target PR", ui.c(f"{spec.PR_target:.2f}", "val", bold=True)),
        ("Mass flow", f"{spec.massflow:.2f} kg/s"),
        ("Max RPM", f"{spec.RPM_max:,.0f}"),
        ("Max U_tip", f"{spec.U_tip_max:.0f} m/s"),
        ("Max r_tip", f"{spec.r_tip_max_mm:.0f} mm"),
        ("Max stages", f"{spec.n_stages_max}"),
        ("Inlet", f"T0={spec.T0_in:.1f} K   P0={spec.P0_in:,.0f} Pa"),
        ("Budget", ui.c(budget, "accent")),
        ("Fidelity", ui.c(args.fidelity, "brand")),
    ]
    if spec.fixed_vars:
        panel_rows.append(("Fixed vars", ui.c(
            "  ".join(f"{k}={v:g}" for k, v in spec.fixed_vars.items()),
            "accent")))
    if args.seed is not None:
        panel_rows.append(("Seed", f"{args.seed}"))
    print(ui.panel("Specification", panel_rows))
    print()
    print(ui.rule("Autonomous optimization", "brand"))

    if designer is not None:
        result = designer.run(n_init=args.n_init, rounds=args.rounds,
                              batch_size=args.batch_size,
                              checkpoint_path=ckpt)
    else:
        result = design(spec, rounds=args.rounds, n_init=args.n_init,
                        batch_size=args.batch_size, fidelity=fidelity,
                        n_workers=args.workers, checkpoint_path=ckpt,
                        log=logger, seed=args.seed)

    emit_deliverables(args, spec, result, args.outdir, ckpt=ckpt,
                      logger=logger)
    return _exit_code(args)


if __name__ == "__main__":
    sys.exit(main())
