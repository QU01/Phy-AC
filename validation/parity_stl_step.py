#!/usr/bin/env python3
"""
QUASAR Phy-AC · validation/parity_stl_step.py
=============================================
Comprueba que la CAPA 5c (C#/PicoGK, ruta de fabricacion) y la VIA
CADQUERY (ruta de re-CAD) describen LA MISMA MAQUINA.

Por que existe
--------------
Las dos rutas leen el mismo contrato y construyen por su cuenta. Entre
las fases 8.2 y 10 se separaron sin que nada avisara: el C# seguia
hundiendo la raiz del alabe en el cuerpo (`RootSinkMm`), no conocia el
abeto, no recortaba contra la vena y no tenia tornilleria de disco. Un
usuario que comparase el STL con el STEP encontraba la discrepancia en el
primer vistazo; la suite no. Esto lo comprueba.

Dos niveles, porque cuestan cosas muy distintas:

  perfil   matematica pura (perfil de abeto, longitud de raiz y pendiente
           de plataforma por fila). Sin PicoGK, sin GPU, milisegundos.
           Es lo que caza que los dos puertos se separen.
  piezas   construye cada pieza en las dos rutas y compara VOLUMEN y
           radio exterior. Necesita el runtime de PicoGK y CadQuery, y
           tarda minutos. Es lo que caza que describan maquinas
           distintas.

La descomposicion en piezas NO es la misma a proposito: el STL parte la
maquina en piezas IMPRIMIBLES (un blisk por etapa) y el STEP en piezas de
re-CAD (disco, alabes y casco por separado). La paridad se compara por
GRUPOS equivalentes, declarados en `PART_GROUPS`.

Uso:
    python validation/parity_stl_step.py <axial_compressor.json> [voxel]
"""

from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import geometry_generator as gg     # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DLL = os.path.join(ROOT, "AxialCompressorDesigner.Example", "bin", "Debug",
                   "net9.0-windows", "AxialCompressorDesigner.Example.dll")

# Tolerancias DECLARADAS, fijadas contra medida (2026-08-16, maquina de 2
# etapas a voxel 0.6 mm): Shaft -0.1%, Rotor -1.1%, Casing -0.7% en
# volumen; radios exteriores dentro de 0.4 mm. El margen de 5% cubre la
# discretizacion en voxeles y las diferencias de muestreo entre las
# revoluciones de las dos rutas (64 bandas en PicoGK, 40 puntos en
# CadQuery). El radio exterior se aprieta mas: es una cota, no un
# volumen, y es lo que decide si la pieza cabe.
TOL_VOLUME_FRAC = 0.05      # 5% en volumen por conjunto
TOL_RADIUS_MM = 1.5         # 1.5 mm en radio exterior (~2-3 voxeles)

# El STEP no modela los filetes de raiz y el STL si, asi que la
# comparacion se hace con el filete APAGADO en la ruta STL. No es un
# apano: es comparar lo mismo con lo mismo. Y de paso dejo anotado lo que
# midio esta comparacion — el paso de filete de la capa 5c anade 117 cm3
# a un anillo de carcasa de 488 cm3, un 24%, cuando un filete de 2 mm
# sobre las raices de ~70 vanos no puede pasar de unos pocos cm3. Ese
# `Fillet` (over-offset morfologico) esta anadiendo mucho mas material
# del que dice; pendiente de investigar aparte.
PARITY_DISABLES_FILLET = True


def dotnet_bin() -> str | None:
    exe = shutil.which("dotnet")
    if exe:
        return exe
    std = r"C:\Program Files\dotnet\dotnet.exe"
    return std if os.path.exists(std) else None


def csharp_available() -> tuple[bool, str]:
    """(disponible, motivo). No compila: eso lo hace el que llama."""
    if dotnet_bin() is None:
        return False, "no se encontro 'dotnet'"
    if not os.path.exists(DLL):
        return False, f"no existe {os.path.relpath(DLL, ROOT)} (dotnet build)"
    return True, ""


def csharp_profile(contract_path: str) -> dict:
    """Vuelca el perfil de abeto y la geometria de raiz de la capa 5c.

    Matematica pura del lado C#: no arranca PicoGK ni necesita GPU.
    """
    out = subprocess.run([dotnet_bin(), DLL, contract_path, "--profile-json"],
                         capture_output=True, text=True, cwd=ROOT)
    if out.returncode != 0:
        raise RuntimeError(f"--profile-json fallo: {out.stderr[-2000:]}")
    return json.loads(out.stdout)


def python_profile(contract: dict) -> dict:
    """El MISMO volcado desde la via CadQuery, para comparar."""
    ft = contract["assembly"]["firtree"]
    asm = contract["assembly"]
    rows = []
    for st in contract["stages"]:
        rot = st["rotor"]
        ax, slope = gg._root_axial_and_slope(contract, rot)
        rows.append(dict(
            stage=st["index"],
            z_center_mm=float(rot["z_center_mm"]),
            r_hub_c_mm=gg._r_on_line(contract["annulus"]["hub"],
                                     rot["z_center_mm"]),
            root_axial_mm=float(ax),
            platform_slope=float(slope),
            n_blades=int(rot["n_blades"])))
    args = (ft["w_top_mm"], ft["depth_mm"], ft["n_teeth"], ft["taper_deg"],
            ft["lobe_frac"])
    return dict(
        firtree_enabled=bool(ft["enabled"]),
        n_sections=gg.N_FIRTREE_SECTIONS,
        slot_axial_factor=1.05,
        tie_bolt_count=int(asm["tie_bolt_count"]),
        tie_bolt_d_mm=float(asm["tie_bolt_d_mm"]),
        rim_relief_mm=float(asm.get("rim_relief_mm", 0.10)),
        root_profile=[list(p) for p in gg.firtree_profile(*args, 0.0)],
        slot_profile=[list(p) for p in
                      gg.firtree_profile(*args, ft["clearance_mm"])],
        rows=rows)


def compare_profiles(cs: dict, py: dict, tol_mm: float = 1e-3) -> list[str]:
    """Diferencias entre los dos volcados. Vacia = las dos vias comparten
    la misma retencion de alabe."""
    errs: list[str] = []
    for k in ("firtree_enabled", "n_sections", "tie_bolt_count"):
        if cs[k] != py[k]:
            errs.append(f"{k}: C# {cs[k]} vs Python {py[k]}")
    for k in ("slot_axial_factor", "tie_bolt_d_mm", "rim_relief_mm"):
        if abs(float(cs[k]) - float(py[k])) > 1e-6:
            errs.append(f"{k}: C# {cs[k]} vs Python {py[k]}")

    for name in ("root_profile", "slot_profile"):
        a, b = cs[name], py[name]
        if len(a) != len(b):
            errs.append(f"{name}: {len(a)} puntos en C# vs {len(b)} en Python")
            continue
        worst = max(max(abs(pa[0] - pb[0]), abs(pa[1] - pb[1]))
                    for pa, pb in zip(a, b))
        if worst > tol_mm:
            errs.append(f"{name}: se separan {worst:.2e} mm "
                        f"(tolerancia {tol_mm:g} mm)")

    if len(cs["rows"]) != len(py["rows"]):
        errs.append(f"filas de rotor: {len(cs['rows'])} vs {len(py['rows'])}")
        return errs
    for rc, rp in zip(cs["rows"], py["rows"]):
        for k in ("z_center_mm", "r_hub_c_mm", "root_axial_mm"):
            if abs(float(rc[k]) - float(rp[k])) > tol_mm:
                errs.append(f"etapa {rp['stage']} {k}: "
                            f"{rc[k]} vs {rp[k]}")
        if abs(float(rc["platform_slope"])
               - float(rp["platform_slope"])) > 1e-5:
            errs.append(f"etapa {rp['stage']} platform_slope: "
                        f"{rc['platform_slope']} vs {rp['platform_slope']}")
        if int(rc["n_blades"]) != int(rp["n_blades"]):
            errs.append(f"etapa {rp['stage']} n_blades: "
                        f"{rc['n_blades']} vs {rp['n_blades']}")
    return errs


# ---------------------------------------------------------------------------
# Paridad por PIEZAS (cara: construye las dos rutas)
# ---------------------------------------------------------------------------
def csharp_parts(contract_path: str, outdir: str,
                 voxel: float = 1.0) -> dict:
    """Volumen y radio exterior por pieza de la capa 5c."""
    os.makedirs(outdir, exist_ok=True)
    rep = os.path.join(outdir, "parity_5c.json")
    if PARITY_DISABLES_FILLET:
        with open(contract_path, encoding="utf-8") as f:
            c = json.load(f)
        c["assembly"]["blade_fillet_r_mm"] = 0.0
        contract_path = os.path.join(outdir, "contract_nofillet.json")
        with open(contract_path, "w", encoding="utf-8") as f:
            json.dump(c, f)
    out = subprocess.run(
        [dotnet_bin(), DLL, contract_path, outdir, str(voxel),
         "--parity-json", rep],
        capture_output=True, text=True, cwd=ROOT)
    if out.returncode != 0 or not os.path.exists(rep):
        raise RuntimeError(f"--parity-json fallo: "
                           f"{(out.stderr or out.stdout)[-2000:]}")
    with open(rep, encoding="utf-8") as f:
        return json.load(f)["parts"]


# Las dos rutas parten la maquina en PLANOS DISTINTOS a proposito: el STL
# corta donde van las juntas apernadas reales (un blisk por etapa) y el
# STEP donde estan las piezas de re-CAD (disco, alabes y casco del tambor
# por separado, con los cascos ocupando solo los huecos entre bandas de
# disco). Comparar pieza a pieza compararia esa decision, no la maquina.
# Se comparan los tres CONJUNTOS que si tienen que coincidir.
GROUPS = {
    "Shaft": (["Shaft"], ["Shaft"]),
    "Rotor": (["RotorStage"], ["Disc", "Rotor", "DrumSpacer"]),
    "Casing": (["StatorRing"], ["CasingRing", "Stator", "IGV_", "OGV_"]),
}


def _sum_group(items, prefixes) -> tuple[float, float]:
    vol = 0.0
    r_max = 0.0
    for name, vol_i, r_i in items:
        if not any(name.startswith(p) for p in prefixes):
            continue
        vol += vol_i
        r_max = max(r_max, r_i)
    return vol, r_max


def group_csharp(parts: dict) -> dict:
    items = [(n, float(d["volume_mm3"]), float(d["r_max_mm"]))
             for n, d in parts.items()]
    out = {}
    for g, (stl_pref, _) in GROUPS.items():
        vol, r = _sum_group(items, stl_pref)
        if vol > 0.0:
            out[g] = dict(volume_mm3=vol, r_max_mm=r)
    return out


def cadquery_groups(contract: dict, k_span: int = 2,
                    k_pts: int = 2) -> dict:
    """Los MISMOS tres conjuntos desde la via CadQuery.

    La tornilleria (TieBolt / CasingBolt) queda fuera a proposito: el STL
    lleva los taladros pero no los pernos, que no se imprimen.
    """
    cq = gg._cq()
    if cq is None:
        raise RuntimeError("cadquery no instalado")
    parts = gg._asm_parts(gg._cq_detailed_machine(cq, contract, k_span,
                                                  k_pts))
    items = []
    for name, shp in parts:
        leaf = name.rsplit("/", 1)[-1]
        try:
            bb = shp.BoundingBox()
            items.append((leaf, float(shp.Volume()),
                          max(abs(bb.xmin), abs(bb.xmax),
                              abs(bb.ymin), abs(bb.ymax))))
        except Exception:
            continue
    out = {}
    for g, (_, step_pref) in GROUPS.items():
        vol, r = _sum_group(items, step_pref)
        if vol > 0.0:
            out[g] = dict(volume_mm3=vol, r_max_mm=r)
    return out


def compare_parts(cs: dict, cqp: dict,
                  tol_vol: float = TOL_VOLUME_FRAC,
                  tol_r: float = TOL_RADIUS_MM) -> list[str]:
    errs: list[str] = []
    common = sorted(set(cs) & set(cqp))
    if not common:
        return ["ningun grupo de piezas en comun entre las dos rutas"]
    for g in common:
        a, b = cs[g], cqp[g]
        va, vb = float(a["volume_mm3"]), float(b["volume_mm3"])
        d = abs(va - vb) / max(vb, 1e-9)
        if d > tol_vol:
            errs.append(f"{g}: volumen {va / 1000:.1f} vs {vb / 1000:.1f} cm3 "
                        f"({100 * (va - vb) / max(vb, 1e-9):+.1f}%, "
                        f"tolerancia {100 * tol_vol:.0f}%)")
        dr = abs(float(a["r_max_mm"]) - float(b["r_max_mm"]))
        if dr > tol_r:
            errs.append(f"{g}: radio exterior {a['r_max_mm']:.2f} vs "
                        f"{b['r_max_mm']:.2f} mm ({dr:.2f} mm, tolerancia "
                        f"{tol_r:g} mm)")
    missing = sorted((set(cs) | set(cqp)) - set(common))
    if missing:
        errs.append("grupos sin pareja: " + ", ".join(missing))
    return errs


def _main(argv: list[str]) -> int:
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    if len(argv) < 2:
        print(f"uso: python {os.path.basename(__file__)} "
              f"<axial_compressor.json> [voxel_mm] [--profile-only]")
        print("  --profile-only: solo la parte barata (perfil de abeto y "
              "geometria de raiz); no necesita GPU ni CadQuery.")
        return 2
    contract_path = os.path.abspath(argv[1])
    profile_only = "--profile-only" in argv
    rest = [a for a in argv[2:] if not a.startswith("--")]
    voxel = float(rest[0]) if rest else 1.0
    with open(contract_path, encoding="utf-8") as f:
        contract = json.load(f)

    ok, why = csharp_available()
    if not ok:
        print(f"x capa 5c no disponible: {why}")
        return 2

    print("— perfil de abeto y geometria de raiz")
    errs = compare_profiles(csharp_profile(contract_path),
                            python_profile(contract))
    if errs:
        for e in errs:
            print(f"  x {e}")
        return 1
    print("  ok las dos vias comparten la retencion de alabe")
    if profile_only:
        return 0

    print(f"— piezas (voxel {voxel} mm"
          + (", filete de raiz apagado en las dos rutas"
             if PARITY_DISABLES_FILLET else "") + ")")
    outdir = os.path.join(os.path.dirname(contract_path), "_parity")
    cached = os.path.join(outdir, "parity_5c.json")
    if os.environ.get("PHYAC_PARITY_REUSE") and os.path.exists(cached):
        with open(cached, encoding="utf-8") as f:
            cs = group_csharp(json.load(f)["parts"])
        print(f"  (reusando {os.path.relpath(cached, ROOT)})")
    else:
        cs = group_csharp(csharp_parts(contract_path, outdir, voxel))
    cqp = cadquery_groups(contract)
    errs = compare_parts(cs, cqp)
    for g in sorted(set(cs) & set(cqp)):
        va = float(cs[g]["volume_mm3"]) / 1000.0
        vb = float(cqp[g]["volume_mm3"]) / 1000.0
        print(f"  {g:<16} STL {va:9.1f} cm3   STEP {vb:9.1f} cm3   "
              f"{100 * (va - vb) / max(vb, 1e-9):+6.1f}%   "
              f"r {cs[g]['r_max_mm']:7.2f} / {cqp[g]['r_max_mm']:7.2f} mm")
    if errs:
        for e in errs:
            print(f"  x {e}")
        return 1
    print("  ok las dos rutas describen la misma maquina")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv))
