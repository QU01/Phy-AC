"""QUASAR Phy-AC — Capa 0: pipeline de datos públicos agregados.

Política Phy-* (heredada de Phy-HL/Phy-Wing): SOLO agregados pequeños con
whitelist hard-coded y manifiesto SHA-256 — nunca tiers masivos de campos.

Para compresores axiales NO existe un dataset público masivo (DATED en
Zenodo es de centrífugos). Los agregados públicos relevantes son tablas
pequeñas de informes NASA (dominio público):

  - Punto de diseño y specs de NASA Stage 35 / Rotor 37 / Rotor 67 / E³ HPC
    → viven hard-coded con cita en validation/machines.py (patrón Phy-CC).
  - Anclas de correlación de cascadas (Lieblein θ/c vs D_eq, desviación de
    Carter) → data/cascade/*.csv versionados EN el repo con procedencia en
    cabecera; se regeneran con --build y se verifican contra el manifiesto.

REMOTE_FILES está VACÍO deliberadamente: las URL de turbmodels.larc.nasa.gov
rotaron (verificado 2026-07: redirigen a nasa.gov), así que todo agregado se
versiona en el repo y la validación funciona offline. La maquinaria de
descarga+manifiesto se conserva por si aparece un mirror estable.

El dataset de ENTRENAMIENTO del surrogate no se descarga: se autogenera
muestreando el meanline (flywheel dataset.csv de cada corrida — ver
neural_optimizer.py), la misma estrategia con la que se construyó DATED.

Uso CLI:
  python data_pipeline.py                  # build + validate
  python data_pipeline.py --build          # regenera las anclas de cascada
  python data_pipeline.py --validate       # verifica manifiesto + tablas
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# ----------------------------------------------------------------------------
# Whitelist remota (vacía — ver docstring) y registro de tablas en repo
# ----------------------------------------------------------------------------

REMOTE_FILES: dict[str, str] = {}

#: Tablas agregadas versionadas en el repo. Cada CSV lleva su procedencia
#: en las líneas de comentario (#) de cabecera.
IN_REPO_FILES = (
    "cascade/lieblein_theta_c.csv",     # θ/c vs D_eq (Lieblein 1959)
    "cascade/carter_deviation.csv",     # m de Carter vs stagger (regla clásica)
    "cascade/howell_wdf.csv",           # work-done factor vs nº de etapa
)

DEFAULT_DATA_DIR = Path(__file__).resolve().parent / "data"


# ----------------------------------------------------------------------------
# Manifiesto SHA-256
# ----------------------------------------------------------------------------

def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, dest: Path, retries: int = 3, timeout: int = 60) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": "phy-ac/0.1"})
    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                dest.write_bytes(r.read())
            return
        except Exception as exc:
            last_err = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"No se pudo descargar {url}: {last_err}")


def fetch(data_dir: Path = DEFAULT_DATA_DIR, force: bool = False) -> dict:
    """Descarga los archivos de la whitelist remota (hoy: ninguno).

    Idempotente; devuelve/actualiza el manifiesto de descargas.
    """
    raw = data_dir / "raw"
    manifest_path = raw / "fetch_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    if not REMOTE_FILES:
        return manifest
    raw.mkdir(parents=True, exist_ok=True)
    for name, url in REMOTE_FILES.items():
        dest = raw / name
        if dest.exists() and not force:
            continue
        print(f"[fetch] {name} <- {url}")
        _download(url, dest)
        manifest[name] = {
            "url": url, "sha256": _sha256(dest), "bytes": dest.stat().st_size,
            "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        }
    manifest_path.write_text(json.dumps(manifest, indent=2))
    return manifest


def write_manifest(data_dir: Path = DEFAULT_DATA_DIR) -> Path:
    """Escribe data/manifest.json con el SHA-256 de las tablas en repo."""
    man = {}
    for rel in IN_REPO_FILES:
        p = data_dir / rel
        if p.exists():
            man[rel] = {"sha256": _sha256(p), "bytes": p.stat().st_size}
    out = data_dir / "manifest.json"
    out.write_text(json.dumps(man, indent=2))
    return out


def verify_manifest(data_dir: Path = DEFAULT_DATA_DIR) -> list[str]:
    """Devuelve la lista de problemas (vacía si todo verifica)."""
    problems = []
    man_path = data_dir / "manifest.json"
    if not man_path.exists():
        return [f"falta {man_path}"]
    man = json.loads(man_path.read_text())
    for rel in IN_REPO_FILES:
        p = data_dir / rel
        if not p.exists():
            problems.append(f"falta {rel}")
        elif rel not in man:
            problems.append(f"{rel} sin entrada en el manifiesto")
        elif _sha256(p) != man[rel]["sha256"]:
            problems.append(f"{rel}: SHA-256 no coincide con el manifiesto")
    return problems


# ----------------------------------------------------------------------------
# Build: anclas de correlación (tabuladas desde las correlaciones publicadas)
# ----------------------------------------------------------------------------

def _write_csv(path: Path, header_comments: list[str], cols: list[str],
               rows: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # LF determinista (csv.writer emite \r\n por defecto en Windows):
    # los hashes del manifest deben ser idénticos en cualquier plataforma
    # y coincidir con la normalización eol=lf de .gitattributes.
    with open(path, "w", newline="", encoding="utf-8") as f:
        for line in header_comments:
            f.write(f"# {line}\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(cols)
        for row in rows:
            w.writerow([f"{v:.6g}" for v in row])


def build(data_dir: Path = DEFAULT_DATA_DIR) -> None:
    """Regenera las anclas de cascada desde las correlaciones publicadas.

    Son ANCLAS DE REGRESIÓN: fijan la forma de las correlaciones que
    physics_core/blade_profiles implementan, para que un cambio accidental
    de constantes dispare la validación (mismo rol que las anclas de
    Phy-CC/validate.py). No son mediciones digitalizadas.
    """
    import physics_core as pc

    deq = np.linspace(1.0, 2.2, 25)
    theta_c = np.array([pc._lieblein_theta_c(d) for d in deq])
    _write_csv(
        data_dir / "cascade" / "lieblein_theta_c.csv",
        ["Ancla de regresión: espesor de momento de estela vs difusión "
         "equivalente.",
         "Correlación de Lieblein (1959), fit clásico "
         "theta/c = 0.0045/(1 - 0.95·ln(Deq)),",
         "según NASA SP-36 cap. VII y Aungier, 'Axial-Flow Compressors' "
         "(2003).",
         "Generado por data_pipeline.build(); verificado por validate.py."],
        ["Deq", "theta_over_c"], np.column_stack([deq, theta_c]))

    # Regla de Carter: m = 0.23·(2a/c)² + ζ/500 (arco circular, a/c=0.5 →
    # m = 0.23 + stagger/500), Carter (1950) vía Dixon & Hall §3.
    stagger = np.linspace(0.0, 60.0, 13)
    m = 0.23 + stagger / 500.0
    _write_csv(
        data_dir / "cascade" / "carter_deviation.csv",
        ["Ancla de regresión: coeficiente m de la regla de Carter vs "
         "stagger (deg).",
         "m = 0.23·(2a/c)² + stagger/500 con línea de comba de arco "
         "circular (a/c = 0.5).",
         "Carter (1950); Dixon & Hall, 'Fluid Mechanics and Thermodynamics "
         "of Turbomachinery' §3.",
         "Generado por data_pipeline.build(); verificado por validate.py."],
        ["stagger_deg", "m_carter"], np.column_stack([stagger, m]))

    n_stage = np.arange(1, len(pc._HOWELL_WDF) + 1, dtype=float)
    _write_csv(
        data_dir / "cascade" / "howell_wdf.csv",
        ["Ancla de regresión: work-done factor medio vs número de etapa.",
         "Howell & Bonham (1950), tabla clásica reproducida en Dixon & "
         "Hall §5 y Saravanamuttoo et al.",
         "Generado por data_pipeline.build(); verificado por validate.py."],
        ["stage", "work_done_factor"],
        np.column_stack([n_stage, np.array(pc._HOWELL_WDF)]))

    write_manifest(data_dir)
    print(f"[build] {len(IN_REPO_FILES)} anclas regeneradas + manifest.json")


# ----------------------------------------------------------------------------
# Acceso
# ----------------------------------------------------------------------------

def load_table(rel: str, data_dir: Path = DEFAULT_DATA_DIR
               ) -> tuple[list[str], np.ndarray]:
    """Lee un CSV de agregados: (columnas, matriz float)."""
    p = data_dir / rel
    with open(p, encoding="utf-8") as f:
        rows = [r for r in csv.reader(
            (ln for ln in f if not ln.startswith("#")))]
    cols, data = rows[0], np.array(rows[1:], dtype=float)
    return cols, data


# ----------------------------------------------------------------------------
# Validación
# ----------------------------------------------------------------------------

@dataclass
class Check:
    name: str
    ok: bool
    detail: str = ""


def validate(data_dir: Path = DEFAULT_DATA_DIR) -> list[Check]:
    checks: list[Check] = []

    def add(name: str, ok: bool, detail: str = "") -> None:
        checks.append(Check(name, bool(ok), detail))

    problems = verify_manifest(data_dir)
    add("manifiesto SHA-256 verifica", not problems, "; ".join(problems))
    if problems:
        return checks

    cols, t = load_table("cascade/lieblein_theta_c.csv", data_dir)
    add("lieblein: θ/c monótona creciente y en rango",
        bool(np.all(np.diff(t[:, 1]) >= 0) and t[0, 1] > 0.003
             and t[-1, 1] < 0.10),
        f"θ/c ∈ [{t[:,1].min():.4f}, {t[:,1].max():.4f}]")

    cols, c = load_table("cascade/carter_deviation.csv", data_dir)
    add("carter: m ∈ [0.20, 0.40] y creciente con stagger",
        bool(np.all(np.diff(c[:, 1]) > 0) and c[0, 1] >= 0.20
             and c[-1, 1] <= 0.40),
        f"m ∈ [{c[:,1].min():.3f}, {c[:,1].max():.3f}]")

    cols, w = load_table("cascade/howell_wdf.csv", data_dir)
    add("howell: wdf decreciente y ∈ [0.80, 1.00]",
        bool(np.all(np.diff(w[:, 1]) < 0) and w[:, 1].min() >= 0.80
             and w[:, 1].max() <= 1.00),
        f"wdf ∈ [{w[:,1].min():.3f}, {w[:,1].max():.3f}]")
    return checks


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    args = ap.parse_args(argv)

    do_all = not (args.fetch or args.build or args.validate)
    if args.fetch:
        fetch(args.data_dir, force=args.force)
    if args.build or do_all:
        build(args.data_dir)
    if args.validate or do_all:
        checks = validate(args.data_dir)
        n_fail = sum(not c.ok for c in checks)
        for c in checks:
            mark = "OK " if c.ok else "FAIL"
            print(f"  [{mark}] {c.name}" + (f"  ({c.detail})" if c.detail else ""))
        print(f"[validate] {len(checks) - n_fail}/{len(checks)} checks OK")
        return 1 if n_fail else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
