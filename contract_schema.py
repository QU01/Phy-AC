"""
QUASAR Phy-AC · contract_schema.py
==================================
El contrato `axial_compressor.json` es la ÚNICA frontera entre la capa 5a
(Python) y sus consumidores: la capa 5c (C#/PicoGK), la vía CadQuery y
cualquier integrador externo. Este módulo publica esa frontera como un
JSON Schema versionado y la comprueba.

Por qué existe
--------------
Hasta la fase 9 el contrato ganó campos nuevos —abeto, tirantes, exponente
de vórtice, sangrado, margen en caudal— conservando el identificador
`phyac-axial-1`. Un consumidor escrito contra la versión anterior no tenía
forma de saber que le faltaba información: leía lo que entendía y aplicaba
defaults al resto, en silencio. Con `phyac-axial-2` la versión MAYOR sube
con el contrato y el consumidor puede rechazar lo que no entiende.

Regla de versionado
-------------------
El entero final de `phyac-axial-N` es la versión MAYOR. Añadir un campo
OPCIONAL no la sube; volver un campo obligatorio, quitarlo, renombrarlo o
cambiar su significado o sus unidades, sí. Un consumidor lee la mayor y
para con un mensaje claro si no coincide con la que soporta.

Validador
---------
`validate()` implementa el subconjunto de JSON Schema draft 2020-12 que el
esquema usa: type, properties, required, additionalProperties, items,
minItems/maxItems, enum, const, minimum/maximum (y sus variantes
exclusivas), oneOf y $ref local a $defs. Es deliberadamente propio: el
proyecto depende solo de numpy y una comprobación que se salta sola
cuando falta una dependencia opcional no es una comprobación.
"""

from __future__ import annotations

import json
import math
import os
import sys

SCHEMA_ID = "phyac-axial-2"
SCHEMA_MAJOR = 2

SCHEMA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "schemas")


def schema_path(major: int = SCHEMA_MAJOR) -> str:
    """Ruta del JSON Schema publicado para una versión mayor."""
    return os.path.join(SCHEMA_DIR, f"phyac-axial-{major}.schema.json")


def load_schema(major: int = SCHEMA_MAJOR) -> dict:
    with open(schema_path(major), encoding="utf-8") as f:
        return json.load(f)


def schema_major(contract: dict) -> int | None:
    """Versión MAYOR declarada por un contrato, o None si no es uno.

    Es lo que un consumidor debe mirar ANTES de leer nada más.
    """
    sid = contract.get("schema")
    if not isinstance(sid, str) or not sid.startswith("phyac-axial-"):
        return None
    try:
        return int(sid.rsplit("-", 1)[1])
    except (IndexError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Validador (subconjunto de draft 2020-12)
# ---------------------------------------------------------------------------
_TYPES = {
    "object": dict,
    "array": list,
    "string": str,
    "boolean": bool,
    "null": type(None),
}


def _is_type(value, name: str) -> bool:
    if name == "integer":
        # bool es subclase de int en Python; un booleano NO es un entero.
        return isinstance(value, int) and not isinstance(value, bool)
    if name == "number":
        return (isinstance(value, (int, float))
                and not isinstance(value, bool)
                and math.isfinite(value))
    py = _TYPES.get(name)
    if py is None:
        return True
    if py is bool:
        return isinstance(value, bool)
    return isinstance(value, py)


def _resolve(ref: str, root: dict) -> dict:
    """Resuelve un $ref local del estilo `#/$defs/nombre`."""
    if not ref.startswith("#/"):
        raise ValueError(f"$ref no local no soportado: {ref}")
    node = root
    for part in ref[2:].split("/"):
        node = node[part.replace("~1", "/").replace("~0", "~")]
    return node


def _check(node, schema: dict, root: dict, path: str, errs: list[str],
           limit: int) -> None:
    if len(errs) >= limit:
        return

    if "$ref" in schema:
        _check(node, _resolve(schema["$ref"], root), root, path, errs, limit)
        return

    def bad(msg: str) -> None:
        errs.append(f"{path or '<raíz>'}: {msg}")

    # --- type -------------------------------------------------------------
    t = schema.get("type")
    if t is not None:
        names = [t] if isinstance(t, str) else list(t)
        if not any(_is_type(node, n) for n in names):
            bad(f"se esperaba {'|'.join(names)}, hay "
                f"{type(node).__name__}")
            return

    # --- const / enum -----------------------------------------------------
    if "const" in schema and node != schema["const"]:
        bad(f"se esperaba el valor constante {schema['const']!r}, "
            f"hay {node!r}")
    if "enum" in schema and node not in schema["enum"]:
        bad(f"{node!r} no está entre {schema['enum']}")

    # --- oneOf ------------------------------------------------------------
    if "oneOf" in schema:
        sub_errs = []
        n_ok = 0
        for i, sub in enumerate(schema["oneOf"]):
            e: list[str] = []
            _check(node, sub, root, path, e, limit)
            if e:
                sub_errs.append(f"[alternativa {i}] " + "; ".join(e[:2]))
            else:
                n_ok += 1
        if n_ok != 1:
            bad(f"debe cumplir exactamente una alternativa (cumple {n_ok}): "
                + " | ".join(sub_errs))

    # --- numéricos --------------------------------------------------------
    if isinstance(node, (int, float)) and not isinstance(node, bool):
        for key, op, txt in (("minimum", lambda a, b: a >= b, "≥"),
                             ("maximum", lambda a, b: a <= b, "≤"),
                             ("exclusiveMinimum", lambda a, b: a > b, ">"),
                             ("exclusiveMaximum", lambda a, b: a < b, "<")):
            if key in schema and not op(node, schema[key]):
                bad(f"{node!r} no cumple {txt} {schema[key]}")

    # --- objetos ----------------------------------------------------------
    if isinstance(node, dict):
        props = schema.get("properties", {})
        for req in schema.get("required", []):
            if req not in node:
                bad(f"falta la propiedad obligatoria «{req}»")
        extra = schema.get("additionalProperties", True)
        for k, v in node.items():
            if k in props:
                _check(v, props[k], root, f"{path}.{k}" if path else k,
                       errs, limit)
            elif extra is False:
                bad(f"propiedad no declarada «{k}»")
            elif isinstance(extra, dict):
                _check(v, extra, root, f"{path}.{k}" if path else k,
                       errs, limit)
            if len(errs) >= limit:
                return

    # --- arrays -----------------------------------------------------------
    if isinstance(node, list):
        if "minItems" in schema and len(node) < schema["minItems"]:
            bad(f"{len(node)} elementos, mínimo {schema['minItems']}")
        if "maxItems" in schema and len(node) > schema["maxItems"]:
            bad(f"{len(node)} elementos, máximo {schema['maxItems']}")
        item_schema = schema.get("items")
        if isinstance(item_schema, dict):
            for i, v in enumerate(node):
                _check(v, item_schema, root, f"{path}[{i}]", errs, limit)
                if len(errs) >= limit:
                    return


def validate(contract: dict, schema: dict | None = None,
             limit: int = 40) -> list[str]:
    """Devuelve la lista de incumplimientos; vacía si el contrato es válido.

    `limit` corta el informe: un contrato roto de raíz genera miles de
    errores idénticos y ninguno de ellos añade información.
    """
    root = schema if schema is not None else load_schema()
    errs: list[str] = []
    _check(contract, root, root, "", errs, limit)
    return errs


# ---------------------------------------------------------------------------
# CLI: python contract_schema.py <ruta al axial_compressor.json>
# ---------------------------------------------------------------------------
def _main(argv: list[str]) -> int:
    # En una consola Windows con cp1252 el «✓» revienta el print. El
    # resto del proyecto hace lo mismo (ver test_phyac.py).
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    if len(argv) != 2:
        print(f"uso: python {os.path.basename(__file__)} "
              f"<axial_compressor.json>")
        return 2
    with open(argv[1], encoding="utf-8") as f:
        contract = json.load(f)

    major = schema_major(contract)
    if major is None:
        print(f"✗ «{argv[1]}» no declara un esquema phyac-axial-N")
        return 1
    if major != SCHEMA_MAJOR:
        print(f"✗ contrato en esquema mayor {major}; este validador "
              f"soporta {SCHEMA_MAJOR}. Regenera la geometría con la "
              f"versión actual de Phy-AC.")
        return 1

    errs = validate(contract)
    if errs:
        print(f"✗ {len(errs)} incumplimiento(s) de {SCHEMA_ID}:")
        for e in errs:
            print(f"  · {e}")
        return 1
    print(f"✓ contrato válido contra {SCHEMA_ID}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv))
