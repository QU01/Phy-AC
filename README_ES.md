# QUASAR Phy-AC — Diseño autónomo de compresores axiales

<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-Apache%202.0-blue.svg" alt="Licencia: Apache 2.0"></a>
  <img src="https://img.shields.io/badge/python-%E2%89%A53.10-blue.svg" alt="Python >= 3.10">
  <img src="https://img.shields.io/badge/.NET-9-blue.svg" alt=".NET 9">
</p>

<p align="center">
  <img src="assets/gallery/assembly_meridional.png" alt="Corte meridional del rotor y la carcasa de QUASAR Phy-AC, renderizado desde los STL multi-parte verificados" width="820">
</p>
<p align="center"><i>Corte meridional de las piezas imprimibles generadas por la capa 5c (C#/PicoGK) desde un diseño verificado de QUASAR Phy-AC — eje, discos-álabe por etapa y casco del hub (izquierda); anillos de carcasa con vanos de estátor y bridas apernadas (derecha).</i></p>

> **Espec → geometría verificada → piezas imprimibles.** El ingeniero
> escribe seis números (PR objetivo, flujo másico, RPM máx, velocidad de
> punta máx, radio de punta máx, condiciones de entrada) y recibe, en el
> CPU de una laptop y sin proponer geometría alguna: el frente de Pareto
> factible físicamente verificado, el diseño multietapa recomendado, el
> desglose de pérdidas y carga por etapa, las condiciones de frontera de
> CFD, la lista de materiales, **y los STL imprimibles de cada pieza —
> eje, discos-álabe por etapa, y anillos de carcasa con bridas apernadas
> y puerto de sangrado.**
>
> Documentos: [docs/Quasar_PhyAC_Science.md](docs/Quasar_PhyAC_Science.md)
> (ecuaciones, pseudocódigo, referencias) ·
> [docs/VALIDATION.md](docs/VALIDATION.md) (validación física contra
> máquinas NASA). Phy-AC es el hermano axial de **Phy-CC** (compresores
> centrífugos) y comparte su arquitectura: *prior físico calibrado →
> ensemble residual con puerta de incertidumbre → NSGA-II restringido →
> informe autocontenido*.

## Visualización

El `report.html` de cada corrida embebe figuras matplotlib generadas por
[`visualization.py`](visualization.py) (capa 5b, dependencia opcional)
desde el mismo record físico verificado del que salen las tablas: mapa
fuera de diseño, carga por etapa, annulus meridional y secciones de álabe.

<table>
<tr>
<td width="50%"><img src="assets/gallery/annulus.png" alt="Annulus meridional"></td>
<td width="50%"><img src="assets/gallery/sections.png" alt="Secciones del rotor 1"></td>
</tr>
<tr>
<td align="center"><sub><b>Annulus meridional</b> — línea media constante; rotores en azul, estátores en verde; la altura de álabe cae con la compresión</sub></td>
<td align="center"><sub><b>Secciones del rotor 1</b> (hub/media/punta, staggered) — twist de free vortex: el stagger crece y la comba cae hacia la punta; DCA sobre M≈0.8, NACA-65 debajo</sub></td>
</tr>
<tr>
<td width="50%"><img src="assets/gallery/stage_loading.png" alt="Carga por etapa"></td>
<td width="50%"><img src="assets/gallery/map.png" alt="Mapa fuera de diseño"></td>
</tr>
<tr>
<td align="center"><sub><b>Carga por etapa</b> — ψ, DF de Lieblein y margen de bombeo de Koch; línea roja = SM mínimo exigido por g</sub></td>
<td align="center"><sub><b>Mapa fuera de diseño</b> — speedlines PR(ṁ) con línea de bombeo y choke (proxy de tendencia L0; geometría fija o schedule VSV auto con --map-vsv)</sub></td>
</tr>
</table>

## Pipeline

```
espec (6 números)
   │   phyac_cli.py · neural_optimizer.design()          [Python, capas 1–4]
   ▼
diseño óptimo verificado (vector θ + record físico)
   │   geometry_generator.py                             [capa 5a]
   ▼
geometry/axial_compressor.json · bom.csv · BCs CFX ────→ report.html [capa 5b]
   │   PhyACImport.FromContractJson(...) + PicoGK        [capa 5c, C#/.NET 9]
   ▼
STLs Shaft + RotorStage{i} + StatorRing{i}  (cada pieza, imprimible)
```

Las capas 1–5b son Python (núcleo: solo NumPy). La capa 5c es la librería
C# **AxialCompressorDesigner** (.NET 9 + paquete NuGet
[PicoGK](https://picogk.org) 2.2.0), que convierte el diseño optimizado en
las piezas reales de la máquina — construcción de turborreactor, partida
donde van las juntas apernadas reales:

| Pieza | Archivo exportado | Contenido |
|---|---|---|
| **Eje motriz** | `<Name>_Shaft.stl` | Eje con muñones de montaje en ambos extremos y barreno central |
| **Disco-álabe** ×N | `<Name>_RotorStage{i}.stl` | Disco (barreno de ajuste, alma, rim siguiendo la línea de cubo, ranuras de abeto brochadas, círculo de tirantes) + casco del hub entre bandas de disco + álabes de rotor de la etapa i, cada uno con su raíz de abeto |
| **Anillo de carcasa** ×N | `<Name>_StatorRing{i}.stl` | Segmento del casco + vanos de estátor; brida apernada en CADA extremo de cada anillo (se apernan entre sí); el de la etapa de sangrado lleva el puerto |
| **Rotor / Carcasa** (vistas) | `<Name>_Rotor.stl` / `<Name>_Casing.stl` | Uniones en posición de marcha |
| **Ensamble** (vista) | `<Name>.stl` | Todo (solo inspección) |

Los planos de corte quedan a mitad del hueco entre el borde de fuga del
estátor de una etapa y el borde de ataque del rotor siguiente; cada disco
calza deslizante sobre el eje.

## Arranque rápido

El CLI tiene interfaz colorizada con arte ASCII, paneles y progreso en
vivo. Lo más fácil es el **asistente interactivo** (Enter acepta los
defaults):

```bash
python phyac_cli.py --interactive        # o -i
```

```bash
# 1. Diseño aerodinámico — humo rápido (~2 min, fidelidad L0):
python phyac_cli.py --pr 4.0 --mdot 25 --quick --fidelity L0 --outdir runs/smoke

#    Corrida completa (el through-flow L1 verifica los ganadores):
python phyac_cli.py --pr 4.0 --mdot 25 --rpm-max 18000 --utip-max 460 \
    --rtip-max 400 --rounds 5 --n-init 320 --outdir runs/pr4

#    Con pares de calibración CFD (cierra el lazo estilo Noyron):
python phyac_cli.py --pr 4.0 --mdot 25 --hifi-pairs cfx_pairs.json --outdir runs/cal
```

### Controles del ingeniero en el lazo

```bash
# Fijar variables de diseño y optimizar el resto (repetible), semilla propia:
python phyac_cli.py --pr 4 --mdot 25 --fix n_stages=5 --fix phi1=0.60 --seed 42

# Evaluar UN diseño dado (sin optimizar) con los entregables completos
# (las 10 variables de diseño — n_stages,RPM,HTR,phi1,psi_mid,psi_slope,
#  Rx,sigma_r,sigma_s,AR — el punto de operación sale de los flags):
python phyac_cli.py --eval-theta "4,12500,0.62,0.55,0.32,-0.10,0.60,1.20,1.10,2.20" \
    --pr 4 --mdot 25 --outdir runs/eval --map    # --map-vsv: schedule VSV auto
#   (12 valores = esos 10 + phi_slope,Rx_slope; 13/15 = θ completo.
#    Modo experto: --per-stage dist.json sobreescribe phi/psi/Rx por etapa)

# Reanudar una corrida desde su checkpoint (el dataset sustituye al LHS):
python phyac_cli.py --outdir runs/pr4 --resume --rounds 3

# Inspeccionar el frente de Pareto verificado de una corrida terminada y
# regenerar los entregables (reporte/geometría/STL) de cualquier punto:
python phyac_cli.py --outdir runs/pr4 --list-pareto
python phyac_cli.py --outdir runs/pr4 --pareto-pick 3 --stl
```

Control de interfaz: `--no-color` (texto plano, respeta `NO_COLOR`),
`--quiet`/`-q` (solo hitos). El estilo vive en
[`cli_style.py`](cli_style.py) (sin dependencias). Ejemplos completos:
`python phyac_cli.py --help`.

### 2. STL imprimibles desde el diseño optimizado (capa 5c)

#### Opción A: generación automática directa (recomendada)

```bash
# Todas las piezas con vóxel default de 0.5 mm:
python phyac_cli.py --pr 4.0 --mdot 25 --quick --stl

# Resolución custom (0.3–0.4 mm para piezas finales):
python phyac_cli.py --pr 4.0 --mdot 25 --quick --voxel 0.4
```

#### Opción B: invocación manual

```powershell
& "C:\Program Files\dotnet\dotnet.exe" build AxialCompressorDesigner.sln
& "C:\Program Files\dotnet\dotnet.exe" `
  "AxialCompressorDesigner.Example\bin\Debug\net9.0-windows\AxialCompressorDesigner.Example.dll" `
  runs\pr4\geometry\axial_compressor.json output\pr4 0.4   # json, outDir, vóxel mm
```

Salidas de la fase Python en `<outdir>/`: `report.html` (autocontenido;
frente de Pareto, tabla por etapa, carta de Smith, desglose de pérdidas,
triángulos, annulus, márgenes estructurales, BCs de CFD y trazabilidad —
más la sección matplotlib si está instalado; se omite con
`--no-figures`), `figures/`, `geometry/` (`axial_compressor.json` — el
contrato schema `phyac-axial-2` que consume la capa 5c —, `annulus.csv`,
`stage_summary.csv`, `bom.csv`, con `--step` el STEP de re-CAD de la
máquina — `parts/*.step` por pieza física con un álabe de muestra +
conteos en `parts/README.txt`, ensamble nombrado con
`--step-mode assembly`, la máquina jerárquica con las piezas reales
separadas con `--step-mode detailed`, o cada corona patronada con
`--step-mode full`; con `--check-interference` se comprueba por
booleanas y por pares la máquina ya montada — imprime el par culpable y
el volumen compartido y sale con código 3 si algo se solapa —,
`blade_stage0.step` si CadQuery está
disponible), `dataset.csv` (flywheel de datos), `phyac_run.json`
(checkpoint trazable), `phys_cache.jsonl` (caché persistente de
evaluaciones L1), `map.csv` con `--map`.

Salidas de la fase C#: un STL por pieza más las vistas de unión
(binarios, en mm).

## Mapa de módulos

| Módulo | Capa | Rol |
|---|---|---|
| `physics_core.py` | 1 | Núcleo físico multi-fidelidad: meanline L0 de stage-stacking (θ 15-D, correlaciones Lieblein/Howell/Koch, corrección de Reynolds por fila, $g(	heta)$ ×9, mapa fuera de diseño), calibración afín L2, caché, features físicos. |
| `scm_core.py` | 1 (L1) | Through-flow por curvatura de líneas de corriente: equilibrio radial completo con el término de curvatura meridional sobre 5–11 líneas, continuidad por tubo de corriente, cierre por ángulos del álabe, Euler por línea, pérdidas resueltas en el span. Sin dependencias externas. |
| `structures_core.py` | 1s | Núcleo estructural L0s: biblioteca de materiales con derating térmico, solver de disco rotatorio 1-D por etapa (validado contra Timoshenko), márgenes de raíz de álabe (K_t de Peterson), AN², estallido y Campbell (paso de álabes), métricas Goodman + cribado de flutter. **Restricción dura del optimizador** vía `DesignSpec.material`. |
| `neural_optimizer.py` | 2–4 | Deep ensemble con incertidumbre + quality gate, NSGA-II con dominancia restringida de Deb, adquisición LCB + k-means, orquestador `design(spec)` — portado de Phy-CC (núcleo agnóstico al dominio). |
| `blade_profiles.py` | 5a | Secciones NACA-65 / DCA sobre comba de arco circular, incidencia de diseño Lieblein/Aungier, desviación de Carter, punto fijo de ángulos metálicos. |
| `geometry_generator.py` | 5a | Secciones spanwise por free vortex (13 estaciones/fila), colocación axial exacta de filas, contrato `axial_compressor.json`, CSVs de annulus/etapas/BOM, BCs CFX, STEP opcional. |
| `data_pipeline.py` | 0 | Política de agregados públicos: anclas de correlación en el repo + manifiesto SHA-256 (sin descargas masivas — ver Science §6). |
| `report_generator.py` | 5b | Informe HTML autocontenido (SVG en Python puro: Pareto, carta de Smith, pérdidas, triángulos, annulus, trazabilidad). |
| `visualization.py` | 5b | Figuras matplotlib (mapa, carga por etapa, annulus, secciones) embebidas como PNG base64. Dependencia opcional con degradación elegante. |
| `AxialCompressorDesigner/` | 5c | Librería C# + PicoGK: importa `axial_compressor.json` (`PhyACImport`) y construye eje, discos-álabe y anillos de carcasa como STL imprimibles. |
| `AxialCompressorDesigner.Example/` | 5c | Ejecutable: CLI `axial_compressor.json → STLs`. |
| `phyac_cli.py` | producto | CLI end-to-end: espec → diseño → geometría → informe → dataset [→ STLs vía --stl/--voxel]. |
| `contract_schema.py` | 5a | JSON Schema publicado de `phyac-axial-2` + validador sin dependencias (también CLI: `python contract_schema.py <contrato>`). |
| `test_phyac.py` | VV&UQ | Suite de verificación: 171 checks (triángulos, conservación, continuidad de g, perfiles, esquema del contrato, solver de disco, núcleo del optimizador, through-flow L1, equilibrio radial, interferencias del ensamble). |
| `validation/` | VV&UQ | Campaña de validación vs NASA Stage 35, Rotor 37/67 y GE/NASA E³ HPC → `RESULTS.md`. |

## API Python (capas 1–5b)

```python
from neural_optimizer import DesignSpec, design
spec = DesignSpec(PR_target=4.0, massflow=25.0, RPM_max=18_000,
                  U_tip_max=460.0, r_tip_max_mm=400.0,
                  n_stages_max=8, material="Ti-6Al-4V")
out = design(spec)              # → {theta, record, pareto_front, history}
```

## API C# (capa 5c)

La superficie pública son **tres tipos**: la clase de parámetros, el punto
de entrada estático y el importador Phy-AC.

```csharp
using AxialCompressorDesigner;

// Desde una corrida Phy-AC (camino de producto):
var p = PhyACImport.FromContractJson(@"runs\pr4\geometry\axial_compressor.json");
string stlPath = AxialCompressorBuilder.Build(
    p, outputDir: @"C:\ruta\salida", showViewer: false);
// stlPath = ensamble; los STL _Shaft, _RotorStage{i}, _StatorRing{i},
// _Rotor y _Casing se crean al lado
```

La capa C# es deliberadamente **data-driven**: toda la matemática de
álabes (perfiles, ángulos metálicos, free vortex) vive en Python; el
contrato lleva las secciones de comba resueltas por estación de span y el
C# solo transforma, cose y voxeliza (ver
`.agent/axial-compressor-pattern.md`).

### Mapeo axial_compressor.json → Parámetros

| Bloque del contrato | Parámetro(s) C# | Regla |
|---|---|---|
| `annulus.hub` / `annulus.tip` | `HubLine` / `TipLine` | polilíneas [z, r] mm directas |
| `annulus.tip_clearance_mm` | `TipClearanceMm` | holgura de marcha en el extremo libre |
| `stages[i].rotor/stator` | un `RowParams` cada uno | nº de álabes, hueco axial, secciones |
| `sections[j].camber_points` | `SectionParams.CamberPoints` | superficie media de comba, marco de cuerda, centroide en el origen; engrosada con `voxMeshShell` (½·`thickness_mm`) |
| `sections[j].stagger_deg` | `SectionParams.StaggerDeg` | **firmado** (rotor +, estátor −); el builder aplica la rotación firmada uniformemente |
| `structural.drum_inner_r_mm` | `DrumInnerRadiusMm` | radio del eje bajo los discos |
| `assembly.*` | parámetros de eje/discos/bridas/sangrado | detalle de construcción (ver tabla) |

### Parámetros

Todas las distancias en milímetros. Defaults del bloque `assembly` del
contrato entre paréntesis.

| Parámetro | Default | Significado |
|---|---|---|
| `VoxelSizeMm` | 0.5 | Resolución del campo de vóxeles (0.8 borrador, 0.3–0.4 final) |
| `TipClearanceMm` | contrato | Holgura de marcha punta de rotor / hub de estátor |
| `RootSinkMm` | 1.5 | Hundimiento de la raíz del álabe en su cuerpo (soldadura) |
| `ShaftStubMm` | 30 | Muñón del eje a cada lado (asiento de montaje) |
| `ShaftBoreFrac` | 0.35 | Barreno central / radio del eje |
| `HubWallMm` | 5 | Pared del casco del hub (tambor, no tocho macizo) |
| `DiskWebFrac` | 0.5 | Espesor del alma del disco / cuerda axial (acotado 4–14) |
| `FlangeBoltCount` / `FlangeBoltDMm` | 12 / 6 | Círculo de pernos de ambas bridas |
| `FlangeWMm` / `FlangeTMm` | 12 / 8 | Ancho radial / espesor axial de brida |
| `BleedStage` | etapa media | Etapa (1-based) cuyo anillo lleva el puerto de sangrado |
| `BleedHoleDMm` / `BleedBossDMm` / `BleedBossHMm` | 18 / 32 / 14 | Barreno, diámetro y altura del boss del sangrado |

## Geometría (capa 5c)

Convención: **Z es el eje de rotación**; el gas fluye en +Z. Las filas se
construyen como **loft sólido del perfil real**: el contorno cerrado de
cada sección (distribución de espesor NACA-65/DCA de `points`) se vuelve
una malla estanca por álabe — paredes regladas entre costillas + tapas
hub/punta trianguladas LE→TE — voxelizada directa, conservando el radio
de LE y la asimetría presión/succión. Las filas con espesor < 2 vóxeles
(y los contratos legacy sin `points`) caen a la receta v3 de Phy-CC
(lámina de comba + `voxMeshShell`) con **WARNING en el log** cuando el
clamp de espesor actúa. Las filas **IGV y OGV** que la física asume están
en el contrato y cuelgan del primer/último anillo de carcasa. Las raíces
se hunden en su cuerpo; los extremos libres se retraen la holgura para
que el redondeo del vóxel no se coma el gap de marcha. Los huecos axiales
de fila se calculan con la envolvente **exacta** de la comba rotada de
cada sección (el hub corre a mucho menos stagger que la punta), así que
las piezas nunca se solapan — un test de regresión lo cuida.

## Estructura del proyecto

Ver [README.md](README.md#project-structure) (inglés) — el árbol es
idéntico.

## Dependencias

- **Capas 1–5b (Python ≥ 3.10)**: núcleo **solo NumPy** — L1 incluida:
  desde la fase 11 el through-flow por curvatura de líneas de corriente
  es propio (`scm_core.py`), así que no hay dependencia externa que
  instalar ni que pueda faltar en silencio. Opcionales con degradación
  elegante: `matplotlib` (figuras del informe), `cadquery` (STEP).
- **Capa 5c (C#)**: .NET 9 SDK + runtime x64. PicoGK se resuelve como
  paquete NuGet (`PicoGK` v2.2.0).

Licencias de terceros: ver [README.md](README.md#third-party-licenses).

## Verificación y validación

```bash
python test_phyac.py               # verificación: 171 checks
python validation/validate.py      # validación: máquinas NASA → RESULTS.md
python data_pipeline.py            # anclas de datos: rebuild + SHA-256
python contract_schema.py runs/x/geometry/axial_compressor.json   # contrato
python validation/parity_stl_step.py runs/x/geometry/axial_compressor.json 0.6   # STL vs STEP
python validation/bench_scm.py     # banco de L1: cobertura, coste, malla, residual
```

Con los extras instalados, `PHYAC_REQUIRE_STEP=1` y `PHYAC_REQUIRE_L1=1`
convierten en FALLO lo que si no sería un salto silencioso — es como lo
corre el CI, para que un runner sin CadQuery no pase en verde sin haber
comprobado la geometría.

**Verificación** (¿resolvemos bien las ecuaciones?): identidades de
triángulos, conservación de Euler/entalpía, límite isentrópico,
continuidad del vector g a través del choke, invariantes de
perfiles/contrato (conteos iguales, staggers firmados, sin solape axial),
solver de disco vs Timoshenko exacto, reglas de dominancia de Deb, puerta
del ensemble y reproducibilidad bit a bit (semilla 71).

**Validación** (¿las ecuaciones correctas?): califica el meanline contra
máquinas NASA medidas **sin recalibración por máquina** — el θ de cada
máquina reproduce su annulus publicado y su trabajo medido; se califica la
predicción de pérdidas → (η, PR). Tabla vigente en
[validation/RESULTS.md](validation/RESULTS.md):

| Máquina | ΔPR | Δη |
|---|---|---|
| NASA Stage 35 (etapa transónica) | +1.2% | +1.7 pts |
| NASA Rotor 37 (rotor transónico) | −1.2% | −1.6 pts |
| NASA Rotor 67 (fan transónico) | −1.2% | −2.5 pts |
| GE/NASA E³ HPC (10 etapas) | −5.6% | −1.6 pts |

## Estado actual y límites declarados (v0.1)

- **Parametrización por etapa** (fase 8): el θ 15-D lleva pendientes
  lineales frontal→trasero de φ y Rx (además de la de ψ de siempre) —
  el patrón de primer orden de las máquinas reales; la tolerancia del
  E³ se endureció 8%→6% con sus pendientes declaradas. Las
  distribuciones arbitrarias por etapa son el override experto
  `per_stage` (`--eval-theta --per-stage archivo.json`, solo meanline),
  deliberadamente fuera del espacio de búsqueda.
- **Modelo de choque L0**: promedio de choque normal en dos puntos con
  K_SHOCK = 0.70 calibrado en Rotor 37/67; sobre M_punta ≈ 1.5 es
  extrapolación.
- **L1 (`scm_core.py`)** resuelve la ecuación de equilibrio radial
  COMPLETA —con el término de curvatura— sobre 9 líneas de corriente, con
  continuidad por tubo de corriente y los ángulos metálicos del álabe
  congelados de la ley de torbellino de diseño. Es un peldaño de verdad:
  con vórtice libre reproduce el meanline con +0.1% en PR (verificación),
  y con torbellino controlado se separa varios por ciento porque ahí la
  forma cerrada del meanline es solo aproximada. Rechaza en vez de
  adivinar: una estación que no puede pasar el gasto, o un perfil de Cm
  que el limitador numérico sigue sujetando al converger, lanzan y el
  punto degrada a L0 con el motivo en `source`. Sin calificar todavía
  contra medida — el mismo hueco que el mapa fuera de diseño.
- **Álabes**: comba de arco circular con espesor NACA 65-010 / biconvexo,
  impresos como loft sólido del perfil (filas < 2 vóxeles caen a lámina
  de espesor uniforme con WARNING en el log); para álabes custom fieles
  a CFD, el camino es BladeGen/AGF de TD3.
- Materiales con valores típicos de handbook; sustituir por specs
  certificadas para diseño real. El STL de ensamble es solo inspección —
  imprimir las piezas.

## Historia

- **2026-08-17 — el mapa queda calificado contra medida (fase 12 ·
  F-02)**: la campaña de validación solo calificaba el PUNTO DE DISEÑO,
  mientras que el margen de bombeo se había convertido —desde la fase 9—
  en la restricción dura que más recorta el espacio. El AGARD AR-355
  §2.1.4.1 da dos números medidos del Rotor 37 al 100% como TEXTO, no
  como figura: `ṁ_choke = 20.93 kg/s` y un gasto de near-stall
  determinado experimentalmente en `ṁ/ṁ_choke = 0.925`. Entre los dos
  acotan todo el rango de gasto de la máquina. Un tipo de caso nuevo
  `OFFDESIGN` en `machines.py` califica ambos. **El criterio de bombeo
  PASA**: el modelo sitúa el stall al 90.4% del choke frente al 92.5%
  medido (−2.2%, objetivo ±3%) — el número que más pesaba con menos
  respaldo ya tiene alguno. **El de choke NO**: +6.6% de más, y de forma
  sistemática — el test de Mach AXIAL `MX_CHOKE = 0.78` declara choke
  mucho después de que la garganta del pasaje de un rotor transónico se
  haya bloqueado de verdad (M_rel de punta 1.49). El mapa queda desplazado
  hacia gastos altos por el lado del choke aunque el punto de bombeo esté
  bien. Arreglarlo pide un criterio de garganta con Mach relativo
  (G-09/F-07) y revalidar las cuatro máquinas, así que no se hace de paso;
  mientras tanto `--strict` corre contra una guarda interina declarada.
  Verificación 165 → 171 checks.

- **2026-08-16 — L1 pasa a ser un peldaño de verdad: through-flow propio
  por curvatura de líneas de corriente (fase 11 · F-01/H2)**:
  `scm_core.py` sustituye a `turbo-design`. Resuelve la **ecuación de
  equilibrio radial COMPLETA**, con término de curvatura meridional,
  sobre 9 líneas de corriente repartidas por fracción de gasto, con
  continuidad por tubo de corriente, los ángulos metálicos del álabe
  congelados de la ley de torbellino de diseño, **Euler por línea de
  corriente** (el trabajo deja de ser uniforme en el span) y pérdidas
  calculadas donde ocurren — choque solo donde el Mach relativo lo
  justifica, holgura en el 25% exterior del span, débito de pared de
  Koch & Smith en las bandas de pared. Lo que sustituye era
  `turbo-design` 1.4.2 con tres parches, forzado a `num_streamlines=1`
  porque con más su ODE de equilibrio radial se colgaba: L1 era otro
  meanline, lanzado en un subproceso con timeout, y muerto en silencio
  por una dependencia que el paquete no declara. **Verificación**: con
  curvatura nula el ODE integrado numéricamente reproduce la forma
  cerrada de la fase 9.1 con desvío 2e-16 en vórtice libre; sobre el θ de
  referencia L1 queda a +0.13% en PR de L0 con vórtice libre (donde el
  meanline es exacto) y a −6.5% con torbellino controlado (donde no lo
  es). Rechaza en vez de adivinar: una estación bloqueada, o un perfil de
  Cm que el limitador numérico sigue sujetando al converger, lanzan y el
  punto degrada a L0 con el motivo anotado. **El extra `l1` desaparece**
  — todo el lado Python vuelve a ser solo NumPy, L1 incluida, y corre en
  proceso a ~3 s por máquina. Un banco de pruebas
  (`validation/bench_scm.py` → `validation/BENCH_SCM.md`, 80 diseños por
  LHS) lo caracteriza: **85% de cobertura** con vórtice libre, 3.8 s de
  mediana por máquina (≈0.85 s por etapa), 0.33% de dispersión mediana
  entre 5 y 13 líneas de corriente, 15 iteraciones exteriores de mediana.
  Su hallazgo central es que la cobertura CAE con el número de etapas
  (100% en una, 71–75% en seis a ocho) por una razón estructural — L0
  dimensiona el annulus con su Cx uniforme, L1 resuelve un perfil, y el
  álabe de ángulo fijo convierte esa diferencia en trabajo, que cambia la
  densidad, que cambia la siguiente estación; sobre siete u ocho etapas se
  compone hasta ±27% de PR. Una ventana declarada `PR_WINDOW` (±15%)
  rechaza esos puntos en vez de devolver el número. Verificación
  153 → 171 checks.

- **2026-08-16 — la capa 5c vuelve a describir la misma máquina que la
  vía CadQuery (fase 10 · G-01)**: el STL (ruta de fabricación) y el STEP
  (ruta de re-CAD) se habían separado. Portado a C#: la **retención de
  abeto** (`FirTree.cs`, puerto línea a línea de `firtree_profile` y
  `_cq_firtree_solid` — perfil, loft de 5 secciones, brochado inclinado,
  plataforma que sigue la línea de cubo), el **recorte de cada fila
  contra la vena** en vez de posarla en los radios de su sección, las
  **holguras de marcha** con la convención del contrato, el **disco**
  (barreno de ajuste, alma, rim siguiendo la línea de cubo con su rebaje,
  ranuras brochadas, círculo de tirantes común, casco del tambor que se
  detiene en la banda del disco) y **`tie_bolt_count`/`tie_bolt_d_mm`**.
  Dos cosas corregidas por el camino: los anillos de carcasa del STL solo
  llevaban brida en los dos extremos de la máquina, así que los anillos
  impresos no se podían apernar entre sí; y la carcasa de CadQuery no
  tenía puerto de sangrado, que el contrato declara desde la fase 7.
  `validation/parity_stl_step.py` mide ya las dos rutas una contra otra
  — **eje −0.1 %, rotor −1.1 %, carcasa −0.7 %** en volumen y radios
  exteriores dentro de 0.4 mm (máquina de 2 etapas, vóxel 0.6 mm, filete
  de raíz apagado en las dos). La mitad barata de esa comparación corre
  en cada pasada de la suite y en el CI. Y midió algo que conviene
  señalar: el paso de filete de raíz de la 5c añade el **24 % del
  volumen de un anillo de carcasa** para un filete de 2 mm — ver
  docs/VALIDATION.md. Verificación 149 → 153 checks.

- **2026-08-16 — contrato versionado, chequeo de interferencias del
  ensamble y un CI que los cubre (fase 10)**: tres cosas que existían
  sobre el papel y que ninguna máquina comprobaba. (a) **El contrato se
  versiona**: `phyac-axial-2` con JSON Schema publicado
  (`schemas/phyac-axial-2.schema.json`), un validador sin dependencias
  (`contract_schema.py`, también CLI) ejercitado contra seis contratos
  rotos a propósito, y un lector de la capa 5c que **rechaza** una
  versión mayor que no conoce en vez de rellenar los huecos con defaults.
  La fase 9 había añadido abeto, tirantes, `vortex_n`, `bleed` y
  `sm_flow` conservando la etiqueta `phyac-axial-1`: un consumidor v1 no
  tenía forma de enterarse. (b) **Chequeo de interferencias del
  ensamble** (`--check-interference`, código de salida 3): booleana por
  pares sobre la máquina montada, reportando el par culpable y el volumen
  compartido. Encontró de inmediato tres defectos que el test por
  vértices sobre filas sueltas no podía ver — cascos de tambor
  atravesando de lado a lado el rim de los discos (53 cm³ cada uno),
  tirantes cruzando el alma de los discos delanteros porque cada disco
  calculaba su propio círculo de taladros, y rims de disco metidos dentro
  del álabe porque el rim era un cilindro mientras la línea de cubo sube.
  Los tres corregidos; el ensamble queda limpio con tolerancia de 1 mm³.
  (c) **El CI cubre por fin lo que cambió**: dependencias fijadas, un job
  con los extras `step` y `l1` donde `PHYAC_REQUIRE_STEP`/`_L1`
  convierten la ausencia de una dependencia opcional en FALLO en vez de
  en un salto silencioso, y un smoke del CLI de punta a punta que valida
  el contrato emitido. Ese job destapó que L1 llevaba muerto en
  silencio: `turbo-design` no declara `requests`, y `_scm_solve` lanzaba
  su worker con multiprocessing, que en Windows reimporta el `__main__`
  del que llama — así que cualquier script sin guard agotaba el timeout y
  degradaba a L0 sin decir nada. Ambos arreglados. Verificación 112 → 149
  checks.

- **2026-08-16 — física del rango y gas real (fase 9)**: la contabilidad
  de pérdidas se rehízo sobre mecanismo citado en vez de constantes
  ajustadas. (a) **Gas caloríficamente imperfecto**: cp(T), γ(T) y
  rendimientos EXACTOS para gas imperfecto vía la función de entropía
  phi(T) = ∫cp/T dT. (b) **Pérdidas por ENTROPÍA** (Dixon & Hall §5.5,
  ec. 5.4–5.9): la etapa se marcha con las presiones totales reales y el
  equivalente en trabajo de una pérdida es T₀₃·Δs — el viejo ω̄·½W₁²
  subestimaba las etapas traseras calientes un 15–25% cada una; además ω̄
  se refiere ya a la cabeza dinámica COMPRESIBLE (P₀−p), la definición de
  Koch & Smith. (c) **Correlación de stall de Koch 1981** en lugar de la
  constante `CH_STALL_MAX`: Ch_ef,stall desde el parámetro de difusión
  L/g₂ de la cascada con correcciones de Reynolds, holgura y espaciado
  axial, y con el factor de cabeza dinámica efectiva 𝔉_ef. (d) **Endwall
  de Koch & Smith 1976**: el bloqueo del annulus EMERGE de la carga
  (2δ*/g ∝ x³ + 2(ε/g)x) y trae su propio débito de rendimiento,
  jubilando el arrastre de annulus de Howell y el apaño `K_ENDWALL = 1.4`.
  **El mapa pasa a ser una herramienta**: el choke limita el gasto (las
  speedlines se hacen verticales en vez de desplomarse por debajo de
  PR = 1), la línea de trabajo de Dixon (ec. 5.26b) aporta el denominador
  que le faltaba a "margen de bombeo", y el margen en gasto es
  **restricción dura** (`SM_FLOW_MIN = 15%`). Además: el **sangrado**
  existe en la física (no solo como agujero en la carcasa) y el puerto se
  coloca detrás de la etapa que stallea primero; el **IGV** por fin paga
  su pérdida y su longitud; la **reacción de raíz** es la 9ª restricción
  dura; el surrogate se cortocircuita en fidelidad L0 (estaba aprendiendo
  a copiar una de sus propias entradas); se reporta la solidez fabricada
  frente a la optimizada; y los álabes de rotor reciben **raíz de abeto**
  con su ranura brochada en el disco por la vía STEP (el último punto
  donde turbodesigner iba por delante). Verificación 80 → 112 checks; las
  cuatro máquinas NASA siguen PASS con tres constantes libres menos
  (docs/VALIDATION.md §3).

- **2026-07-17 — STEP de ensamble (fase 8.2)**: `--step` exporta STEP
  de re-CAD en coordenadas de máquina (eje/hub/carcasa de revolución
  con bridas apernadas desde el contrato, álabe de muestra por fila con
  loft spline incl. IGV/OGV, `parts/README.txt` con los conteos para el
  patrón; `--step-mode assembly` para un cq.Assembly nombrado).
  CadQuery sigue siendo opcional con degradación silenciosa.
- **2026-07-17 — θ por etapa (fase 8)**: vector de diseño 15-D con
  pendientes lineales de φ/Rx añadidas AL FINAL (los índices del punto
  de operación no se mueven); los θ legacy de 13, checkpoints y anclas
  se paddean con pendientes 0 **bit-exacto** (anclas SIN re-congelar);
  override experto `per_stage` + `--per-stage`; la validación E³
  declara sus pendientes (φ cae, Rx sube hacia atrás) y su tolerancia
  de PR se endurece 8%→6% (−5.55%→−4.80%).
- **2026-07-17 — dinámica de álabes (fase 7)**: frecuencias propias de
  viga rotante (Southwell), **margen de Campbell frente al paso de
  álabes como 5º componente duro de g_struct**, K_t de Peterson del
  filete en la restricción de raíz (el MISMO radio de 2 mm que la capa
  5c ahora IMPRIME como fillet de raíz por closing morfológico
  restringido), Goodman (σ_alt admisible) y cribado de flutter V* como
  métricas reportadas, `figures/campbell.png` y bloque de dinámica en el
  reporte. Factibilidad estructural medida intacta (65% → 65% en
  LHS-500); anclas aero intactas.
- **2026-07-17 — off-design físico**: bucket de incidencia dependiente
  del Mach (semiancho 10° → 3.5° entre M 0.2 y 0.8, Aungier 2003; rama
  de choke 1.5× más ancha), desviación progresiva fuera de diseño
  Δδ = 0.30·i⁺ (Creveling 1968 — la ψ lograda cae hacia stall), y
  schedule automático de **estátores variables** en el mapa
  (`--map-vsv`, Δ = 50°·(1−N/Nd) tope 35° decayendo a media máquina) —
  sin VSVs las speedlines de velocidad parcial de máquinas PR ≳ 4 se
  declaran extrapolación no física. Punto de diseño exactamente
  invariante (anclas intactas).
- **2026-07-17 — holgura de punta por fila**: ε pasa a ser ABSOLUTA en
  mm (crece relativa a los álabes traseros que encogen — el efecto
  físico que un ε/h fijo borraba), con los regímenes de Sakulkaew 2013
  (óptimo bajo 0.8%, ~1.6 pts/1% lineal, punta descargada sobre 3.4%);
  física, contrato de geometría y validación comparten la misma ε (la
  holgura de marcha publicada de cada máquina NASA). Recalibración
  global documentada K_ENDWALL 1.0→1.4; el |Δη| máximo de las cuatro
  máquinas baja de 2.5 a 1.6 pts; REF_AX4 re-congelada.
- **2026-07-17 — álabes de perfil sólido**: la capa 5c hace loft del
  contorno cerrado de cada sección (`points` — la distribución de
  espesor NACA-65/DCA real) a un sólido estanco por álabe, conservando
  el radio de LE y la asimetría presión/succión que la lámina de espesor
  uniforme destruía; la lámina queda como fallback para filas < 2
  vóxeles y contratos legacy, ahora con **WARNING en el log** cuando el
  clamp de espesor actúa (antes era una distorsión silenciosa).
- **2026-07-17 — IGV/OGV**: las filas de guía que la física asume ya
  existen en el contrato y en las piezas 3D — el IGV (axial → pre-swirl
  α₁, delante del rotor 1) cuelga del primer anillo de carcasa y el OGV
  (α₁ residual → axial) del último; la desviación de Carter se
  generaliza con el signo del camber (χ₂ = β₂ − sgn(θ)·δ): las filas
  aceleradoras sobre-giran el metal más allá del ángulo de flujo
  objetivo.
- **2026-07-16 — corrección de Reynolds**: f_Re por fila en las pérdidas
  de fricción (nominal Koch & Smith 1976 a Re_c = 10⁶, exponente
  Wassell/Schäffler, rama laminar bajo 2×10⁵, sin crédito sobre 10⁶ —
  deltas de validación NASA sin cambio, ancla REF_AX4 re-congelada).
- **2026-07-16 — controlabilidad**: CLI con el ingeniero en el lazo
  (`--fix`, `--seed`, `--eval-theta`, warm start `--resume`,
  `--list-pareto`/`--pareto-pick`); `.gitattributes` corrige el fallo
  del hash de anclas por CRLF en clones Windows frescos.
- **2026-07-11 — v0.1**: proyecto creado desde el patrón Phy-CC (θ 13-D
  con radio de punta derivado de continuidad; conversión
  arrastre→pérdida del endwall y normalización de Koch corregidas durante
  la validación NASA; spike de integración TD3 axial con tres patches;
  construcción multi-parte de turborreactor — eje, discos-álabe, anillos
  de carcasa, bridas, puerto de sangrado — tras la revisión de STL que
  encontró solapes rotor/estátor por colocación a stagger medio, ahora
  exacta y con test de regresión).

## Licencia

Copyright 2026 Quasar Solutions Lab. Licenciado bajo
**[Apache License, Version 2.0](LICENSE)**. El software se distribuye
"TAL CUAL", sin garantías ni condiciones de ningún tipo — ver
[LICENSE](LICENSE).
