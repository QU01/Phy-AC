# 🗺️ Quasar Phy-AT — Hoja de ruta de desarrollo

**Estado**: planificación (la investigación previa está completa — ver
[PhyAT_Investigacion.md](PhyAT_Investigacion.md) y los seis informes de
`docs/phyat_research/`).
**Objetivo**: sistema autónomo de diseño inverso de **turbinas axiales** con
la lógica de Phy-AC — *spec de números → frente de Pareto verificado →
diseño recomendado → contrato de geometría → informe autocontenido →
STEP/STL imprimibles* — con escalera de fidelidades e incertidumbre
calibrada.

---

## 0. Definición de éxito (v0.1)

El equivalente exacto del "Spec → verified geometry → printable parts" de
Phy-AC:

> El ingeniero escribe **siete números** (relación de expansión objetivo ER,
> gasto, T₀₄, P₀₃, RPM máx, radio de punta máx, vida objetivo t_life) y
> recibe, en un laptop CPU y sin proponer geometría alguna: el frente de
> Pareto factible verificado, el diseño multietapa recomendado, el desglose
> de pérdidas y carga por etapa, las condiciones de contorno CFD, la BOM, y
> los STL imprimibles de cada pieza — eje, discos álabeados con abeto,
> segmentos de tobera y anillos de carcasa.

**Criterios de aceptación medibles** (los "targets" que cierran la v0.1):

| # | Criterio | Objetivo | Fuente del listón |
|---|---|---|---|
| A1 | Coste L0 | ≤ 1 ms/punto (NumPy/JAX, CPU) | Phy-AC: 0.5 ms; informe 01 §8.6 |
| A2 | Caudal en validación Kofskey | **< 2.5 % en el 100 % de los puntos** | métricas publicadas de TurboFlow (informe 01 §6.2) |
| A3 | Par/η en validación Kofskey | ±5 % en ≥ 90 % de los puntos (100 %/90 %/110 % de ω) | ídem |
| A4 | Ángulo de salida | ±5° en ≥ 95 % de los puntos | ídem |
| A5 | E³ HPT (refrigerada) | Δη ≤ ±2 pts sobre la definición termodinámica declarada | informe 01 §8.4 |
| A6 | Coste L1 (SCM) | ≤ 5 s/máquina, ≥ 80 % de cobertura en LHS | Phy-AC: 3 s, 85 % |
| A7 | Suite de verificación | ≥ 150 checks, todos deterministas, semilla reproducible | Phy-AC: 177 |
| A8 | Paridad STL↔STEP | volúmenes por conjunto dentro de ±1.5 % | Phy-AC: −0.1/−1.1/−0.7 % |
| A9 | Interferencias del ensamble | 0 solapes > 1 mm³ | Phy-AC fase 10 |
| A10 | Cero recalibración por máquina en la campaña de validación | invariante, no métrica | protocolo Phy-AC |

---

## 1. Invariantes heredados (no negociables)

1. **NumPy/JAX puro en el núcleo** — cero dependencias duras nuevas
   (matplotlib y CadQuery siguen opcionales con degradación anunciada).
   *Novedad respecto a Phy-AC*: el L0 se escribe con `jax.numpy` desde el
   día 1 (informe 06) — `import jax.numpy as jnp` con fallback declarado a
   NumPy si JAX no está instalado; el modelo es idéntico, JAX solo añade
   `grad`.
2. **g(θ) continuo y finito para θ degenerados** — la dominancia de Deb lo
   exige; toda saturación con `smooth_min`/interpolantes C² (Aungier), y
   todo clamp marcado en el código como lo que es.
3. **Fallo honesto** — excepciones tipadas con la razón registrada; nunca
   degradación silenciosa; distinción entre fallo numérico (degradar) e
   infactibilidad física (informar al optimizador).
4. **Anclas de regresión congeladas** — `--freeze-anchors` es un acto
   deliberado y citado; bit-exactitud con semilla fija.
5. **Contrato versionado con validador sin dependencias** — el consumidor
   rechaza la versión MAYOR que no entiende.
6. **VV&UQ como partes del producto** — verificación (¿resolvemos bien las
   ecuaciones?) separada de validación (¿las ecuaciones correctas?), campaña
   contra máquinas medidas sin recalibración.
7. **Los cinco "no trasplantar" de la investigación**: rangos de ψ/φ,
   f_Re de compresor, choke-como-fallo, `T_METAL_FRAC=1.00`, Carter.

---

## 2. Vista general de fases

Las estimaciones son **semanas-persona de esfuerzo** (no calendario), con
la incertidumbre propia de un plan previo al primer commit. Dependencias en
la columna "Dep.".

| Fase | Entregable principal | Est. | Dep. | Puerta de salida (gate) |
|---|---|---|---|---|
| **F0** Preparación | fuentes primarias resueltas, datos Kofskey en repo, esqueleto | 1–2 | — | DOIs [S]/[?] verificados; manifiesto SHA-256 de datos |
| **F1** L0 meanline | `physics_turbine.py` + `test_phyat.py` (núcleo) | 4–6 | F0 | Validación Kofskey fase 1 PASA (A2–A4) |
| **F2** Estructural 1s | `structures_turbine.py` | 2–3 | F1 | Timoshenko bit-exacto con ΔT=0; anclas de creep documentadas |
| **F3** Capas 2–4 | `neural_optimizer` adaptado + conformal | 2–3 | F1, F2 | puerta de calidad + cobertura conformal en benchmark sintético |
| **F4** Geometría 5a + contrato | `turbine_profiles.py`, `geometry_generator`, `phyat-axial-1` | 4–5 | F1 | validador + 6 contratos rotos rechazados; export TurboGrid |
| **F5** L1 SCM | `scm_turbine.py` + `scm_common.py` | 5–7 | F1, F4 | V1–V6 verificados; bench ≥80 % cobertura; Kofskey E-7776 |
| **F6** Capa 5c C# | `AxialTurbineDesigner` | 5–7 | F4 | paridad STL↔STEP (A8); interferencias (A9) |
| **F7** L2/L3 + flywheel | BCs de turbina, MULTALL interno, calibración | 3–4 | F1, F4, F5 | ≥15 pares semilla; guardarraíles de calibración activos |
| **F8** Producto | `phyat_cli.py`, informe HTML, `Quasar_PhyAT_Science.md`, CI | 3–4 | F1–F7 | run end-to-end en laptop < 5 min (quick) |
| **F9** Validación caliente | campaña E³ HPT + endurecimiento | 2–3 | F7, F8 | A5; límites declarados en README |
| — | **Total v0.1** | **31–44** | | los 10 criterios de aceptación |

El **camino crítico** es F0→F1→F4→F5/F6; F2, F3 y F7 pueden solaparse con
él. El primer hito visible ("smoke run": spec → diseño L0 verificado →
informe mínimo) llega al cerrar F1+F3, aproximadamente al 30 % del esfuerzo.

---

## 3. Fase 0 — Preparación (1–2 semanas)

La fase que evita pagar intereses después.

**F0.1 · Resolver las fuentes de segunda mano.** Los informes marcan
[S]/[B]/[M]/[?] todo lo no verificado contra fuente primaria. Con red
abierta, resolver los DOI críticos y archivar los PDF/valores:
- Kacker & Okapuu (1982) — el ±1.5 % y las tablas f_hub/y_te/Δφ².
- Benner 2006 I/II y 1997 — coeficientes completos (el código de TurboFlow
  ya los da; confirmar contra el paper).
- Pritchard 85-GT-219 — **los defaults de los 6 parámetros defaultables**.
- Gauntner TM-81453 — constantes del algoritmo de refrigeración.
- B del tip-clearance sin shroud (≈0.47 [?]) y curvas Larson-Miller por
  material (IN-718, MAR-M247, CMSX-4 de referencia).
- E³ LPT: localizar el número de CR en NTRS antes de comprometerlo.

**F0.2 · Datos de validación en repo.** Portar las geometrías Kofskey de
los YAML de TurboFlow (MIT) a `validation/machines.py` + los datos
experimentales, con manifiesto SHA-256 (patrón `data_pipeline.py`).
Descargar SPLEEN C1 (Zenodo) y las tablas de LS89/PAK-B que se usarán como
anclas L2.

**F0.3 · Esqueleto del repo.** Decisión: monorepo junto a Phy-AC o repo
hermano (recomendado: **repo hermano `Phy-AT`** con el mismo layout, como
Phy-AC lo es de Phy-CC). Ficheros vacíos con docstrings de contrato,
`pyproject.toml`, CI mínima, `.gitattributes` (la lección CRLF de Phy-AC).

**F0.4 · Congelar las convenciones.** Un documento corto
(`docs/CONVENCIONES.md`) que fija de una vez: ángulos desde el eje,
positivos antihorario; ψ = Δh₀/U²; Y referido a la dinámica de SALIDA;
estaciones 1/2/3; unidades SI internas + mm en el contrato. **Cada
correlación importada se audita contra este documento.**

**Gate F0**: fuentes críticas resueltas y archivadas; datos con manifiesto;
convenciones firmadas.

---

## 4. Fase 1 — L0 meanline (`physics_turbine.py`) (4–6 semanas)

El corazón del sistema. Orden interno (cada subfase con sus tests antes de
la siguiente):

**F1.1 · Gas** — `cp(T, FAR)`, `R(FAR)`, funciones de estado h/φ/s hasta
1900 K (extensión del patrón de la fase 9 de Phy-AC). Tests: rangos, límites
γ(T), consistencia con el JANAF de Phy-AC en aire frío.

**F1.2 · Triángulos y stacking** — (φ, ψ, R) → triángulos; marcha de
entropía en expansión; gasto variable por estación (¡desde el inicio!);
identidades tanα−tanβ=1/φ a 1e-6; conversión ax↔tan en **una** función
testeada.

**F1.3 · Cadena de pérdidas** — AMDC-KO (ajustes de Aungier) + Benner
(penetración Z_TE/H + secundaria + incidencia B97) + Dunham-Came (holgura)
+ f_Re de KO. Cada bloque con **hook de corrección aprendida** (mejora M10).
Tests: reproducir los valores del código de TurboFlow en 3–5 puntos por
correlación (el "oráculo" barato).

**F1.4 · Desviación y choking** — gauging arccos(o/s) + Aungier C²;
Mach crítico analítico corregido por pérdidas; ley de Stodola para
off-design; mapa (ER, ṁ√T/p) con la primera tobera como limitador.
Test: tobera bloqueada analítica (V6).

**F1.5 · Refrigeración** — Gauntner (x_factor como `cool_tech`), mezcla
de Hartsel, contabilidad Young & Wilcock (estátor/rotor separados, campo
`efficiency_definition`). Suavizado de la activación (T_gas<T_metal).

**F1.6 · Restricciones g** — grupos A (validez de correlaciones), B
(aerodinámicas, incl. margen de choking), C (mecánicas — placeholder hasta
F2), D (spec). Test de continuidad de g a través del choke y de la
activación de refrigeración (barrido de 81 puntos, como Phy-AC).

**F1.7 · Validación fase 1 (Kofskey)** — `validation/validate.py` con el
protocolo Phy-AC (θ reconstruido de la geometría publicada, cero
recalibración): D-7625 (estátor estrangulado), D-6967 1 etapa (rotor
estrangulado), D-6967 2 etapas, E-7776 (barrido de garganta). Congelar las
primeras **anclas de regresión** al pasar.

**Gate F1**: A1–A4 cumplidos; ~60–80 checks de verificación; anclas
congeladas y citadas.

---

## 5. Fase 2 — Estructural 1s (`structures_turbine.py`) (2–3 semanas)

- **F2.1** Biblioteca de materiales calientes **por pieza** (disco/álabe/
  NGV/carcasa): IN-718, Waspaloy, IN-713LC, MAR-M247, CMSX-4 de referencia;
  derating de σ **y de E**; T_max; coeficientes Larson-Miller con fuentes
  documentadas en `validation/`.
- **F2.2** `solve_disk` heredado + término termoelástico αE·dT/dr (misma
  tridiagonal). **Test de regresión: ΔT=0 reproduce Timoshenko bit a bit.**
  Perfil T(r) paramétrico (bore→rim).
- **F2.3** Las 9 g: fluencia y burst a **1.20·N** (14 CFR 33.27, citado en
  el código), AN² derivado del creep, raíz con K_t + carga de shroud,
  **creep de álabe y de rim** (a velocidad de diseño, con t_life), T_metal
  por la cadena térmica (OTDF/RTDF → T₀_rel → η_c → TBC), Campbell (nº NGV
  y nº inyectores; λ₁ empotrado-apoyado si shroud), poste del abeto.
- **F2.4** Métricas reportadas: margen de temperatura, Goodman, flutter,
  vida a la T de pico del OTDF; LCF declarado como no modelado.

**Gate F2**: Timoshenko bit-exacto; LHS(500) sin excepciones y con g
finito; factibilidad estructural medida y documentada.

---

## 6. Fase 3 — Capas 2–4 (optimizador) (2–3 semanas, solapable con F4)

- **F3.1** Port del core de `neural_optimizer.py` (es agnóstico de dominio
  — el mismo trasplante Phy-CC→Phy-AC). Adaptador de dominio: rangos del θ
  de 18-D, `fix_operating_point`, embedding de turbina (M5 del informe 04:
  Zweifel×2, M₂/M₃rel, margen de garganta continuo, R_hub mín, deflexión
  máx, Re mín, τ/h, AN², t_TE/o, swirl de salida, Smith (φ,ψ), log ER y
  η_L0).
- **F3.2** Mejoras UQ baratas: **conformal split** sobre el ensemble (M1,
  ~50 líneas; q_α recalibrado por ronda y cobertura empírica en el log),
  K=8–10 miembros (M2, medir antes/después), recalibración isotónica como
  métrica del informe (M3).
- **F3.3** Adquisición: LCB+k-means heredado como default; ΔHV esperado por
  MC (M6) y sesgo a la frontera g≈0 (M8) tras un A/B en benchmarks
  sintéticos — solo se adoptan si ganan.
- **F3.4** Calibración L2: afín ponderada y regularizada hacia la identidad,
  estratificada por régimen; extensión afín+2 features (M4) con puerta LOO.

**Gate F3**: puerta de calidad + cobertura conformal verificadas en función
sintética; primer **run end-to-end L0** (spec → Pareto verificado) en
laptop.

---

## 7. Fase 4 — Geometría 5a y contrato (4–5 semanas)

- **F4.1** `turbine_profiles.py`: **Pritchard 11-par** con renderizador
  dual (`PRITCHARD11` analítico / `BSPLINE_G2` opcional); correlaciones de
  los 6 defaults (Zweifel→paso, o=s·cosα₂ corregido, δ_ug, ε_i, r_LE,
  t_TE/o); invariante de no autointersección; g geométricas (δ_ug≤20°,
  o/s∈[0.15,0.85], t_TE≤0.15·o, M≤1.15 o `CD_SUPERSONIC_REQUIRED`).
- **F4.2** 3D: apilado por centroide + compound lean paramétrico de 3
  números + sweep; ley de vórtice con `vortex_n`; garganta **por sección**;
  17/21 estaciones con densificación coseno; ≥100 puntos/sección con
  clustering LE/garganta/TE.
- **F4.3** **Contrato `phyat-axial-1`**: schema JSON publicado + validador
  sin dependencias (CLI), con los bloques nuevos (`thermal`, `interfaces`
  con requires/provides, `manufacturing`, `provenance`, `profile_params`,
  `index_markers`); 6 contratos deliberadamente rotos como tests; anillos,
  BOM, CSVs.
- **F4.4** Exportadores: TurboGrid (`hub/shroud/profile.curve`), CadQuery
  STEP (reutilizando `_cq_blade_trimmed` y compañía), BCs CFD (paquete
  completo del informe 03 §8.3 — perfiles radiales, Tu+L_t, térmica de
  pared, refrigeración por hilera, purgas con swirl ratio, definición de η).

**Gate F4**: contrato validado; sin solapes axiales (test de envolvente
rotada); STEP de una etapa inspeccionado; export TurboGrid cargable.

---

## 8. Fase 5 — L1 SCM (`scm_turbine.py`) (5–7 semanas)

El orden del informe 02 §8.9, cada paso con su verificación:

1. **Esqueleto**: `scm_common.py` extraído de `scm_core.py` (reposicionado,
   curvatura, mezcla spanwise, streamtubes, estado termodinámico);
   estaciones LE/GARGANTA/TE (3n+1); rotalpía en rotor; cierre
   arccos(A_t/A_o) sin desviación. → **V1, V2 (T21 reutilizado), V6.**
2. **Curvatura con flare**: eliminar `CURV_MAX` (suavizado de r'' en su
   lugar). → **V3 (actuator disc)** — el test del término de curvatura.
3. **Choking**: M_crit corregido por pérdidas por línea; límite de gasto
   por streamtube con redistribución spanwise; selección de rama
   sub/supersónica; registro de filas bloqueadas y fracción de span.
   → **validar contra Kofskey E-7776** (barrido de garganta).
4. **Pérdidas spanwise**: KO por línea (sin f_hub) + bandas Z_TE/H de
   Benner + Moustapha; fuga sin shroud en banda + defecto de trabajo; fuga
   con shroud como fuente/sumidero. → validar contra D-6967.
5. **Desviación de Aungier** por línea.
6. **Refrigeración en tres lazos** (primitiva `inject()`). → **V4 (balance
   de energía cerrado como assert)** antes de cualquier prestación.
7. **Mezcla spanwise** recalibrada (Lewis 1994); **V5** (independencia de
   N_sl con flare).
8. Hueco de compound lean (F_r con λ≡0, documentado).
9. **Taxonomía de fallo** completa (TurbineChoked, MassFlowInfeasible con
   el gasto alcanzable, LimitLoading, NegativeReactionAtHub,
   EnergyImbalance…), y la ventana L0↔L1 sobre el **gasto reducido**.
10. `bench_scm_turbine.py` (80 LHS): cobertura, coste, sensibilidad de
    malla, residuales → `BENCH_SCM.md`.

**Gate F5**: V1–V6 en verde; A6 (≤5 s, ≥80 % cobertura); las cuatro
máquinas Kofskey resueltas también a L1.

---

## 9. Fase 6 — Capa 5c C# (`AxialTurbineDesigner`) (5–7 semanas, solapable con F5)

- **F6.1** Fork estructurado de la base de Phy-AC (~60 % reutilizable):
  `PhyATImport` (rechazo de MAYOR desconocida), `BladeRow` (loft sólido —
  camino normal con TE romo), recorte contra la vena, `RotorDrum.voxDisc`,
  `Casing` con bridas dobles, `LatticeUtils`, `ParityReport`.
- **F6.2** Piezas nuevas, por prioridad: **`TurbineFirTree`** (multi-lóbulo,
  poste del disco, agujeros de refrigerante) y **`NozzleSegment`** (bandas +
  feather seals + huecos térmicos) → `TipShroud` (rieles, Z-notch, panal en
  carcasa) → `InterstageSeal`, `CoverPlate`, `RimCavity`.
- **F6.3** **`ManufacturabilityCheck`**: rasgos de sello < 2 vóxeles se
  **rechazan o marcan** `feature_suppressed` — nunca clamp silencioso.
- **F6.4** Paridad y ensamble: `parity_stl_step.py` adaptado; chequeo de
  interferencias par a par (exit 3); el orden de piezas empieza en un
  NozzleSegment (planos de corte desplazados media etapa).
- **F6.5** Bloque `manufacturing` del contrato materializado en el README
  de export (no es pieza caliente; sin refrigeración interna; qué demuestra
  y qué no).

**Gate F6**: A8 y A9; una máquina de 2 etapas impresa "en seco"
(voxel 0.6) revisada visualmente; suite T20-equivalente (perfil de abeto
Python↔C# punto a punto).

---

## 10. Fase 7 — L2/L3 y data flywheel (3–4 semanas, solapable)

- **F7.1** `HiFiCalibration` de turbina: KEYS = (**Γ primero**, η_tt, ER,
  α_exit, Λ); ponderación por σ del par; regularización hacia la identidad;
  estratificación HPT/LPT y bandas de Re; guardarraíles (|a−1|>0.15 →
  rechazo como error de definición; congelación fuera de la cáscara de los
  pares; suelo de 0.25 %).
- **F7.2** **MULTALL como L3 interno opcional**: generador de
  `meangen.in`/`stagen.in` desde el contrato → ejecución → parser de salida
  a pares (con perfiles radiales). Contrato de retorno con y⁺/malla/flag de
  convergencia obligatorios (sin ellos, el par no entra).
- **F7.3** Pares semilla: Aachen 1.5 (smoke test del pipeline) → SPLEEN C1
  → LS89 → E³ (geometría+malla+BC de data.gov) → objetivo ≥15–25 pares
  para activar el residual no lineal.
- **F7.4** Calibración por gradiente de los coeficientes de pérdidas
  (`jax.grad` sobre L0/L1) con regularización hacia los valores publicados —
  la mejora que sustituye a la recta de 2 parámetros; puerta LOO.

**Gate F7**: lazo completo demostrado — un par MULTALL generado, ingerido,
y la calibración medida contra un hold-out.

---

## 11. Fase 8 — Producto (3–4 semanas)

- **F8.1** `phyat_cli.py`: wizard interactivo, `--quick`, `--fidelity`,
  `--fix/--seed/--eval-theta/--resume/--list-pareto/--pareto-pick`,
  `--stl/--voxel/--step`, `--map` (característica ER-ṁ), `--hifi-pairs`.
- **F8.2** `report_generator`/`visualization`: diagrama de Smith de turbina
  con el diseño superpuesto (guardrail visual), triángulos, anillo
  meridional, desglose de pérdidas por bloque, característica con la fila
  estranguladora marcada, Campbell, márgenes de creep (margen de
  temperatura en primer plano), trazabilidad.
- **F8.3** `Quasar_PhyAT_Science.md`: ecuaciones, pseudocódigo y
  referencias — solo con fuentes ya verificadas (F0.1); límites declarados
  (wake recovery, clocking fuera de alcance, techo de ganancia 3D 0.5–1.0
  pt como aviso).
- **F8.4** CI: suite + validación en estricto; `PHYAT_REQUIRE_STEP/_L1`;
  smoke end-to-end con validación del contrato emitido.

**Gate F8**: run quick < 5 min en laptop; informe autocontenido; CI verde.

---

## 12. Fase 9 — Validación caliente y endurecimiento (2–3 semanas)

- **F9.1** Campaña **E³ HPT** (CR-168289 + hardware 19850002687):
  refrigerada, con `efficiency_definition` declarada; tolerancia ±2 pts.
- **F9.2** Fase 3 de validación si hay tiempo: series de alta temperatura
  (19680006274, 19720024133), potencia libre (19790009688).
- **F9.3** Endurecimiento: fuzzing de θ degenerados (LHS 500 con g finito),
  `--strict`, documentación de límites v0.1 en el README (el bloque
  "Current Status and Declared Limits" de Phy-AC).

**Gate F9 = release v0.1**: los 10 criterios de aceptación de §0.

---

## 13. Backlog post-v0.1 (deliberadamente fuera)

| Ítem | Disparador para abordarlo |
|---|---|
| Compound lean activo en L1 (F_r con λ≠0) y en el optimizador | cuando la 5a emita lean y haya pares 3D que lo midan |
| Endwall no axisimétrico (`endwall.non_axisym`) | bloque ya reservado en el contrato; requiere pares CFD dedicados |
| Perfiles convergente-divergentes (MoC) | demanda real de M₂ > 1.15; hoy `CD_SUPERSONIC_REQUIRED` |
| Q3D de transición propio (tipo MISES) en L2 | cuando el lapso de Re de LPT domine el error de calibración |
| Co-kriging recursivo (Le Gratiet) en L2 | ≥15–50 pares hi-fi acumulados |
| SU2 con adjunto como vía de optimización de forma | cuando haga falta gradiente sobre la geometría 3D |
| Interfaz Phy-CB activa (`interfaces.upstream.system = "phycb-annular-1"`) | cuando Phy-CB exista — hoy defaults declarados de OTDF/RTDF/swirl |
| LPT high-lift (Zweifel > 1.15, Coull–Hodson / Zhu–Sjolander) | demanda de LPT; hoy restricción Zweifel ≤ 1.15 |
| Multi-spool / acoplamiento con Phy-AC (turbina que mueve al compresor) | tras v0.1, junto al matching de ciclo |
| JAX-Fluids como banco de priors 2D; JAX-FVM/DiFVM como L1.5 | madurez de esas librerías (informe 06: vigilar) |

---

## 14. Riesgos del plan y mitigaciones

| Riesgo | Prob. | Impacto | Mitigación |
|---|---|---|---|
| Coeficientes de correlaciones mal transcritos | media | alto | F0.1 (fuentes primarias) + tests-oráculo contra TurboFlow en F1.3 |
| Convención de ángulos inconsistente entre módulos | alta si no se ataja | alto | F0.4 + función única + tests de identidad — desde el primer commit |
| El choking rompe el solver L1 (doble raíz) | media | alto | es el paso 3 de F5 con validación dedicada (E-7776); plan B: formulación Newton-acoplada (T-AXI/MTFLOW) |
| La deriva L0↔L1 sobre el gasto es peor de lo previsto | media | medio | `MassFlowInfeasible` informativo + ventana de gasto reducido + re-dimensionado de garganta en el lazo |
| JAX como dependencia incomoda al usuario de laptop | baja | medio | fallback NumPy transparente; JAX solo se exige para calibrar por gradiente |
| El esfuerzo de la 5c se subestima (NGV segmentados) | media | medio | F6 solapada con F5; el STL de demostración puede salir con shroud/sellos en v0.2 |
| E³ HPT no pasa a ±2 pts | media | medio | es la fase 2 de validación: se diagnostica por bloques de pérdida (hooks M10) antes de tocar constantes globales |
| Alcance-fluencia ("ya que estamos…") | alta | alto | el backlog §13 existe para decir que no; cada fase tiene gate escrito |

---

## 15. Recursos

- **Datos**: Kofskey (NTRS + YAML de TurboFlow, MIT), E³ HPT (CR-168289 +
  data.gov), SPLEEN C1 (Zenodo), LS89, PAK-B/T106, Aachen (tutorial SU2).
- **Software de referencia**: TurboFlow y AxialOpt (oráculos de
  correlaciones, MIT), MULTALL (L3 interno), NASA TD3 (contraste
  Craig-Cox/Traupel), pyturbo-aero/ParaBlade (referencias de perfil),
  CadQuery, PicoGK 2.2.0, .NET 9.
- **Hardware**: laptop CPU para todo el lazo; GPU opcional solo si el
  usuario corre CFD externo moderno (Turbostream/SU2) — nunca requerida.
- **Personas**: el plan asume el patrón Phy-AC — un desarrollador principal
  con revisión; las semanas de §2 son esfuerzo, no calendario.

---

## 16. Primeros 10 pasos concretos (arranque)

1. Crear el repo `Phy-AT` (o el directorio hermano) con esqueleto y CI.
2. Escribir `docs/CONVENCIONES.md` y congelarlo.
3. Resolver y archivar los 6 DOI críticos de F0.1.
4. Portar las 3 geometrías Kofskey + datos experimentales con manifiesto.
5. `physics_turbine.py`: gas cp(T,FAR) + tests de rango.
6. Triángulos (φ,ψ,R) + identidades a 1e-6 + conversión de ángulos única.
7. Y_p de KO sobre ajustes de Aungier + test-oráculo contra TurboFlow.
8. Benner completo (penetración + secundaria + incidencia) + oráculo.
9. Gauging + Aungier C² + Mach crítico + test de tobera bloqueada.
10. Primera pasada de `validate.py` contra Kofskey D-7625 — y mirar el
    caudal: **si no sale a <2.5 %, parar y depurar antes de seguir**.

---

*Hoja de ruta derivada de la investigación completa
([PhyAT_Investigacion.md](PhyAT_Investigacion.md), informes
[01](phyat_research/01_meanline_perdidas_turbinas.md) ·
[02](phyat_research/02_scm_throughflow_turbinas.md) ·
[03](phyat_research/03_sota_cfd_altas_fidelidades.md) ·
[04](phyat_research/04_ml_pinns_uq_optimizacion.md) ·
[05](phyat_research/05_geometria_estructuras_turbina.md) ·
[06](phyat_research/06_gpu_diferenciable_xlb_jaxfluids_warp.md)) y del
patrón de fases documentado en el historial de Phy-AC.*
