# Informe: Generación de geometría y núcleo estructural para turbinas axiales — sistema **Phy-AT (Quasar)**

**Análogo de Phy-AC (compresores axiales) · capas 5a / 5c / 1s** · Fecha: 2026-08-19

---

## 0. Método, alcance y honestidad sobre las fuentes

**Leyenda de verificación:**

| Marca | Significado |
|---|---|
| **[V]** | Verificado: documento o código fuente leído en la sesión de investigación |
| **[B]** | Buscador: resumen de búsqueda con URL de la fuente primaria; fiable en lo cualitativo, **números a re-verificar** |
| **[M]** | Memoria / referencia estándar: **comprobar antes de codificar como constante** |

Dominios bloqueados en la sesión: ASME, NTRS, ScienceDirect, MDPI, arXiv, TU Delft, DTU Orbit, archive.org, picogk.org, GPPS. Presupuesto de búsqueda agotado antes de cubrir dos temas menores (endwall contouring paramétrico y gradiente térmico de disco), documentados desde [M].

**Base Phy-AC leída directamente del repo [V]**: `contract_schema.py`, `schemas/phyac-axial-2.schema.json`, `blade_profiles.py`, `geometry_generator.py`, `structures_core.py`, `.agent/axial-compressor-pattern.md`, `docs/Quasar_PhyAC_Science.md` §8.

---

## 1. Parametrización de perfiles de turbina

### 1.1 Por qué la familia NACA-65/DCA de Phy-AC **no se puede reutilizar**

`blade_profiles.py` construye espesor NACA 65-010/DCA sobre comba de arco circular con Lieblein + Carter [V] — específico de cascadas de **difusión** (comba 10–45°, t/c 0.04–0.10, TE afilado). Una cascada de turbina es lo contrario:

- **Acelera**: tolera giros de 60° (NGV) a 100–120° (rotor de impulso) [M].
- **Gruesa**: t_max/c ≈ 0.15–0.25, máximo adelantado (x/c ≈ 0.3–0.4) [M].
- **TE grueso** (estructura + refrigeración): espesores de TE medidos de **5–20 % del ancho de garganta**; un TE redondo del ~11 % de la garganta produce **60 % más pérdida** que el afilado [B — NASA TM 19720007336].
- **El ángulo de salida lo fija la GARGANTA**, no el metal del TE — Carter carece de sentido.
- Succión con dos zonas: **guiada** (hasta la garganta) y **no guiada / uncovered turning** aguas abajo [V — es la estructura del método de Pritchard].

**Lo que sí se reutiliza de Phy-AC es el *contrato de salida* del perfil** (polígono cerrado CCW + camber + stagger firmado), no su generación.

### 1.2 Pritchard (1985) — el modelo de **11 parámetros** (la referencia canónica)

**Fuente**: L. J. Pritchard (Williams International), *"An Eleven Parameter Axial Turbine Airfoil Geometry Model"*, **ASME 85-GT-219** (1985) [B]. Del abstract: dado el número de álabes y, a cada radio, cuerda axial y tangencial, garganta, giro no guiado, radios de LE y TE, ángulos metálicos de entrada/salida y cuña de entrada, **ambas superficies quedan descritas por funciones analíticas**; **6 de los parámetros tienen valores por defecto**.

**Los 11 parámetros, verificados contra la implementación de referencia** [`DavidPoves/11-Parameters-Turbine-Blade-Generator`](https://github.com/DavidPoves/11-Parameters-Turbine-Blade-Generator) [V]:

| # | Parámetro | Símbolo | ¿Defaultable? |
|---|---|---|---|
| 1 | radio del cilindro | r | no (del span) |
| 2 | número de álabes | Z | no |
| 3 | cuerda axial | c_x | no |
| 4 | cuerda tangencial | c_t | sí |
| 5 | ángulo metálico de entrada | β_i | no |
| 6 | ángulo metálico de salida | β_o | no |
| 7 | garganta | o | **sí** (de o/s) |
| 8 | giro no guiado | δ_ug | **sí** |
| 9 | semiángulo de cuña de entrada | ε_i | **sí** |
| 10 | radio de LE | r_LE | **sí** |
| 11 | radio de TE | r_TE | **sí** |

**Construcción [V]**: 5 puntos elementales → arcos de LE y TE → **cierre iterativo de la garganta** (una iteración escalar sobre ε_o, tol 1e-5) → **cúbicas** en presión y succión con ángulos metálicos como condiciones de contorno → **arco circular de giro no guiado** entre garganta y TE. Succión = [LE] + [cúbica guiada] + [arco no guiado] + [TE]; presión = [LE] + [cúbica] + [TE]. **Todo analítico, sin optimización.**

**Debilidad documentada [B]**: G1 pero **no G2** en LE y garganta → picos de velocidad en CFD fino. Irrelevante para geometría impresa.

### 1.3 Parametrizaciones B-spline / NURBS

**(a) ParaBlade** — Agromayor, Anand, Müller, Pini, Nord, *"A Unified Geometry Parametrization Method for Turbomachinery Blades"*, **Computer-Aided Design 133 (2021)** [B]; [repo](https://github.com/NAnand-TUD/parablade) [V]. 10 variables por sección (stagger, θ_in/out, radios LE/TE, dist_in/out, 6+6 espesores); comba = B-spline cúbica de 4 puntos de control (**G2**); espesor como offset normal; 3D con canal meridional de 4 B-splines + leyes spanwise B-spline; sensibilidades por paso complejo; problema inverso [V]. Cubre "reaction and impulse turbine blades".

**(b) pyturbo-aero (NASA Glenn)** — [repo](https://github.com/nasa/pyturbo-aero), licencia NASA-1.3 [V]. `Airfoil2D(alpha1, alpha2, axial_chord, stagger)`; comba Bézier de 3 puntos; espesores ss/ps con flag **`aft_loaded`** (turbina vs compresor); **TE con arco y cuñas independientes** (`te_create(radius, wedge_ss, wedge_ps)`); `match_thickness()` iguala curvatura en LE (resuelve el no-G2, con minimización); **`channel_get()` calcula la garganta contra el álabe vecino** [V]. 3D: apilado LE/TE/**centroid**, `add_sweep()`, `add_lean()` (puntos de control en la Bézier de apilado), splines PCHIP spanwise (correspondencia punto a punto 1:1 — la invariante del loft 5c), export STL/JSON. **Sin STEP, fillets, shroud ni tip clearance** [V]. Aviso del README: *"solo para exploración de diseño; el componente final siempre en CAD"*.

**(c) T-Blade3 (GTSL/UC, usado por NASA)** — [repo](https://github.com/GTSL-UC/T-Blade3) [V]. Comba por **B-spline de la segunda derivada, integrada analíticamente dos veces** → continuidad de curvatura por construcción; espesor Wennerstrom / B-spline cuártica / NACA modificada con t_TE y su derivada; LE/TE por B-splines; Fortran, integra con ESP para el sólido CAD [V].

**(d) CST/Kulfan** — [Kulfan 2007](http://servidor.demec.ufpr.br/CFD/bibliografia/aerodinamica/kulfan_2007.pdf) [B]. Clase × forma; radio de LE, espesor de TE y boat-tail ligados a los extremos de la función de forma; usado por Clark para poblar bases de datos de HPT [B]. **Sus coeficientes no son cantidades de ingeniería** — no sirve como generador directo sin optimización; sí como representación secundaria.

### 1.4 Comparativa

| Método | ¿Entradas del meanline? | Continuidad | ¿Familia válida sin optimizar? | Coste |
|---|---|---|---|---|
| **Pritchard 11-par** | **Sí, todas** | G1 | **Sí** (analítico + 1 iteración escalar) | ~µs |
| ParaBlade camber+thickness | parcial | **G2** | sí con espesores de familia | ~ms |
| pyturbo-aero | mayoritariamente | G2 con minimización | sí | ~ms |
| T-Blade3 | no (curvatura) | G2/G3 | con plantillas | ~ms |
| CST | no | C∞ | no | µs+ajuste |

### 1.5 Diseño de la **garganta** (o/s) — el corazón del perfil de turbina

**Regla del gauge angle** [B]: α_gauge = arccos(o/s) (desde el eje). Desviación δ = α_gauge − α₂,real. Correlaciones: AM 1951 (α₂ = f(o/s, s/e) con e = radio del dorso post-garganta); Dunham & Came; Kacker & Okapuu; **Aungier 2006** (función a trozos o/s→ángulo) [B]; actualización moderna *A Reliable Update of the AM Profile and Secondary Correlations* [B]. En bloqueo: **el ángulo sale de conservación de masa**, no de correlación [B].

**Implicación arquitectónica**: en Phy-AC el flujo es (β₁,β₂)→(χ₁,χ₂)→perfil. **En Phy-AT el orden se invierte**: (α₂ objetivo, s) → **o = s·cosα₂ corregido** → y *o* es una **entrada de Pritchard**. La garganta es la variable primaria; el metal de salida es casi decorativo.

**Paso por Zweifel (1945)**: Z_w = 2(s/b_x)cos²α₂(tanα₁+tanα₂); coeficiente **0.8–1.1** según aplicación [B]; default Z_w ≈ 0.8. **Zweifel sustituye a σ_r/σ_s como variable de diseño.**

### 1.6 TE, cuñas, unguided turning

- **t_TE**: rango medido 5–20 % de la garganta; +60 % de pérdida al 11 % [B]. Default: t_TE/o ≈ 0.05 (demostrador no refrigerado); 0.08–0.12 con eyección de TE.
- **Cuñas** [M]: ε_i 5–15°; ε_o **no libre** (sale del cierre de garganta) [V].
- **Unguided turning** [B]: restricción de primer orden; dato verificado: **19.47°** en un diseño NASA real (Flegel, [NTRS 20140012427](https://ntrs.nasa.gov/api/citations/20140012427/downloads/20140012427.pdf)). Default ~8°, **restricción dura δ_ug ≤ 20°**.

### 1.7 Perfiles supersónicos convergente-divergentes

Para M₂ > ~1.2–1.3, garganta + divergente por método de las características (Deich; Goldman/NASA) [M]. **Recomendación v1: NO implementar MoC.** Si M₂,is > 1.15 en alguna sección → warning en el contrato y `profile: "CD_SUPERSONIC_REQUIRED"` con `supported: false`. Declarar el límite en vez de generar un perfil que la física no respalda.

### 1.8 Decisión

**Ganador: Pritchard**, por (1) **la garganta y el giro no guiado son ENTRADAS** — en cualquier otro método hay que generar y luego *medir* la garganta (optimización 2D encubierta); (2) **analítico y determinista** (µs, sin modos de fallo de convergencia); (3) sus 6 defaults son exactamente la lista de correlaciones de diseño a escribir — el análogo de "Lieblein + Carter".

**Mitigación del no-G2**: **dos renderizadores sobre los mismos 11 parámetros** — `PRITCHARD11` (cúbicas analíticas; suficiente para vóxel y STEP) y `BSPLINE_G2` opcional (B-splines tipo ParaBlade para la vía CFD). El contrato no cambia: `profile_params` idéntico, solo cambia `profile`.

---

## 2. Del 2D al 3D

### 2.1 Apilado (stacking)

| Concepto | Efecto documentado |
|---|---|
| **Straight** (radial, por **CG**) | baseline; CG minimiza el momento flector centrífugo [M] |
| **Lean** | gradiente de presión spanwise; empuja la secundaria hacia el endwall |
| **Compound lean / bow** | *"reduce mid-span losses and forces the secondary flow region closer to the end wall"* [B] |
| **Sweep** | *"more rear-loaded at the hub and fore-loaded at the casing, reducing secondary penetration at the hub"* [B] |

Advertencia a declarar [B]: *"lo que mejora midspan empeora endwall y viceversa"*; receta industrial: **midspan recto + compound lean en ambos endwalls** ([AD Technology](https://blog.adtechnology.com/axial-turbine-stacking-best-practices-secondary-flow-suppression)). Académico: Bagshaw et al. (2008), DOI 10.1243/09576509JPE477 (**reverse compound lean prepara la fila para el endwall perfilado**).

**Recomendación**: apilado por **centroide** (reutilizando `polygon_section_props` de Phy-AC [V]) + ley de lean de 3 parámetros (`lean_hub_deg`, `lean_tip_deg`, `lean_mid_frac`) + `sweep_deg` opcional, **emitidos al contrato como línea de apilado ya evaluada** (polilínea de desplazamientos por estación) — la 5c sigue siendo "tonta".

### 2.2 Twist

Vórtice libre exacto como punto de partida [M]; con HTR < 0.6 → reacción negativa en cubo → **vórtice controlado** (Phy-AC ya lleva `vortex_n` [V]: reutilizar el mecanismo, recalibrar el rango). Restricción derivada: **reacción de hub ≥ 0.05–0.10** [M] como g. La garganta **varía con el radio**: `throat_mm` y `o_over_s` **por sección**.

### 2.3 Estaciones de span

Phy-AC: 13 (menos → pliegues visibles en el loft [V]). Turbinas exigen más (más twist; plataforma y shroud). **Recomendación**: 5 secciones de control (leyes B-spline) para el optimizador; emisión de **17 estaciones (HPT) / 21 (LPT)** con densificación coseno en hub/tip; **≥100 puntos por sección** con clustering en LE/garganta/TE (el círculo de TE necesita 12–20 puntos propios); mantener la invariante de conteo idéntico y winding consistente [V].

### 2.4 Endwall contouring

[M — pendiente de verificación]: **axisimétrico** (r_hub(z) — Phy-AC ya tiene polilíneas meridionales [V]) vs **no axisimétrico** (r(z,θ) — envolvente axial × armónicos de Fourier en θ; ref. canónica Harvey, Rose, Taylor, Shahpar, Hartland & Gregory-Smith (2000), J. Turbomach., Parts I/II [M]). **v1: solo axisimétrico + reservar el bloque `endwall.non_axisym` con `enabled: false`** (el hueco evita subir la versión MAYOR después).

### 2.5 Filetes, plataformas, tip shrouds

- **Filetes**: mantener el parámetro único (`blade_fillet_r_mm`, fuente única física+impresión como Phy-AC [V]) + opcionales LE/TE; **subir default a 3–4 mm** (creep + gradiente térmico).
- **Plataformas**: forman el endwall, sellan la cavidad de rim y alojan dampers [M]; la de Phy-AC (sigue la línea de cubo [V]) es reutilizable + overhangs (*angel wings*) + ranuras de damper opcionales.
- **Tip shroud con laberinto** [B — patentes GE US7762779B2, US10513934B2]: plataforma, **Z-notches** (enclavamiento entre álabes adyacentes, reduce tensión), **rieles de sello** (1–3), cavidades. Parametrización mínima: `enabled, platform_thickness_mm, overhang_ax_mm, overhang_tang_frac, n_fins, fin_height_mm, fin_thickness_mm, fin_z_frac[], znotch_enabled, znotch_angle_deg, honeycomb_depth_mm, clearance_mm`.
- **Consecuencia estructural**: σ_root += m_shroud·ω²·r_tip/A_root — puede ser **20–40 % de la tensión de raíz** [M]. Va en el g (§4.6).

---

## 3. La turbina completa como piezas reales

### 3.1 Diferencias mecánicas vs el compresor de Phy-AC

| # | Aspecto | Compresor | Turbina |
|---|---|---|---|
| 1 | Orden de filas | rotor→estátor | **NGV→rotor** (la primera fila es estacionaria) |
| 2 | Fijación | dovetail 1 lóbulo | **fir-tree multi-lóbulo** |
| 3 | Retención axial | ranura+placa | placa de bloqueo + **cover plates** |
| 4 | Estátores | vanos en anillos | **anillos de tobera SEGMENTADOS** con bandas |
| 5 | Punta | holgura libre | libre **o shroud** + laberinto + panal |
| 6 | Aire secundario | sangrado hacia fuera | **refrigeración hacia dentro**: purgas, alimentación, sellos |
| 7 | Térmico | moderado | **severo**; creep dominante |

### 3.2 Fir-tree de turbina vs dovetail

[B — ScienceDirect topics; PMC10384575]: el compresor usa **dovetail de un lóbulo**; la turbina **fir-tree de pares de lóbulos** — más área de cizalladura y **redundancia si agrieta un lóbulo interior**; ranuras por **brochado axial**. Referencia analítica: *Fir tree fastening — I: Deflection analysis*, Int. J. Mech. Sci. (1982) [B]. Cambios sobre `FirTree.cs` de Phy-AC [V]: (1) 2–4 pares de lóbulos + **g de anchura mínima del poste del disco** (sección crítica); (2) flancos planos ~30–45° [M]; (3) el brochado inclinado de Phy-AC se reutiliza; (4) **cuello/shank** (concentrador + conductos); (5) **agujeros radiales de alimentación de refrigerante**; (6) **tensión de apoyo en flancos** como comprobación nueva.

### 3.3 NGV segmentados

[B — patentes US7249928, EP2886800A1]: segmentos (singletes/dobletes/tripletes) con **bandas interior y exterior arqueadas**, cavidad de purga, y **strip/feather seals** en ranuras entre segmentos. Razón: un anillo continuo a 1500 K rompe por dilatación restringida [M]. Bloque `nozzle_ring` por etapa: `n_segments, vanes_per_segment, inner/outer_band_t_mm, band_overhang_ax_mm, feather_seal_slot_w/d_mm, thermal_gap_mm, mount_type`.

### 3.4 Sellos inter-etapa y cavidades de purga

[B — US8740554B2; Aerospace 6(5):60]: purga hacia la **rim cavity** para impedir ingestión de gas caliente; **knife-edge seals** hacia la cavidad aguas abajo; **cover plate con sello inter-etapa**. Bloques `interstage_seal` (laberinto: n_knives, knife_h/t/pitch_mm, r_seal_mm, clearance_mm, honeycomb_depth_mm) y `rim_cavity` (axial_gap_mm, overlap_mm angel wing, purge_slot_n/w, purge_frac).

### 3.5 Piezas del constructor 5c de Phy-AT

| Pieza | Análogo Phy-AC | Estado |
|---|---|---|
| `Shaft` | Shaft | Reutilizar (+curvic opcional) |
| `TurbineDisc{i}` | RotorStage{i} | Reutilizar voxDisc + FirTree adaptado |
| `Blade{i}` (suelto, opcional) | — | **Nuevo** — demuestra el ajuste del abeto |
| `NozzleSegment{i}` | StatorRing{i} | **Nuevo** (partición circunferencial) |
| `CasingRing{i}` | StatorRing (casco) | Reutilizar revolucionado + bridas |
| `CoverPlate{i}` | — | **Nuevo** |
| `InterstageSeal{i}` | — | **Nuevo** |
| `TipShroudRing` / panal | — | **Nuevo** (si shroud) |
| `TieBolt`/`CurvicCoupling` | tirantes | Reutilizar |
| Vistas + bom.csv + ParityReport | ídem | Reutilizar |

### 3.6 Impresión 3D de demostradores — lo que un sistema honesto declara

1. **NO son piezas calientes**: el real es monocristal/DS colado con núcleo cerámico, film cooling por EDM/láser y TBC — nada de eso se representa. Un L-PBF IN-718 es policristalino con porosidad y **no tiene la resistencia a fluencia** del material de los márgenes.
2. **Sin refrigeración interna**: la pieza es maciza; la T_metal calculada asume una refrigeración que la pieza impresa no tiene.
3. **Demuestra**: forma, ajuste, montaje, interferencias, BOM. **No demuestra**: vida, integridad, térmica, aerodinámica real.
4. Imprimir en polímero/acero para verificación dimensional; si IN-718, declarar "no apto para rotación ni exposición térmica".
5. **Trazabilidad de vóxel**: cuchillas de laberinto (0.2–0.4 mm) y ranuras de feather seal — **la 5c debe RECHAZAR (no clampear en silencio) rasgos de sello < 2 vóxeles**, o marcarlos `feature_suppressed` en el ParityReport.

**Ventaja inesperada**: los perfiles de turbina son gruesos con TE romo → el modo de fallo "espesor < 2 vóxeles" de Phy-AC prácticamente desaparece; **el loft sólido será el camino normal**.

---

## 4. Núcleo estructural de turbina (capa 1s)

### 4.1 Materiales calientes

Biblioteca propuesta (**[M] salvo indicado — verificar contra MMPDS/ASM antes de codificar**):

| Material | Forma | Uso | T capacidad |
|---|---|---|---|
| IN-718 | forja | **discos**, cover plates, carcasas | ~925 K (ya en Phy-AC [V]) |
| Waspaloy / IN-720 | forja | discos más calientes | ~1000 K |
| IN-713C/LC | colada | álabes de baja exigencia | ~1100 K |
| MAR-M247 | colada/DS | álabes y NGV | ~1250 K — *"creep 10× mayor que IN-713C a 980 °C/150 MPa"* [B — Metals 11(1):152] |
| CMSX-4 | monocristal | álabes HPT | ~1350 K |
| René N5/80 | SX/colada | álabes GE | ~1350/1200 K |
| MCrAlY + 7-8YSZ (TBC) | recubrimiento | sobre los anteriores | ΔT adicional |

**Cambios estructurales**: añadir por material los **coeficientes de Larson-Miller** y T_max de servicio; **material por PIEZA** (disco/álabe/NGV/carcasa), no uno global; **derating de E**, no solo de σ.

### 4.2 Creep — la restricción dominante

**Larson-Miller** [B — austenite.org]:

$$P_{LM} = T\,(C + \log_{10} t_r)\qquad (T\ \text{en K},\ t_r\ \text{en h},\ C\approx20\ \text{para Ni [M]})$$

*"Colapsa miles de ensayos en una curva maestra; el caballo de batalla del diseño a alta temperatura"*; *"creep limitado (<1 %) deseado para álabes"* [B]; LM más robusto frente a dispersión experimental que Norton-Bailey u Omega [B].

**Implementación L0s (analítica, <1 ms)**:
1. Curva maestra por material: log₁₀σ_rupt = a − b·P_LM + c·P_LM² (recopilación documentada en `validation/`, como Phy-AC valida el disco contra Timoshenko).
2. P_LM = T_m(C + log₁₀ t_life_h) → σ_rupt_allow.
3. **g_creep = SF·σ_root(a velocidad de DISEÑO)/σ_rupt − 1 ≤ 0**, SF≈1.25 [M]. **Sin sobrevelocidad** (carga de larga duración) — distinto del criterio de fluencia de Phy-AC.
4. Segunda g de creep para el **rim del disco**.
5. Métrica: **margen de temperatura** (K que se puede subir T_m manteniendo g≤0) — el número que un diseñador de turbinas mira primero.

**Vidas típicas [M — verificar]**: industrial 25k–50k h; civil hot section 5k–20k h; militar 1k–3k h; demostrador 200–1k h. **`t_life_h` es variable de diseño explícita del contrato.**

### 4.3 AN² para turbinas

σ_root = (π/1800)·ρ_b·A·N²·K_taper (K_taper 0.5–0.7 [M]; Phy-AC usa 0.55 [V]). Límites típicos [M]: **3–5×10¹⁰ in²·rpm²** (equiaxial/forjado), ~6×10¹⁰ (DS/SX). Contexto verificado: programa NETL "High AN² Last Stage Blade" [B]. **Diferencia clave: AN²_max debe depender de T_metal y t_life → derivarlo de la propia curva de creep** en vez de tabularlo — mejora conceptual sobre Phy-AC.

### 4.4 Temperatura del gas vs temperatura del metal

El `T_METAL_FRAC = 1.00` de Phy-AC [V] sería **catastrófico** en turbina (T₀=1700 K → siempre infactible). Cadena térmica correcta:

1. **T_gas local con perfil de combustor**: T_gas(r) = T₀₄_mean + RTDF(r)·(T₀₄−T₀₃); pico +OTDF·(T₀₄−T₀₃).
2. **Rotor: temperatura total RELATIVA** — T₀_rel = T₀ − U·C_u/c_p (decenas a >100 K menos que el NGV; **debe estar en el modelo**).
3. Temperatura de recuperación ≈ T₀_rel × factor ~0.89 [M].
4. **Efectividad de refrigeración** [B — NETL handbook 4.2.2.1]:
   $$\eta_c = \frac{T_{gas} - T_{metal}}{T_{gas} - T_{coolant}}\ \Rightarrow\ T_{metal} = T_{gas} - \eta_c(T_{gas}-T_{coolant});\qquad \eta_c \approx 0.50\ \text{"conventional"}$$
   Niveles: 0 / 0.30 (convección) / 0.50 (convección+película) / 0.60–0.70 (avanzada) [B/M].
5. **TBC**: resta adicional (§4.7).

**Contrapartida termodinámica**: el refrigerante se descuenta del gasto que trabaja y se reinyecta con pérdida de mezcla — el contrato lleva `cooling.frac_per_row` y el meanline lo usa.

### 4.5 Discos: burst **con gradiente térmico** y sobrevelocidad regulatoria

`solve_disk` de Phy-AC (Thomas O(n), validado vs Timoshenko [V]) se **extiende con el término termoelástico αE·dT/dr** — misma tridiagonal, ~20 líneas; **test de regresión: con ΔT=0 reproduce Timoshenko bit a bit**. ΔT barreno→rim de 200–400 K [M]; perfil T(r) = T_bore + (T_rim−T_bore)·((r−r_b)/(r_rim−r_b))^p, p≈2 [M].

**Sobrevelocidad**: Phy-AC usa 1.05 [V]. Cifra regulatoria [B]:
> **14 CFR 33.27**: el rotor no debe estallar operado **5 minutos al 120 %** de la velocidad máxima permisible ([eCFR](https://www.ecfr.gov/current/title-14/chapter-I/subchapter-C/part-33/subpart-B/section-33.27)); FAA armonizó 115→**120 %** con EASA CS-E ([Federal Register 2011-18002](https://www.federalregister.gov/documents/2011/07/18/2011-18002/airworthiness-standards-rotor-overspeed-requirements)); AC 33.27-1A.

**Phy-AT debe usar `OVERSPEED_BURST = 1.20` con la cita en el código.**

### 4.6 Campbell / HCF en turbinas

1. **La excitación dominante es el conteo de NGV aguas arriba** + puntales + **inyectores del combustor** (el OTDF circunferencial ES una excitación de orden n_injectors) + fila anterior.
2. **Álabes con shroud NO son vigas empotradas-libres**: el Z-notch acopla circunferencialmente (diámetros nodales); modelarlo empotrado-libre **subestima gravemente** f₁. Con shroud: λ₁² de empotrado-**apoyado** (≈15.42 vs 3.516) [M], declarado como aproximación de primer orden.
3. **E(T) derrateado** (IN-718 cae ~20 % de 293 a 923 K [M]) — Phy-AC no lo hace [V].
4. Dampers y fricción de Z-notch no modelados → declarar que el margen de Campbell es la única defensa modelada.

### 4.7 TBC

[B — NETL 4.4.2]: bond coat MCrAlY/aluminide (~0.13 mm) + top coat **7-8YSZ ~0.25 mm**; *"YSZ redujo la T máxima del álabe un 18 %"* [B]. **Advertencia**: circula un dato incompatible (1–5 mm, ΔT 373–573 K) — **no usarlo**. Modelo L0s: ΔT_TBC = q″·t/k con k_YSZ ≈ 1.0–1.7 W/mK [M] → **~80–150 K** para 0.25 mm. El TBC añade masa en punta, puede desprenderse, y **el demostrador impreso no lo lleva** — declararlo.

---

## 5. Software y librerías para la vía STEP/STL

### 5.1 CadQuery/OCCT

La vía CadQuery de Phy-AC (`_cq_blade_trimmed`, `_cq_firtree_solid`, `_cq_disc`, `assembly_interferences`, `export_step` [V]) **funciona mejor con TE grueso** que con el afilado del compresor (el punto singular del TE afilado es donde OCCT degenera). Dos precauciones: conteo idéntico de puntos con correspondencia consistente (ya invariante [V]); mayor twist → 17–21 estaciones (§2.3).

### 5.2 PicoGK / LEAP 71

[B — leap71.com]: PicoGK open-source 2023; **Noyron** = Large Computational Engineering Model propietario; CEMs especializados: **RP (cohetes), EA (actuación), HX (intercambiadores)**. Hito: **hot-fire de un motor cohete diseñado por Noyron (jun 2024)**. **Nota de honestidad: NO hay evidencia de un "turbofan paramétrico" de LEAP 71** — no construir argumentos sobre esa premisa. Lecciones: (1) **el modelo ES el código**; (2) kernel implícito/vóxel para robustez booleana; (3) la familia de CEMs valida Phy-AC/AT/CB; (4) **la validación es física, no gráfica** (paridad STL↔STEP, interferencias, bancos contra soluciones exactas — lo que Phy-AC ya hace [V]).

### 5.3 Generadores open-source (inventario)

| Proyecto | Turbinas | Salida | Utilidad |
|---|---|---|---|
| [nasa/pyturbo-aero](https://github.com/nasa/pyturbo-aero) [V] | **Sí** | STL, JSON | ★★★★★ referencia de perfil 2D y apilado |
| [ParaBlade](https://github.com/NAnand-TUD/parablade) [V] | **Sí** | NURBS + sensibilidades | ★★★★ continuidad G2, canal meridional |
| [T-Blade3](https://github.com/GTSL-UC/T-Blade3) [V] | Sí | .dat + UDP ESP | ★★★ comba por curvatura |
| [11-Parameters-Turbine-Blade-Generator](https://github.com/DavidPoves/11-Parameters-Turbine-Blade-Generator) [V] | Sí (Pritchard) | puntos 2D | ★★★★★ implementación legible del método elegido |
| [OpenOrion/turbodesigner](https://github.com/OpenOrion/turbodesigner) [V] | **NO** — "plans to support axial turbines; no turbine implementation exists yet" | STEP | ★★★ el ancestro de la 5a de Phy-AC; **la turbina hay que escribirla** |

### 5.4 Formatos

[B — Capvidia/Datakit]: AP203 (B-rep), AP214 (+GD&T/color), **AP242** (fusiona ambos + PMI semántico + teselación). Recomendación: **AP242 si OCCT lo soporta; AP214 fallback** (lo que Phy-AC emite hoy [V]). Emitir STL + STEP y **medir paridad de volumen** (reutilizar `parity_stl_step.py` + `ParityReport.cs` [V]). Declarar la diferencia sistemática (redondeo de ½ vóxel).

### 5.5 Requisitos de malladores CFD

**TurboGrid** [B — User's Guide]: necesita `hub.curve`, `shroud.curve`, `profile.curve` (secciones (x,y,z) por span), holgura de punta como % de span o distancia. Un exportador `contrato → {hub,shroud,profile}.curve` es inmediato si el contrato lleva polilíneas meridionales, secciones 3D por span, n_blades y holgura — Phy-AC ya tiene casi todo [V]. Phy-AT añade **la garganta** (para refinamiento) y **shrouded/unshrouded** (topología de malla distinta). AutoGrid5: mismo conjunto de datos [M].

---

## 6. El contrato `phyat-axial-1`

### 6.1 Principios heredados de Phy-AC [V]

Frontera única 5a↔consumidores; **versión MAYOR en el identificador** (el consumidor rechaza lo que no entienda); **validador propio sin dependencias** ("una comprobación que se salta sola cuando falta una dependencia opcional no es una comprobación"); `nSCHEMA_MAJOR` en C# sube a la vez y la suite lo comprueba.

### 6.2 Los cuatro consumidores

(a) constructor vóxel/STL; (b) vía CAD/STEP; (c) mallado CFD (perfiles 3D por span, garganta, holguras, RPM, interfaces); **(d) Phy-CB aguas arriba** — interfaz declarada en ambos sentidos: plano, radios, **perfil radial T₀ (RTDF), OTDF y nº de inyectores** (= orden de excitación del NGV), swirl α(r), P₀(r), fracción/temperatura del aire de refrigeración reclamado, presupuesto de caída de presión.

### 6.3 Bloques propuestos (esqueleto)

```
schema: "phyat-axial-1"
design_vector: { n_stages, RPM, HTR_in, phi1, psi_mid, psi_slope, Rx_mean, Rx_slope,
                 Zw_nozzle, Zw_rotor,      ← Zweifel SUSTITUYE a sigma_r/sigma_s
                 AR, T04, P03, massflow,
                 cool_frac_total, t_life_h }  ← NUEVOS
derived: { ER, eta_poly/isen, power_W, U_tip, AN2, length_mm, r_tip_in_mm,
           M_abs_max, M_rel_max, M_exit_throat_max, reaction_hub_min,
           vortex_n, gas_model, ds_machine, cooling{per_row,total_frac,eta_c,mdot_exit},
           choked_rows[], feasible, source }
annulus: { hub[[z,r]], tip[[z,r]], per_row_clearance_mm[] }   ← por fila (shrouded ≠ un-)
endwall: { axisym_contour_enabled, non_axisym{enabled:false, harmonics[], axial_envelope[]} }
stages[]: { index, nozzle:<row>, rotor:<row>, interstage_seal:<seal>, rim_cavity:<cavity> }
           ← ORDEN: nozzle PRIMERO
exit_guide_vane: <row>|null
assembly: { shaft*, disc{...}, firtree{n_lobes, flank_angle_deg, neck_w_mm,
            min_post_width_mm, coolant_feed_d_mm/n, ...}, cover_plate{...},
            tie_bolt*/curvic, flange*, blade_fillet_r_mm(+le/te),
            platform{t, overhangs, damper_slot}, tip_shroud{...}, nozzle_ring[]{...} }
thermal:   ← BLOQUE NUEVO DE PRIMER NIVEL (lo consumen estructural, CFD y refrigeración)
  { T03_K, T04_mean_K, OTDF, RTDF, n_injectors,        ← vienen de Phy-CB
    radial_profile[[span,T0/T0mean]],
    per_row[{ T0_abs, T0_rel, T_recovery, eta_cooling, T_coolant, T_metal,
              tbc{enabled,thickness_mm,dT_K} }] }
structural: { materials{disc,blade,nozzle,casing,shroud},   ← por PIEZA
    per_stage[{ T_metal_blade, T_bore, T_rim, sigma_vm_max, sigma_root, k_t,
                P_LM, sigma_rupt_allow, creep_ratio, yield_ratio, burst_margin, AN2 }],
    campbell{}, flutter{}, g_struct[], feasible_struct,
    overspeed_burst: 1.20, cert_ref: "14 CFR 33.27", t_life_h, creep_SF }
interfaces:   ← EL BLOQUE DE REUTILIZACIÓN
  { upstream:   { system:"phycb-annular-1"|null, plane_z_mm, r_hub/tip_mm,
                  requires:[T0_radial_profile, OTDF, RTDF, swirl_alpha_profile,
                            P0_profile, n_injectors, coolant_available_frac] },
    downstream: { system:"phyat-lpt-1"|"exhaust", provides:[T0_profile, alpha_deg, M] },
    coolant_source: { from:"phyac-axial-2", bleed_stage, P0, T0, frac } }
cfd: { boundary_conditions{ P0/T0_in (+perfiles), alpha_in_profile, mdot, RPM,
        P_static_out, Tu, interfaces[] },
       mesh_hints{ per_row[{ n_blades, rotating, shrouded, throat_mm_by_span,
                             o_over_s_by_span, tip_clearance_mm, cavity_present }] } }
manufacturing:   ← el bloque de HONESTIDAD
  { target:"demonstrator_print", min_feature_mm, recommended_voxel_mm,
    not_represented:[internal_cooling_passages, film_holes, TBC,
                     single_crystal_microstructure, damper_pins],
    suppressed_features[], warning:"Geometry/assembly demonstrator only. Not a
    hot-section part; not for rotation or thermal exposure." }
provenance: { profile_method:"pritchard-11", exit_angle_rule, loss_model,
              deviation_model, creep_model:"larson-miller", LM_constant_C:20, citations[] }
run_meta: {...}
```

**`<row>`**: row_id, kind (nozzle|rotor|egv), n_blades, z_le/te/center, rotating, shrouded, zweifel, pitch_mid_mm, stacking{type, lean_deg_by_span[], sweep_mm_by_span[]}, sections[].

**`<section>`** (el bloque más importante): span_frac, r_mm, cuerdas, stagger_deg, metal_in/out_deg, cuñas, r_le/te_mm, **throat_mm, pitch_mm, o_over_s, gauge_angle_deg, unguided_turning_deg**, t_max/c, x_tmax/c, M_in, M_out_is, reaction_local, sigma_local, profile (PRITCHARD11|BSPLINE_G2|CD_SUPERSONIC), **profile_params (los 11 de Pritchard verbatim)**, points[] (cerrado CCW, N idéntico por fila), camber_points[], **index_markers{i_le, i_te_start/end, i_throat_ss, i_ss_end}**.

**Cinco decisiones justificadas**: (1) `profile_params` ADEMÁS de `points` — cualquier consumidor puede re-renderizar a otra resolución o en G2 sin correr la 5a; (2) `index_markers` — un mallador necesita saber qué es dorso/vientre/garganta; (3) garganta **por sección** (el twist la hace variar); (4) `thermal` de primer nivel (lo consumen tres partes); (5) `interfaces` con `requires`/`provides` — el mecanismo por el que Phy-CB sabrá exactamente qué entregar, generalización del `derived.bleed` de Phy-AC [V].

---

## 7. Recomendaciones concretas (capas 5a/5c/1s)

### 7.1 Perfil: **Pritchard 11-par con renderizador dual** (§1.8)

Correlaciones de los 6 defaults:

| Default | Correlación | Verificación |
|---|---|---|
| Paso s | Zweifel = 0.8 (rango 0.8–1.1) | [M]/[B] |
| Garganta o | s·cosα₂ corregido por curvatura del dorso y Mach (Aungier) | [B] |
| δ_ug | default 8°, **duro ≤20°** (ancla: 19.47° NASA) | [B]/[M] |
| ε_i | 5–15° | [M] |
| r_LE | 3–8 % de cuerda | [M] |
| r_TE | t_TE/o ≈ 0.05 (0.08–0.12 con eyección) | [B] |
| ε_o | **no libre** — cierre de garganta | [V] |

g geométricas nuevas: δ_ug ≤ 20° · o/s ∈ [0.15, 0.85] · t_TE ≤ 0.15·o · reaction_hub ≥ 0 · M_out_is ≤ 1.15 (o CD_SUPERSONIC con supported:false) · no autointersección (invariante analítica, como Phy-AC [V]).

Emisión: 17/21 secciones, ≥100 puntos con clustering LE/garganta/TE, conteo idéntico, stagger firmado, camber + index_markers.

### 7.2 Constructor C#: reutilizable ~60 %

**Sin cambios** [V]: PhyATImport (rechazo duro de versiones), BladeRow (loft sólido — funciona MEJOR con perfiles gruesos), recorte contra la vena, voxRevolveZ, LocalFrame, voxDisc (rim que sigue el cubo, bolt circle común, casco que para en la banda), Casing con bridas en ambos extremos, ParityReport + parity_stl_step, convención de holguras + warnings de vóxel, bom.csv.

**Nuevo**: TurbineFirTree (★★★★★), NozzleSegment (★★★★★), TipShroud (★★★★), InterstageSeal (★★★), CoverPlate (★★★), RimCavity (★★), CoolingFeed (★★), **ManufacturabilityCheck** (rechaza o marca `feature_suppressed`, nunca clampea en silencio) (★★★★).

**Orden de piezas**: la primera pieza del tren gaseoso es un NozzleSegment; los planos de corte se desplazan media etapa.

### 7.3 El g estructural de turbina (9 componentes duros)

| # | Restricción | Anclaje |
|---|---|---|
| 1 | Fluencia del disco a sobrevelocidad: σ_vm·(1.20)²/σ_y(T) − 1 ≤ 0 | solver 1-D + término térmico; 14 CFR 33.27 |
| 2 | Burst ≥ 1.22 a 1.20·N_max | método Phy-AC + 33.27 |
| 3 | AN² ≤ AN²_max(material, T_metal, t_life) **derivado de la curva de creep** | §4.3 |
| 4 | Raíz estática con K_t **+ σ_shroud**: K_t(σ_root+σ_sh)(1.05)²/σ_y(T_m) ≤ 1 | Peterson [V] + shroud [M] |
| 5 | **CREEP del álabe** (dominante): 1.25(σ_root+σ_sh)/σ_rupt(P_LM) ≤ 1, **a velocidad de DISEÑO** | Larson-Miller |
| 6 | **CREEP del rim** con T_rim | ídem |
| 7 | **T_metal ≤ T_max(material)** con T_metal = T_gas,rel − η_c(T_gas,rel−T_cool) − ΔT_TBC | η_c≈0.50 [B]; ΔT_TBC 80–150 K |
| 8 | Campbell vs EO = nº NGV **y** EO = nº inyectores; λ₁ empotrado-apoyado si shroud; E(T) derrateado | Southwell + [M] |
| 9 | **Poste del disco / apoyo del abeto**: σ_bearing/σ_y(T_rim) ≤ 1 y w_post ≥ w_min | [B] |

Métricas (no g): margen de temperatura; Goodman con relajación por creep; flutter V*≤1.4; EO bajos k=1..6; **fracción de vida a la T de pico del OTDF** (el hot streak mata álabes); LCF del barreno (declarado NO modelado en L0s).

**Infraestructura**: solve_disk + αE·dT/dr (test: ΔT=0 → Timoshenko exacto); MATERIALS + T_max_K, lm_C, lm_curve, derating de E; materiales por pieza; AN²_max calculado; OVERSPEED_BURST=1.20 citado; N_STRUCT=9.

### 7.4 Riesgos declarados / verificación pendiente

1. El PDF de Pritchard no se pudo leer — los 11 parámetros verificados contra implementación de terceros [V], **conseguir el paper** (especialmente los defaults).
2. Correlaciones de ángulo de salida: cualitativas [B], **ningún coeficiente verificado aquí** (el informe 01 los verifica vía TurboFlow).
3. Datos de Larson-Miller por material: recopilar de fuentes primarias y documentar en `validation/`.
4. Límites de AN²: derivarlos del creep (elimina el problema de raíz).
5. ΔT del TBC: contradicción de un orden de magnitud entre fuentes — documentar y calibrar.
6. Endwall no axisimétrico y gradiente térmico de disco: [M] con referencias nominales (Harvey et al. 2000) a verificar.
7. **No hay evidencia del "turbofan paramétrico" de LEAP 71** (lo documentado: Noyron RP/EA/HX).
8. **Ninguna vía open-source cubre Phy-AT end-to-end**: turbodesigner no tiene turbinas [V]; pyturbo-aero no tiene STEP/fillets/shroud [V]; parablade no tiene ensamblaje. **La capa 5c de Phy-AC es el activo más valioso del proyecto — reutilizable en ~60 %.**

---

## Fuentes

(hipervínculos in-line; principales: Pritchard 85-GT-219 + implementación DavidPoves; Agromayor et al. CAD 133 (2021) + parablade; nasa/pyturbo-aero (NASA-1.3); GTSL-UC/T-Blade3; Kulfan 2007; OpenOrion/turbodesigner; NASA TE study 19720007336; Flegel 20140012427; Bagshaw 2008; patentes de tip shroud US7762779B2/US10513934B2, NGV EP2886800A1, interstage US8740554B2; fir-tree PMC10384575; Larson-Miller austenite.org; Metals 11(1):152; NETL handbook 4.2.2.1 y 4.4.2; 14 CFR 33.27 + Federal Register 2011-18002 + AC 33.27-1A + EASA SC; NETL high-AN² LSB; TurboGrid User's Guide; LEAP 71 PicoGK/Noyron; STEP AP203/214/242 Capvidia/Datakit.)
