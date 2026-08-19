<!--
Quasar Phy-AT · Investigación de recursos — turbinas axiales
Documento de investigación previo al arranque del proyecto Phy-AT.
Síntesis de cinco informes de investigación (docs/phyat_research/).
-->

# 🔬 Quasar Phy-AT — Investigación de recursos para el análogo de Phy-AC en turbinas axiales

**Estado**: investigación previa completada (no hay código Phy-AT todavía).
**Objetivo**: reunir la literatura, los métodos, los sistemas y los principios
necesarios para construir un sistema autónomo de diseño inverso de **turbinas
axiales** con la lógica de Phy-AC: *prior físico calibrado → ensemble profundo
residual con puerta de incertidumbre → NSGA-II restringido → informe
autocontenido → contrato de geometría → STEP/STL imprimibles*, con **escalera
de fidelidades** (meanline → SCM → CFD/datos) e **incertidumbre calibrada**.

**Informes completos** (este documento es la síntesis; el detalle, las
ecuaciones y las citas completas viven en los cinco informes):

| Informe | Contenido |
|---|---|
| [01 · Meanline y pérdidas](phyat_research/01_meanline_perdidas_turbinas.md) | Triángulos, Smith, cadena AMDC-KO/Benner con ecuaciones verificadas, Zweifel, off-design/Stodola/choking, refrigeración meanline (Gauntner/Hartsel/Young-Wilcock), campaña de validación NASA, θ y g propuestos |
| [02 · Through-flow / SCM](phyat_research/02_scm_throughflow_turbinas.md) | Equilibrio radial con rotalpía, choking y doble solución, cierre por garganta (o/s), reparto spanwise de pérdidas (Benner), refrigeración en tres lazos, códigos existentes, plan de reutilización de `scm_core.py`, modos de fallo honesto |
| [03 · Estado del arte CFD / fidelidades altas](phyat_research/03_sota_cfd_altas_fidelidades.md) | Cadena industrial, RANS/transición, GPU, efectos que las correlaciones no ven (purga, hot streaks, shroud, wake recovery), benchmarks públicos, calibración KOH, precisiones honestas por fidelidad, paquete de BCs CFD |
| [04 · ML, PINNs, UQ y optimización](phyat_research/04_ml_pinns_uq_optimizacion.md) | Veredicto sobre PINNs/operadores neuronales/generativos, multifidelidad KOH, conformal prediction, qué conservar de Phy-AC y 10 mejoras puntuales priorizadas |
| [05 · Geometría y estructuras](phyat_research/05_geometria_estructuras_turbina.md) | Parametrización de Pritchard, apilado/lean/shroud, arquitectura mecánica (fir-tree, NGV segmentados, sellos), creep/Larson-Miller, 14 CFR 33.27, contrato `phyat-axial-1`, plan de la capa 5c |
| [06 · GPU y CFD diferenciable](phyat_research/06_gpu_diferenciable_xlb_jaxfluids_warp.md) | XLB, JAX-Fluids y NVIDIA Warp auditados **leyendo su código fuente**: veredictos (descartar/vigilar), el ecosistema (SU2, PyFR, JAX-FVM, FluidX3D…), y dónde sí capturar el valor diferenciable (JAX sobre L0/L1, MULTALL como L3 interno) |

> **Nota de método**: los cinco informes fueron elaborados con acceso web
> restringido (proxy de egreso que bloquea los dominios académicos). Cada
> afirmación lleva marcador de fiabilidad ([V] verificado sobre fuente
> primaria o código de referencia, [S]/[B] de resúmenes indexados, [M]/[?] a
> verificar). **Las ecuaciones centrales (pérdidas, desviación, choking,
> penetración de Benner) están verificadas línea a línea contra TurboFlow y
> AxialOpt (MIT)**; los números marcados deben confirmarse contra los PDF
> antes de entrar en `Quasar_PhyAT_Science.md`.

---

## A. El patrón Phy-AC: qué se replica y qué cambia con una turbina

La arquitectura de Phy-AC (README.md, docs/Quasar_PhyAC_Science.md) se
traslada capa a capa. Lo que cambia no es la topología del sistema sino la
física de cada capa:

| Capa | Phy-AC (compresor) | Phy-AT (turbina) | Veredicto |
|---|---|---|---|
| **L0 físico** | stage-stacking con Lieblein/Howell/Koch & Smith; restricción dominante = **bombeo** (Koch 1981) | stage-stacking con **AMDC-KO + Benner**; restricción dominante = **choking** (la primera tobera fija el caudal — Stodola) | misma marcha de entropía con gas imperfecto (fase 9); pérdidas y restricciones **nuevas** |
| **Gas** | aire, cp(T) 250–1000 K | productos de combustión: **cp(T, FAR), R(FAR) hasta ~1900 K**; el gasto y la composición cambian por fila (refrigeración) | extender `physics_core` |
| **L1 SCM** | `scm_core.py`: equilibrio radial + Carter, banda de pared constante | `scm_turbine.py` **hermano** (patrón PCA SC90T/SC90C): rotalpía en rotor, cierre por **garganta** (arccos(A_t/A_o)), banda de pared = **Z_TE/H de Benner**, choking por streamtube, refrigeración en tres lazos | ~50 % de utilidades compartidas (`scm_common.py`) |
| **L2 calibración** | afín a·y+b sobre (PR, η) | afín **ponderada y regularizada** sobre (**Γ capacidad PRIMERO**, η_tt, ER, α_exit, Λ) | misma maquinaria, claves nuevas |
| **Capas 2–4** | ensemble residual + puerta + NSGA-II/Deb + LCB + k-means | **idéntico** — el core es agnóstico de dominio; cambian el embedding (features de turbina) y los rangos de θ | + mejoras M1–M10 (§E) |
| **5a geometría** | NACA-65/DCA + Lieblein/Carter, 13 estaciones | **Pritchard 11 parámetros** (la garganta es ENTRADA), 17–21 estaciones, ≥100 pts/sección | familia de perfiles nueva; contrato nuevo |
| **5c C#/PicoGK** | eje + discos álabeados + anillos de estátor | + **NGV segmentados, tip shrouds, sellos inter-etapa, cover plates** | **~60 % reutilizable** |
| **1s estructural** | fluencia/burst/AN²/raíz/Campbell (5 g) | + **creep (Larson-Miller), T_metal con refrigeración, burst al 120 % (14 CFR 33.27), poste del abeto** (9 g) | solver de disco reutilizado + término térmico |
| **Validación** | NASA Stage 35, R37, R67, E³ HPC | **Kofskey (aire frío, geometrías ya digitalizadas) → E³ HPT** | mismo protocolo: cero recalibración por máquina |

**Las cinco inversiones conceptuales** que distinguen una turbina y que
gobiernan todo el diseño del sistema:

1. **El choque no es un fallo, es el modo de operación normal.** El estátor 1
   de una HPT está estrangulado en diseño. La penalización continua de choke
   de Phy-AC se reinterpreta: feature continua de margen de garganta + el
   sistema debe distinguir "gargantado por diseño" de "gargantado por error"
   (TurboFlow, J. Turbomach. 147(4):041002).
2. **El ángulo de salida lo fija la garganta, no una correlación de
   desviación.** cosα₂ = o/s (+ corrección de Mach). Carter desaparece; la
   garganta pasa a ser la variable de diseño primaria y una ENTRADA del
   perfil (Pritchard).
3. **El caudal lo fija la primera tobera** (ley de la elipse de Stodola):
   el mapa colapsa a casi una curva; el margen de bombeo se sustituye por el
   margen/asignación de choking; L1 puede declarar el gasto de L0 infactible
   — y eso es información, no divergencia.
4. **La temperatura es una variable de primer orden.** Perfil de combustor
   (OTDF/RTDF), temperatura relativa del rotor (T₀_rel = T₀ − U·Cu/cp),
   refrigeración (efectividad η_c ≈ 0.5 convencional), TBC (~80–150 K), y
   **creep como restricción dominante** con la vida `t_life_h` como variable
   explícita del contrato.
5. **No hay stall → sí hay incidencia.** El talón de Aquiles medido de
   Phy-AC (fase 12/F-02: −4.4 pts de η en el extremo de choke por el bucket
   de incidencia) es un problema **resuelto de forma citada** en turbinas
   (Moustapha 1990, Benner 1997). Phy-AT nace con la ley correcta.

---

## B. Fidelidad L0 — meanline de turbina (síntesis del informe 01)

**Decisión central**: cadena de pérdidas **AMDC-KO + Benner**:

$$Y_{tot} = \bigl(Y_p^{KO} + Y_{te}^{KO} + Y_{inc}^{B97}\bigr)\left(1-\tfrac{Z_{TE}}{H}\right) + Y_s^{B06} + Y_{cl}^{DC}$$

con desviación de **Aungier 2006** (C² — la dominancia de Deb exige
continuidad) y choking por **Mach crítico analítico** (sin iteración).
Razones: forma algebraica cerrada (0.5 ms/punto alcanzable en NumPy),
**±1.5 % en η sobre 33 turbinas sin recalibración** (Kacker & Okapuu 1982),
validación abierta reproducible contra 3 turbinas NASA (TurboFlow: caudal
100 % de puntos <2.5 %), y la incidencia de Benner de fábrica. **Craig-Cox y
Traupel quedan descartados como base** (sistemas de ~18 cartas digitalizadas
— verificado en el código de NASA TD3) y se mantienen como contraste offline.

**Triángulos**: c_θ2/U = (1−R)+ψ/2, c_θ3/U = (1−R)−ψ/2 — el espejo exacto del
par de Phy-AC. Identidades de verificación: tanα−tanβ = 1/φ; ψ = φ(tanβ₂−tanβ₃).
**Trampas nº 1 y 2**: convención de ángulos (axial vs tangencial — función
única de conversión + test) y **coeficiente de pérdida referido a la dinámica
de SALIDA** (no de entrada como en compresor).

**θ propuesto (18-D)**: n_stages, RPM, HTR (suelo 0.50 — tabla f_hub de KO),
φ/ψ/R con pendientes (ψ ∈ **0.80–2.60** — copiar los rangos de Phy-AC sería
absurdo), **Zweifel por fila en lugar de solidez** (= F_t de Benner, sin
conversión), AR por fila, t_te/o, α_in, **cool_tech (x_factor de Gauntner)**,
y el punto de operación pinned al final (compatibilidad de checkpoints).
El radio de punta NO es variable (lección M1 de Phy-AC); holgura **absoluta
en mm** (confirmado por los casos NASA).

**Refrigeración meanline**: caudal por el algoritmo de **Gauntner (NASA
TM-81453)** — implementación exacta verificada en pyCycle, con el factor de
perfil PF=0.30 en la primera tobera y el rotor viendo 0.92·T₀ (temperatura
relativa) —; pérdida de mezcla de **Hartsel (1972)**; marco termodinámico de
**Young & Wilcock (2002)** (estátor y rotor por separado). **La definición de
η de una turbina refrigerada vale 2.5 puntos** (E³: 90.0 % termodinámica vs
92.5 % ciclo GE) → validar primero en aire frío.

**g(θ)** en cuatro grupos: A) validez de correlaciones (|β_out| 40–80°,
s/c 0.30–1.10, t_te/o ≤ 0.40, HTR ≥ 0.50…); B) aerodinámicas (Mach, swirl de
salida ≤15°, reacción 0.05–0.85, Zweifel 0.70–1.15, flare ≤12.5°, **margen de
choking** — el análogo del margen de bombeo); C) mecánicas (σ_ct exacta, no
AN² tabulado; T_metal; fracción de refrigerante ≤ 0.20–0.25); D) spec.
Todas continuas y finitas para θ degenerados.

---

## C. Fidelidad L1 — SCM de turbina (síntesis del informe 02)

**`scm_turbine.py` como módulo HERMANO de `scm_core.py`** (PCA mantiene
SC90T/SC90C gemelos desde 1990), con `scm_common.py` compartido
(reposicionado, curvatura, mezcla spanwise, streamtubes).

Cambios estructurales sobre el SCM de Phy-AC:

1. **Rotalpía** I = h + ½W² − ½U² en el rotor (h₀ no se conserva con
   refrigeración; el refrigerante absorbe bombeo centrífugo).
2. **Estación de GARGANTA por fila** (3n+1 estaciones): Denton (1978) — "el
   choking ocurre en la garganta"; Casey & Robinson (2010) — la garganta como
   plano virtual con **límite de gasto POR TUBO DE CORRIENTE y redistribución
   spanwise** (lo que un meanline no puede ver: el cubo de la tobera bloquea
   antes que la punta).
3. **Cierre del ángulo por áreas**: α₂(r) = arccos(A_throat/A_out) − δ(M₂)
   (Aungier), y en supersónico por continuidad con el gasto crítico. Carter
   desaparece; lo que se congela es la garganta.
4. **Doble solución sub/supersónica**: la bisección de `_solve_station_cm`
   está mal planteada cerca del bloqueo (Tiwari, Stein & Lin 2013) — acotar a
   M_m ≤ M_crit (corregido por pérdidas, fórmula cerrada) + rama supersónica
   explícita. NASA descubrió esta limitación DESPUÉS de construir OTAC
   (Hendricks 2016): **el choking es requisito de la v1, no extensión**.
5. **Reparto spanwise**: la banda de pared deja de ser `WALL_BAND_FRAC=0.30`
   y pasa a ser la **profundidad de penetración Z_TE/H de Benner** (función
   de carga, Re y AR: ~0.10 en rotor largo de LPT, >0.35 en tobera corta de
   HP). El choque se paga por línea de corriente sin el factor f_hub (el SCM
   ya conoce el Mach local). Fuga con shroud = **fuente/sumidero de masa que
   puentea la fila**, no Δs en banda.
6. **Refrigeración en tres lazos anidados** (LUAX-T): geometría / entropía /
   refrigeración. El término cinético de la entropía de mezcla depende del
   **vector** de velocidad — segunda gran razón de existir del L1.
7. **`CURV_MAX` se elimina** (en una LPT con flare recortaría física real —
   suavizar r'' en su lugar) y se añade el **test del actuator disc** (el
   único que verifica de verdad el término de curvatura).
8. **Taxonomía de fallo honesto**: `TurbineChoked` (NO es fallo),
   `MassFlowInfeasible` (contradicción L0↔geometría↔L1, etiquetada con el
   gasto máximo alcanzable), `LimitLoading`, `NegativeReactionAtHub` (fallo
   del DISEÑO, no del solver), `EnergyImbalance` (assert duro), etc.

**Verificación**: T21 de Phy-AC se reutiliza sin cambios (el vórtice general
no sabe si la máquina expande); + actuator disc; + tobera bloqueada analítica
(una línea que atrapa la mayoría de errores de choking); + balance de energía
cerrado con refrigeración.

---

## D. Fidelidades altas y estado del arte (síntesis del informe 03)

**La cadena industrial no ha cambiado de topología desde 2010; ha cambiado el
coste por escalón**: un RANS 3D multi-fila con fugas cabe en **~10 min en una
GPU A100** (Turbostream 3 / turbigen) → un "par hi-fi" de calibración cuesta
~1 h de reloj, y **25 pares semilla son un fin de semana de máquina**.

**La advertencia de Denton (GT2010-22540) es la base del lazo de
calibración**: el CFD de turbomáquinas se usa *comparativamente* — los
errores dominantes son BCs desconocidas, geometría desconocida y la hipótesis
de estacionariedad. Exactamente el argumento del data flywheel.

**Escalera propuesta**: L0 (~1 ms) → L1 SCM (~1–5 s) → **L2 = Q3D de
transición (tipo MISES, el lapso de Re de LPT no es opcional) + corrección
afín-más-residual** (~0.1–10 s) → L3 RANS externo con emisión de BCs y
retorno de pares (~10 min–2 h) → L4 (URANS/LES) **fuera del lazo, solo
priors** (CFD Vision 2030: LES no será viable en el lazo ni en 2030).

**Efectos que las correlaciones no ven** (candidatos prioritarios a
calibración, con magnitudes): purga de rim seals (**0.7–1.2 pts de η por 1 %
de gasto**), fuga de punta (~1/3 de la pérdida del escalón), shroud (hasta
1 pt por la cavidad), film cooling, **wake recovery en LPT** (⇒ L0/L1
sobreestiman la pérdida de perfiles high-lift — sesgo sistemático a
declarar), hot streaks/Kerrebrock-Mikolajczak (+10–30 % de carga térmica en
la cara de presión), clocking (0.08–0.4 % — **bajo el ruido del banco: fuera
de alcance declarado**).

**Precisiones honestas a declarar**: L0 ±1.5–2.5 pts η / ±3–5 % capacidad;
L1 ±1–2 pts; L2 calibrado ±0.3–0.5 pts dentro del envolvente; RANS ±1–2 pts
absoluto, ±0.3–0.5 comparativo; **suelo del banco: 0.45 % U95 absoluto /
0.25 % relativo** (QinetiQ TTF) — nunca declarar mejoras por debajo.
**Techo de ganancia 3D creíble: 0.5–1.0 pt de η** (endwall contouring del
Trent 500: +0.59±0.25 % medido) — codificar como aviso del optimizador.

**Calibración (KEYS)**: **Γ = ṁ√T₀/p₀ (capacidad) SIEMPRE PRIMERO** — el
aft-loading mueve la capacidad hasta un 10 %, y calibrar η sin calibrar Γ
absorbe un error de matching como pérdida —, luego η_tt, ER, α_exit, Λ.
Guardarraíles: rechazar |a−1|>0.15 (error de definición, no física); congelar
fuera de la cáscara convexa de los pares; **invariante de definición de η
refrigerada (Young & Wilcock) verificado entre L0 y CFD**.

**El paquete de BCs de turbina** es mucho más rico que el de Phy-AC:
perfiles radiales de P₀/T₀ (no escalares), swirl residual, hot streak con
clocking, **Tu + escala de longitud L_t** (dos usuarios con el mismo Tu y
distinta L_t miden pérdidas distintas), térmica de pared (adiabática /
isoterma con T_wall por superficie), refrigeración por hilera (BCs de
orificio), purgas con **swirl ratio**, mixing plane no reflectante,
γ-Reθ obligatorio si Re<5e5, gas con FAR. **Sin y⁺_max, recuento de malla y
flag de convergencia, el par no entra en el flywheel.**

**Benchmarks públicos priorizados**: Aachen 1.5 (smoke test, tutorial SU2
público) → **SPLEEN C1** (VKI+Safran, Zenodo abierto: cascada LPT M₂ 0.7–0.95,
Re 70k–120k, con estelas y purga — el mejor caso moderno) → VKI LS89
(aerotérmica) → PAK-B/T106 (lapso de Re) → **NASA E³** (geometría+malla+BCs
públicas en data.gov) → MT1 (η de etapa HP con distorsión) → LISA (shroud y
purga por etapa).

### D.1 · GPU y CFD diferenciable: XLB, JAX-Fluids, NVIDIA Warp (síntesis del informe 06)

Auditoría hecha **clonando y leyendo el código fuente de los tres repos**.
Conclusión: **ninguna puede correr hoy un pasaje de turbina axial transónica
con y+≈1** — les faltan piezas físicas de primer orden, no de conveniencia:

- **XLB (Autodesk, Apache-2.0) — DESCARTAR**: es LBM **incompresible
  isotermo** (la clase se llama `IncompressibleNavierStokesStepper`);
  "Supersonic Flows" está en la *wishlist* del README, no en desarrollo; sin
  modelo de pared (la condición de existencia del LBM cartesiano a Re 10⁶ —
  sin él, 10¹²–10¹⁴ celdas), sin marco rotatorio ni periodicidad de paso; y
  su diferenciabilidad **solo funciona en el backend JAX** — el backend Warp
  (el rápido) devuelve **gradientes cero**, según el propio comentario del
  repo. El LBM transónico industrial existe (ProLB corrió el LS89;
  PowerFLOW llega a Mach 2) pero es otro método (HRR térmico + D3Q39 +
  overset) que XLB no tiene.
- **JAX-Fluids (TUM, MIT) — VIGILAR**: el mejor solver compresible
  diferenciable abierto (WENO/TENO, level-set, positivity, checkpointing
  para el backprop, 512 A100 probadas), pero **cartesiano por arquitectura**:
  sin malla body-fitted, sin marco rotatorio, sin periodicidad de paso, sin
  NSCBC, con γ y cp **constantes** (verificado en `ideal_gas.py`), sin RANS
  ni wall model. Rol acotado: banco de física canónica 2D (SBLI, burbuja de
  separación, mezcla de estela) y laboratorio de aprendizaje de cierres —
  **nunca L3, nunca pares de calibración**.
- **NVIDIA Warp (Apache-2.0) — DESCARTAR como base de solver**: no es un
  solver, es un compilador Python→CUDA con AD (con limitaciones serias:
  `*=` rompe el gradiente, los bucles dinámicos no se reproducen en el
  backward). De ~130 publicaciones del ecosistema, ~90 % son
  robótica/gráficos y **cero turbomáquinas**. Escribir el RANS de turbina en
  Warp = reconstruir Turbostream (3–5 años-persona). Nicho residual:
  geometría GPU (BVH/SDF) para 5a/5c si algún día hiciera falta.

**Dónde sí capturar el valor diferenciable** (el hallazgo accionable):

1. **El gradiente que Phy-AT necesita no es el del CFD.** Para calibrar los
   modelos de pérdidas por descenso de gradiente basta ∂y_L1/∂θ_loss — la
   derivada del **modelo barato** respecto a sus coeficientes; y_hifi es una
   constante. → **Escribir el meanline L0 (y opcionalmente el SCM) en
   `jax.numpy` desde el primer día** (coste ~0 en diseño inicial, carísimo
   como retrofit; corre en laptop CPU): sustituye la calibración afín de 2
   parámetros por una calibración física regularizada de 10–30 coeficientes
   que extrapola. Encaja con el marco KOH del informe 04.
2. **Empaquetar MULTALL como L3 interno opcional** (Meangen/Stagen generados
   desde el contrato → parser de la salida a (Γ, η_tt, ER, α_exit, perfiles
   radiales)): decenas de pares desatendidos en CPU, red de seguridad sin
   licencia comercial. *El mayor retorno por unidad de esfuerzo del
   informe 06.*
3. Si algún día hace falta **gradiente sobre la alta fidelidad** (forma):
   **SU2** — ya tiene mixing plane conservativo + NRBC de Giles + adjunto
   discreto validado en turbinas axiales. Días, no años.
4. **Vigilar sin invertir**: JAX-FVM (arXiv:2607.07385) y DiFVM
   (arXiv:2603.15920) — los primeros FV diferenciables sobre **malla no
   estructurada**, el verdadero horizonte que rompería la barrera
   "diferenciable ⇒ cartesiano". FluidX3D queda descartado por licencia
   (solo no comercial).

---

## E. Capas 2–4: ML, UQ y optimización (síntesis del informe 04)

**Veredicto global: la arquitectura de Phy-AC se conserva** — el campo ha
convergido exactamente a su patrón: *física de bajo orden interpretable +
corrección ML de la discrepancia + features físicas correctas* (Senior &
Miller, Whittle Lab, J. Turbomach. 146(4):041007, 2024). Formalmente, Phy-AC
**ya es** un modelo multifidelidad de Kennedy-O'Hagan (ensemble residual =
δ(x); calibración afín = ρ), y la GNN que aprende URANS−RANS (J. Turbomach.
148(1):011003) valida el residual-learning de forma independiente.

**Qué NO adoptar (con evidencia)**:
- **PINNs en capas 2–4**: coste ≥ CFD por caso, sesgo espectral, fallo en
  choques (= el régimen de una HPT transónica), reentrenamiento por condición
  de contorno; y el residuo verdad−L0 **no obedece ninguna PDE** — física en
  la loss sería decorativa. *Sí* como herramienta offline de asimilación de
  medidas para L2 (Hanrahan et al., J. Turbomach. 147(11), 2025).
- **FNO/DeepONet/GNNs**: necesitan 10³–10⁴ simulaciones (presupuesto: 150–500
  totales); PLAID documenta que **FNO se degrada en mallas anisótropas de
  turbomáquina**; MMGP (NeurIPS 2023) muestra que en datos escasos ganan los
  métodos gaussianos/reducidos — lo que refuerza el modelo pequeño sobre
  embedding físico.
- **MFBO, generativos (difusión/GAN/VAE), DRL, GP en lugar del ensemble**:
  razones detalladas en el informe (régimen de datos, dependencias, o falta
  de ventaja clara).

**Mejoras a adoptar (por valor/riesgo)**:
- **M1 · Conformal prediction split sobre el ensemble** (~50 líneas NumPy):
  la puerta de calidad deja de ser heurística — cobertura marginal
  garantizada, con recalibración de q_α por ronda y reporte de cobertura
  empírica (el aprendizaje activo viola la intercambiabilidad — Fannjiang
  et al., PNAS 2022 — y se documenta).
- **M2 · Ensemble K = 8–10** (en regresión M=5 no garantiza calidad de UQ;
  coste ~cero).
- **M5 · Embedding de turbina** (el trabajo de dominio): Zweifel por fila,
  M₂ de estátor y M₃_rel de rotor, **margen de garganta continuo** (sustituye
  al Koch SM), reacción mínima de cubo, deflexión máxima, Re mínimo, τ/h,
  AN², t_TE/o, swirl de salida, coordenadas de Smith, log(TR) y η_L0.
- **M4 · Calibración afín → afín con 2 features físicas** (ridge + LOO, con
  puerta: solo si LOO mejora).
- **M3 · Recalibración isotónica** (~30 líneas PAVA) para el diagrama de
  fiabilidad del informe de run.
- **M6–M8 · Adquisición**: ΔHipervolumen esperado por Monte Carlo (A/B contra
  el LCB actual), sesgo hacia la frontera g≈0 (con ~20 % de región factible),
  TuRBO-lite opcional.
- **M9 · SPLEEN C1 como ancla experimental de L2** — algo que Phy-AC no
  tiene: validación contra medida abierta y citable.
- **M10 · Hooks de corrección aprendida POR BLOQUE de pérdida**
  (perfil/secundaria/TE/fuga), no un factor global (Senior & Miller 2024).

**Riesgo señalado dos veces (informes 01 y 04)**: copiar los rangos de
`DESIGN_VARS` de Phy-AC produciría un espacio físicamente absurdo (ψ de
turbina llega a 2.5 vs 0.45 del compresor). Es lo primero a reescribir.

---

## F. Geometría y estructuras (síntesis del informe 05)

**Perfil: método de Pritchard (ASME 85-GT-219, 1985) — 11 parámetros.**
Es el único método donde **la garganta y el giro no guiado son ENTRADAS**
(cualquier otro obliga a generar-medir-iterar = optimización 2D encubierta),
es analítico y determinista (µs, sin fallos de convergencia), y sus 6
defaults son exactamente la lista de correlaciones de diseño a escribir (el
análogo de Lieblein+Carter). Su debilidad (G1, no G2 en LE y garganta) se
mitiga con **renderizador dual sobre los mismos 11 parámetros**:
`PRITCHARD11` (vóxel/STEP) y `BSPLINE_G2` opcional (CFD). Referencias de
implementación abiertas: pyturbo-aero (NASA, con TE de arco y cuñas, flag
aft_loaded, cálculo de garganta contra el álabe vecino), ParaBlade (G2),
T-Blade3 (comba por curvatura).

**3D**: apilado por centroide (reutilizando `polygon_section_props`) + ley
paramétrica de compound lean (midspan recto — la receta industrial) + sweep
opcional, emitidos como línea de apilado ya evaluada (la 5c sigue "tonta");
17 (HPT) / 21 (LPT) estaciones con densificación coseno; ≥100 puntos por
sección con clustering en LE/garganta/TE; endwall axisimétrico en v1 con el
bloque no-axisimétrico reservado (`enabled: false`).

**Arquitectura mecánica** (las 7 diferencias): NGV primero; **fir-tree
multi-lóbulo** (redundancia, poste del disco como sección crítica, agujeros
de alimentación); **anillos de tobera SEGMENTADOS** con bandas y feather
seals (un anillo continuo a 1500 K rompe por dilatación); tip shroud con
laberinto y Z-notch; cover plates; sellos inter-etapa; cavidades de purga.
**La capa 5c de Phy-AC es reutilizable en ~60 %** (loft sólido — que funciona
MEJOR con perfiles gruesos de TE romo —, discos, carcasa, paridad STL↔STEP);
lo nuevo: TurbineFirTree, NozzleSegment, TipShroud, InterstageSeal,
CoverPlate, RimCavity y un **ManufacturabilityCheck que rechaza (no clampea
en silencio) rasgos de sello < 2 vóxeles**.

**Honestidad del demostrador impreso** (bloque `manufacturing` del contrato):
no es una pieza caliente (sin monocristal, sin refrigeración interna, sin
film holes, sin TBC); demuestra forma, ajuste, montaje e interferencias — no
vida ni térmica.

**Estructural (9 g duros)**: fluencia y **burst al 120 %** (14 CFR 33.27 —
Phy-AC usa 1.05; aquí la cifra es regulatoria y citada), AN²_max **derivado
de la curva de creep** (no tabulado), raíz con K_t + **carga del shroud**
(20–40 % de la tensión de raíz), **creep de álabe y de rim por Larson-Miller**
(P_LM = T(20+log₁₀t_r), a velocidad de DISEÑO, con `t_life_h` como variable),
T_metal admisible vía cadena térmica (OTDF/RTDF → T₀_rel del rotor → η_c →
TBC), Campbell contra nº de NGV **y nº de inyectores** (con λ₁ de
empotrado-apoyado si hay shroud, y E(T) derrateado), y apoyo del abeto.
El solver de disco de Phy-AC se extiende con el término termoelástico
αE·dT/dr (misma tridiagonal; regresión: ΔT=0 → Timoshenko bit a bit).
Materiales calientes por PIEZA: IN-718 (discos), MAR-M247 (álabes/NGV),
CMSX-4 (referencia SX), con curvas de Larson-Miller a recopilar y documentar
en `validation/`.

---

## G. Recursos de software existentes (inventario consolidado)

| Recurso | Licencia | Qué aporta a Phy-AT |
|---|---|---|
| **TurboFlow** ([turbo-sim/turboflow](https://github.com/turbo-sim/turboflow), JOSS 2025) | MIT | **El recurso nº 1**: meanline de turbina equation-oriented con KO/Benner/Moustapha implementados y verificados aquí línea a línea; tres criterios de choking; casos Kofskey con geometría YAML + datos experimentales + script de error. Referencia de implementación y listón de precisión |
| **AxialOpt** ([turbo-sim/AxialOpt](https://github.com/turbo-sim/AxialOpt)) | MIT | AM/DC completos; conjunto de restricciones g con justificación bibliográfica anotada; correlaciones de stagger/espesor de Kacker |
| **NASA turbo-design (TD3)** (skill interna + [nasa/turbo-design](https://github.com/nasa/turbo-design)) | NASA | `TurbineSpool` con TD2/AM/KO/Craig-Cox/Traupel; export AGF para BladeGen/TurboGrid. **Cautela**: Phy-AC ya lo descartó como L1 (fase 11: ODE de equilibrio radial colgaba con >1 streamline) — usarlo como referencia de modelos y contraste, no como motor |
| **pyCycle** ([OpenMDAO/pyCycle](https://github.com/OpenMDAO/pyCycle)) | Apache | Algoritmo de refrigeración de Gauntner con constantes explícitas (verificado) |
| **MULTALL/Meangen/Stagen** (Denton 2017) | dominio público | La arquitectura de referencia meanline→geometría→solver; throughflow y Q3D de contraste |
| **pyturbo-aero** ([nasa/pyturbo-aero](https://github.com/nasa/pyturbo-aero)) | NASA-1.3 | Perfiles de turbina (TE de arco, aft_loaded, garganta) y apilado 3D con lean/sweep |
| **ParaBlade** ([NAnand-TUD/parablade](https://github.com/NAnand-TUD/parablade)) | permisiva | Parametrización G2 y canal meridional — el renderizador `BSPLINE_G2` |
| **T-Blade3** ([GTSL-UC/T-Blade3](https://github.com/GTSL-UC/T-Blade3)) | ver repo | comba por curvatura, integración ESP |
| **Pritchard ref.** ([DavidPoves/11-Parameters…](https://github.com/DavidPoves/11-Parameters-Turbine-Blade-Generator)) | — | implementación legible del método de perfil elegido |
| **turbigen** ([turbigen.org](https://turbigen.org/)) | abierto | la contra-tesis a responder (CFD GPU se salta el throughflow) y el dato de coste (~10 min/RANS en A100) |
| **T-AXI** (U. Cincinnati) | ejecutable | la alternativa Newton-acoplada (MTFLOW) si el lazo Novak resulta frágil |
| **SU2** ([su2code.github.io](https://su2code.github.io/)) | LGPL-2.1 | la vía realista al **adjunto sobre alta fidelidad**: mixing plane conservativo + NRBC de Giles + adjunto discreto (CoDiPack) validado en turbinas axiales |
| **JAX-Fluids** ([tumaer/JAXFLUIDS](https://github.com/tumaer/JAXFLUIDS)) | MIT | VIGILAR — banco de física canónica 2D y laboratorio de cierres diferenciables; cartesiano por arquitectura, nunca L3 (informe 06) |
| **XLB** ([Autodesk/XLB](https://github.com/Autodesk/XLB)) | Apache-2.0 | DESCARTADO — LBM incompresible isotermo, sin wall model ni marco rotatorio; gradientes cero en su backend rápido (informe 06) |
| **NVIDIA Warp** ([NVIDIA/warp](https://github.com/NVIDIA/warp)) | Apache-2.0 | DESCARTADO como base de solver (no es un solver); nicho residual de geometría GPU para 5a/5c (informe 06) |
| **PyFR** ([pyfr.org](https://pyfr.org/)) | BSD-3 | LES/ILES de cascada en GPU (T106C, MTU-T161) — generador de priors L4, fuera del lazo |
| **SMT (MFCK)** | BSD | co-kriging recursivo cuando el flywheel tenga ≥15 pares |
| **turbodesigner** (OpenOrion) | MIT | el ancestro de la 5a de Phy-AC — **sin turbinas** (verificado): la turbina hay que escribirla |
| **Datos**: SPLEEN C1 (Zenodo), NASA E³ (data.gov: geometría+malla+BCs), Kofskey (NTRS + YAMLs de TurboFlow), LS89, PAK-B, T106, MT1, LISA | abiertos | pares semilla de calibración y campaña de validación |

---

## H. El contrato `phyat-axial-1` y la frontera con Phy-CB (pendiente)

El contrato hereda la disciplina de `phyac-axial-2` (schema JSON versionado
con MAYOR en el identificador, validador sin dependencias, consumidor que
rechaza lo que no entiende, `nSCHEMA_MAJOR` sincronizado y testeado) y añade
los bloques de turbina — el esquema completo está en el informe 05 §6.3. Las
decisiones estructurales:

1. **`profile_params` (los 11 de Pritchard) ADEMÁS de `points`** — cualquier
   consumidor re-renderiza a otra resolución o en G2 sin correr la capa 5a.
2. **`index_markers`** por sección (LE, TE, garganta, dorso/vientre) — lo que
   un mallador CFD necesita y el zipper de la 5c no.
3. **Garganta y o/s POR SECCIÓN** (el twist las hace variar con el radio).
4. **`thermal` como bloque de primer nivel** (lo consumen estructural, CFD y
   refrigeración): T₀₃/T₀₄, OTDF/RTDF/n_injectors, perfil radial, y por fila
   T₀_abs/T₀_rel/η_c/T_coolant/T_metal/TBC.
5. **`interfaces` con `requires`/`provides` explícitos** — el mecanismo de
   reutilización entre sistemas.
6. **`manufacturing`** — el bloque de honestidad del demostrador
   (not_represented, suppressed_features, warning).
7. **`provenance`** — método de perfil, regla de ángulo de salida, modelo de
   pérdidas, modelo de creep, citas.

### La frontera Phy-CB → Phy-AT (PENDIENTE — sin acceso a Phy-CB hoy)

El contrato de Phy-AC se reutilizará en Phy-CB (cámaras de combustión). Para
que el de Phy-AT case en la frontera cuando Phy-CB exista, el bloque
`interfaces.upstream` declara ya **qué debe entregar un combustor**:

```
upstream: {
  system: "phycb-annular-1" | null,     # null hasta que Phy-CB exista
  plane_z_mm, r_hub_mm, r_tip_mm,       # el plano de acople
  requires: [
    T0_radial_profile,                  # RTDF — perfil radial de T0
    OTDF,                               # pico circunferencial (dimensiona el NGV
                                        #   y la vida del álabe al hot streak)
    n_injectors,                        # = orden de excitación de Campbell del NGV
    swirl_alpha_deg_profile,            # swirl residual (lean burn: no despreciable)
    P0_profile,                         # perfil de P0 (pérdida de carga del combustor)
    coolant_available_frac              # aire disponible para refrigerar la turbina
  ]
}
coolant_source: { from: "phyac-axial-2", bleed_stage, P0_Pa, T0_K, frac }
```

Con esto, tres consecuencias quedan resueltas por diseño: (a) mientras no
haya Phy-CB, Phy-AT usa defaults declarados (OTDF/RTDF típicos de la
literatura — MT1/FACTOR como referencia) y lo marca en `provenance`; (b) el
día que Phy-CB exista, el casamiento es la validación de un JSON contra un
schema, no una integración ad hoc; (c) el mismo patrón sirve aguas abajo
(`downstream.provides`: perfil de T₀, swirl, Mach — para una LPT o el
escape). El sangrado del compresor (Phy-AC `derived.bleed`) es el precedente
directo: aquí se generaliza a una **interfaz nombrada y versionada**.

---

## I. Plan de validación y máquinas de referencia

**Protocolo idéntico a Phy-AC**: geometría y punto de operación publicados;
el modelo se califica en su predicción de pérdidas; **cero recalibración por
máquina**; anclas de regresión congeladas antes de la campaña.

**Fase 1 — aire frío NASA** (sin ambigüedad de definición de η; geometrías ya
digitalizadas en los YAML de TurboFlow, licencia MIT):

| Máquina | Modo de choking | Tolerancias (base: métricas publicadas de TurboFlow) |
|---|---|---|
| Kofskey 1974 (TN D-7625), 1 etapa | **estátor** | ṁ ±2.5 % · par ±5 % · α_exit ±5° |
| Kofskey 1972 (TN D-6967), 1ª etapa | **rotor** | ṁ ±2.5 % · par ±5 % · α_exit ±5° |
| Kofskey 1972, 2 etapas | — | ṁ ±2.5 % · par ±5 % · α_exit ±2.5° |
| Kofskey 1974 (E-7776), estátor abierto | barrido de garganta | ṁ ±3 % (test de la regla de gauging) |

Mapas a 70/90/100/110 % de ω, PR 1.6–4.5. **El caudal debe salir a <2.5 % en
el 100 % de los puntos o hay un bug** (la garganta lo gobierna todo — el
contraste con el +6.6 % de choke que Phy-AC arrastra en compresor).

**Fase 2 — E³ HPT** (CR-168289 + hardware NTRS 19850002687): refrigerada;
**declarar la definición de η** (la diferencia termodinámica/ciclo-GE vale
2.5 pts); tolerancia ±2 pts. El análogo directo del E³ HPC de Phy-AC.

**Fase 3 — ampliación**: series NASA de alta temperatura (NTRS 19680006274,
19720024133), turbina de potencia libre (19790009688), E³ LPT (localizar el
CR), LISA (Y_cl con/sin shroud, purga), y las anclas de cascada (SPLEEN C1,
LS89, PAK-B/T106) para la calibración L2.

**Verificación** (suite estilo `test_phyac.py`): identidades de triángulos
(1e-6), conversión de convenciones de ángulo, Y_definición vs Y_modelo por
fila, vórtice general vs forma cerrada (T21 reutilizado), actuator disc
(curvatura), tobera bloqueada analítica, rotalpía/h₀ conservadas, balance de
energía con refrigeración cerrado (assert duro), continuidad de g a través
del choke y de la activación de refrigeración, disco térmico vs Timoshenko
(ΔT=0 bit a bit), reproducibilidad de semilla.

---

## J. Hoja de ruta propuesta y riesgos

**Orden de construcción** (dependencias mínimas, cada fase con su
verificación antes de la siguiente):

1. **`physics_turbine.py` (L0)**: gas cp(T,FAR); triángulos (φ,ψ,R);
   AMDC-KO+Benner; Aungier; Mach crítico analítico; marcha de entropía;
   Gauntner+Hartsel; g A–D. **Escrito en `jax.numpy` desde el primer día**
   (informe 06: coste ~0 en diseño inicial, habilita la calibración por
   gradiente de los coeficientes de pérdidas; sigue corriendo en laptop
   CPU). → **Fase 1 de validación (Kofskey).**
2. **`structures_turbine.py` (L1s)**: materiales calientes por pieza +
   Larson-Miller; disco con término térmico; 9 g. → regresión Timoshenko +
   anclas de creep documentadas.
3. **Capas 2–4**: adaptador de dominio (embedding M5, rangos de θ);
   conformal M1 + K=8–10 (M2); el core de `neural_optimizer.py` se porta
   casi verbatim (como Phy-CC→Phy-AC).
4. **`turbine_profiles.py` + `geometry_generator` (5a)**: Pritchard dual;
   contrato `phyat-axial-1` + validador; exportador TurboGrid.
5. **`scm_turbine.py` (L1)**: esqueleto → curvatura/actuator disc → choking →
   pérdidas Benner spanwise → refrigeración en tres lazos. → bench estilo
   `bench_scm.py`.
6. **`AxialTurbineDesigner` (5c)**: reutilización de la base + FirTree de
   turbina + NozzleSegment (prioridad máxima) → shroud/sellos/cover plates →
   paridad STL↔STEP + interferencias.
7. **L2/L3**: emisor de BCs de turbina + parser de retorno; pares semilla
   (Aachen → SPLEEN → LS89 → E³); calibración estratificada; **MULTALL
   empaquetado como L3 interno opcional** (Meangen/Stagen desde el contrato
   → parser de salida → pares desatendidos en CPU — informe 06, el mayor
   retorno por esfuerzo).

**Riesgos principales** (consolidados de los cinco informes):

| Riesgo | Mitigación |
|---|---|
| Copiar rangos/correlaciones de Phy-AC (ψ, f_Re, choke como fallo, T_METAL_FRAC=1) | Los cinco informes lo marcan: **reescribir, no trasplantar** — checklist explícita en §A |
| Convención de signos de ángulos (la fuente nº 1 de bugs en código de turbinas) | función única ax↔tan + tests de identidad desde el día 1 |
| Doble solución sub/supersónica en el SCM | bisección acotada a M_crit + rama supersónica; requisito de v1 (lección OTAC) |
| Definición de η refrigerada (2.5 pts de ambigüedad) | validar primero en aire frío; campo `efficiency_definition` obligatorio en los pares |
| Deriva L0↔geometría↔L1 (el PR_WINDOW de Phy-AC, peor aquí porque la garganta fija el gasto) | `MassFlowInfeasible` etiquetado con el gasto alcanzable; ventana sobre el **gasto reducido** |
| Números de segunda mano ([S]/[B]/[M]/[?] en los informes) | resolver los DOI con red abierta antes de escribir `Quasar_PhyAT_Science.md`; ningún [?] como ancla |
| Rasgos de sello vs vóxel (laberintos de 0.2–0.4 mm) | ManufacturabilityCheck que rechaza o marca, nunca clampea en silencio |
| Frontera Phy-CB inexistente hoy | defaults declarados de OTDF/RTDF/swirl + `interfaces.requires` versionado (§H) — **pendiente explícito** |

---

## K. Referencias clave (selección; listas completas en cada informe)

**Meanline y pérdidas**: Ainley & Mathieson (1951) ARC R&M 2974; Dunham &
Came (1970); **Kacker & Okapuu (1982)** J. Eng. Power 104(1); Moustapha,
Kacker & Tremblay (1990); **Benner, Sjolander & Moustapha (1997, 2006 I/II)**;
Craig & Cox (1970); Traupel (2001); Denton (1993) IGTI Scholar Lecture;
Coull & Hodson (2012/2013); Zhu & Sjolander (2005); **Aungier (2006),
*Turbine Aerodynamics***; Smith (1965); Zweifel (1945); Glassman SP-290;
Dixon & Hall; Moustapha et al., *Axial and Radial Turbines* (2003).

**Through-flow**: Wu (1952) NACA TN 2604; Smith (1966); Novak (1967);
Wilkinson (1969); **Denton (1978)** J. Eng. Power 100(2); Casey & Robinson
(2010); **Tiwari, Stein & Lin (2013)**; Petrović & Wiedermann (2013);
Lewis (1994); Hendricks (2016) AIAA 2016-0119; Stodola / Cooke (1985);
Genrup (2005) / Sammak (2013) — LUAX-T.

**Refrigeración**: **Gauntner (1980) NASA TM-81453**; **Hartsel (1972) AIAA
72-11**; **Young & Wilcock (2002)** J. Turbomach. 124(2) Parts 1–2; Young &
Horlock (2006).

**Alta fidelidad y calibración**: **Denton (2010) GT2010-22540**; Sandberg &
Michelassi (2022) ARFM 54; Hodson & Howell (2005) ARFM 37; Lee, Dawes &
Coull (2021) JGPPS 5; Brandvik & Pullan (2011); **Kennedy & O'Hagan (2000,
2001)**; Le Gratiet & Garnier; Peherstorfer, Willcox & Gunzburger (2018)
SIAM Rev 60; Garzon & Darmofal (2003); Montomoli et al. (Springer).

**ML/UQ/optimización**: **Senior & Miller (2024)** J. Turbomach. 146(4);
Lakshminarayanan et al. (2017); **Gopakumar et al. (2025)** conformal para
surrogates, MLST; Kuleshov, Fenner & Ermon (2018); Fannjiang et al. (2022)
PNAS; Barber et al. (2023) Ann. Statist.; Deb et al. (2002); Daulton et al.
qNEHVI/MORBO; Eriksson et al. TuRBO/SCBO; Casenave et al. MMGP (NeurIPS 2023);
Zou et al. (2024) AI Review 57.

**Geometría y estructuras**: **Pritchard (1985) ASME 85-GT-219**; Agromayor
et al. (2021) CAD 133 (ParaBlade); Kulfan (2007); Bagshaw et al. (2008);
Harvey et al. (2000) endwall; Larson-Miller; **14 CFR 33.27** (burst 120 %);
NETL Gas Turbine Handbook (§4.2.2.1 refrigeración, §4.4.2 TBC); Rosic &
Denton (2008) shroud; Yaras & Sjolander (1992).

**Software y datos**: TurboFlow (JOSS 2025, MIT); AxialOpt (MIT); NASA TD3;
pyCycle; MULTALL (Denton 2017); pyturbo-aero (NASA); ParaBlade; T-Blade3;
turbigen; SMT; **SPLEEN C1 (Zenodo 10.5281/zenodo.7264761)**; **NASA E³
(data.gov, geometría+malla+BCs)**; Kofskey TN D-6967 / TN D-7625 / E-7776;
VKI LS89; PAK-B; T106; MT1; LISA; PicoGK (LEAP 71).
