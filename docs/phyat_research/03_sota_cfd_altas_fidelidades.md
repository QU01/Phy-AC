# Estado del arte en diseño y simulación de alta fidelidad de turbinas axiales (2015–2026)
### Informe de fundamentación para el sistema **Phy-AT (Quasar)**, análogo de Phy-AC

> **Nota metodológica**: el proxy de egress de la sesión de investigación bloquea `WebFetch` a los dominios académicos (arXiv, ASME, MDPI, Zenodo, ScienceDirect, NTRS…), por lo que la evidencia procede de búsquedas web con extracción de resúmenes/abstracts, no de lectura de texto completo. Los valores numéricos marcados con ⚠️ deben re-verificarse contra el PDF original antes de publicarlos.

---

## 0. Resumen ejecutivo

1. **La cadena industrial no ha cambiado de topología desde 2010; ha cambiado el coste por escalón.** Sigue siendo ciclo → meanline (1D) → throughflow (2D axisimétrico) → blade-to-blade Q3D → RANS 3D multi-fila con mixing-plane → URANS/armónicos selectivo → LES puntual. Lo nuevo es que un RANS 3D multi-fila con fugas ya cabe en **~10 minutos sobre una sola GPU**, lo que mueve el RANS desde "verificación final" hacia "bucle de exploración" ([Brandvik & Pullan, 2011](https://asmedigitalcollection.asme.org/turbomachinery/article-abstract/133/2/021025/468994/An-Accelerated-3D-Navier-Stokes-Solver-for-Flows); [turbigen docs](https://turbigen.org/solver.html)).
2. **La advertencia de Denton (2010) sigue vigente y es la base honesta de cualquier lazo de calibración**: el CFD de turbomáquinas debe usarse *comparativamente*, no como predictor cuantitativo absoluto de prestaciones, porque los errores dominantes no son numéricos sino de **condiciones de contorno desconocidas** (perfiles de P0/T0 de entrada), **geometría desconocida** (juegos de punta, radio de borde de ataque) y **suposición de flujo estacionario** ([Denton, GT2010-22540](http://proceedings.asmedigitalcollection.asme.org/proceeding.aspx?articleid=1609811)). Esto es exactamente el argumento a favor del enfoque "data flywheel" de Phy-AC/Phy-AT: no intentes predecir η en absoluto, calibra el modelo barato con pares medidos en tu propio dominio.
3. **La diferencia esencial turbina vs. compresor para el diseño de Phy-AT**: en turbina la **transición laminar-turbulenta y el número de Reynolds** son de primer orden (LPT a crucero baja a Re≈1e5), la **refrigeración y las fugas de purga** cambian el rendimiento en cantidades comparables a todo lo que optimizas aerodinámicamente (0.7–1.2 % de η por cada 1 % de purga), y la **capacidad (gasto corregido) es un output tan importante como η** porque fija el matching del motor. Un Phy-AT que solo calibre η estará calibrando la mitad del problema.
4. **Recomendación central**: escalera L0 (meanline, ~1 ms) → L1 (throughflow SCM, ~1–5 s) → **L2 (capa Q3D/transición + corrección afín-más-residual calibrada con pares, ~0.1–10 s)** → **L3 (RANS 3D multi-fila externo con emisión de BCs y reinyección de pares, ~10 min–2 h)**, con LES/DNS **fuera del lazo**, usados solo para construir priors y para casos-testigo. Declarar límites: ±1.5–2.5 puntos de η en L0, ±0.3–0.5 puntos dentro del envolvente de calibración en L2/L3, y **nunca por debajo del suelo de incertidumbre del banco (0.25–0.45 % U95)**.

---

## 1. Flujo de diseño industrial moderno de turbinas axiales

### 1.1 La cadena canónica

La secuencia que hoy usan GE, Safran, Rolls-Royce, Siemens Energy y MHI (y que los OEM pequeños replican con herramientas abiertas):

| Escalón | Herramienta típica | Qué decide | Coste |
|---|---|---|---|
| Ciclo termodinámico | GasTurb, NPSS, PROOSIS, in-house | ΔH0 por eje, ṁ, T4, sangrados y refrigeración | ms |
| Meanline 1D | in-house, Meangen (Denton), AxSTREAM, Concepts NREC | nº de escalones, φ/ψ/Λ, radios, α/β, Zweifel, nº de álabes | ms |
| Throughflow 2D (SLC o time-marching) | in-house SCM, MULTALL-throughflow, TRACE-TF | distribución radial de trabajo, torsión, bloqueo de pared, mezcla radial | s |
| Blade-to-blade Q3D | **MISES**, códigos Euler+IBL propietarios | forma de perfil, carga (fore/aft-loading), transición, difusión trasera | s |
| RANS 3D multi-fila (mixing plane) | CFX, Fluent, TRACE, elsA, HYDRA, Turbostream, MULTALL | apilado 3D, lean/sweep, contorneado de endwall, juego de punta, cavidades | min–h |
| URANS / balance armónico | TRACE, elsA, HYDRA, CFX | interacción estela-álabe, clocking, forzado aeroelástico | h–días |
| LES / DNS selectivo | HiPSTAR, AVBP, Nektar++, PyFR, YALES2 | física de pérdidas, transferencia de calor, rugosidad, hot streaks | días–semanas en HPC |

Jerarquía descrita en [Denton & Dawes, 1998, Proc. IMechE Part C 212(2)](https://journals.sagepub.com/doi/10.1243/0954406991522211) y en la revisión más citada de la década: **Sandberg & Michelassi (2022), "Fluid Dynamics of Axial Turbomachinery: Blade- and Stage-Level Simulations and Models", *Annual Review of Fluid Mechanics* 54:255–285** ([enlace](https://www.annualreviews.org/content/journals/10.1146/annurev-fluid-031221-105530)).

### 1.2 Referencias de revisión imprescindibles

- **Denton, J.D. (2010). "Some Limitations of Turbomachinery CFD". ASME GT2010-22540** ([ASME](http://proceedings.asmedigitalcollection.asme.org/proceeding.aspx?articleid=1609811)). Taxonomía de errores: (i) BCs desconocidas, (ii) geometría desconocida (juegos, LE), (iii) hipótesis de estacionariedad, (iv) error numérico, (v) error de modelado (turbulencia). Conclusión operativa: **usar el CFD de forma comparativa**.
- **Sandberg & Michelassi (2022)**, *Annu. Rev. Fluid Mech.* 54:255–285.
- **Tyacke, Vadlamani, Trojak, Watson, Ma & Tucker (2019). "Turbomachinery simulation challenges and the future". *Progress in Aerospace Sciences* 110:100554** ([ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0376042119300715)).
- **Tucker & Tyacke (2016). "Future Directions of High-Fidelity CFD for Aero-Thermal Turbomachinery", AIAA 2016-3322** ([PDF](https://cfd.ku.edu/papers/AIAA-2016-3322.pdf)).
- **Hodson & Howell (2005). "Bladerow Interactions, Transition, and High-Lift Aerofoils in Low-Pressure Turbines". *Annu. Rev. Fluid Mech.* 37:71–98** ([Annual Reviews](https://www.annualreviews.org/content/journals/10.1146/annurev.fluid.37.061903.175511)). Referencia obligada de LPT: las estelas de la fila aguas arriba son la fuente dominante de inestacionariedad y disparan la transición de la capa límite del lado de succión.
- **Lee, W.Y., Dawes, W.N. & Coull, J.D. (2021). "The required aerodynamic simulation fidelity to usefully support a gas turbine digital twin for manufacturing". *J. GPPS* 5:15–27, DOI 10.33737/jgpps/132007** ([GPPS](https://journal.gpps.global/The-required-aerodynamic-simulation-fidelity-to-usefully-support-a-gas-turbine-digital,132007,0,2.html)). **El paper más directamente relevante para Phy-AT**: construye explícitamente una escalera de fidelidades sobre un rotor de HPT civil real, desde secciones Q3D con solver Euler hasta RANS 3D multi-paso y multi-escalón, y pregunta *cuán realista debe ser la simulación para que la predicción sea útil*.

### 1.3 Herramientas abiertas que definen el "suelo" del estado del arte

- **MULTALL / Meangen / Stagen (Denton, 2017)**, *J. Turbomach.* 139(12):121001 (GT2017-63993) ([ASME](https://asmedigitalcollection.asme.org/turbomachinery/article-abstract/139/12/121001/378803); [sitio oficial](https://sites.google.com/view/multall-turbomachinery-design)). **Es literalmente la arquitectura que Phy-AT debería replicar**: `Meangen` (meanline 1D) → `Stagen` (generación y apilado de geometría) → `Multall` (Navier-Stokes 3D multi-escalón), y `Multall` también hace Q3D blade-to-blade y throughflow axisimétrico. Denton lo publica explícitamente para "individuos o pequeñas empresas sin sistema in-house" — que es el nicho de Phy-AT.
- **MISES 2.63 (Drela & Youngren)** ([manual MIT](https://web.mit.edu/drela/Public/web/mises/mises.pdf)). Euler estacionario acoplado a ecuaciones integrales de capa límite, con efectos Q3D. Para secciones 2D fuera de separación masiva, **precisión comparable a RANS estacionario con 2–3 órdenes de magnitud menos de coste** — exactamente lo que debe ser un escalón Q3D de turbina.
- **turbigen (Brind & Pullan, Cambridge)** ([docs](https://turbigen.org/)) — sistema abierto moderno de exploración de espacio de diseño con Turbostream 3 en GPU.

---

## 2. CFD RANS para turbinas: modelos, interfaces, precisión y mallado

### 2.1 Turbulencia y transición

**El modelo por defecto de la industria es k-ω SST (Menter, 1994)**, casi siempre acoplado a un modelo de transición. Tres familias vivas 2015–2026:

1. **γ-Reθ (Langtry & Menter, 2009, AIAA J.)** — dos ecuaciones adicionales sobre SST; cubre transición natural, bypass e inducida por separación ([SU2 V&V](https://su2code.github.io/vandv/LM_transition/)).
2. **Energía cinética laminar (LKE) / k-v'²-ω** — fenomenológicos, populares en el mundo LPT.
3. **Modelos basados en correlaciones propietarias** (TRACE de DLR, HYDRA de RR).

**Por qué la transición importa más en LPT que en cualquier otro componente:**
- A crucero, **Re basado en cuerda axial cae hasta ~1e5** ([revisión de efectos de Re en LPT](https://asmedigitalcollection.asme.org/GT/proceedings-abstract/GT2012/44748/1121/365369)).
- Las pérdidas escalan con **Re⁻¹ en régimen subcrítico y Re⁻¹ᐟ² en supercrítico**; la separación puede **aumentar el coeficiente de pérdida en ~300 % para Re < 1e5** ⚠️.
- La caída de rendimiento entre despegue y crucero puede ser **~2 % en motores grandes y hasta ~7 % en motores pequeños a gran altitud** ⚠️.

**Límite honesto a declarar**: incluso γ-Reθ, aunque predice una burbuja de separación corta en el T106A donde los modelos estándar no predicen separación alguna, **produce una burbuja cuyas características difieren mucho de experimento y DNS** ([discusión T106](https://arxiv.org/abs/2004.10967)).

Confirmación cuantificada sobre caso abierto moderno: **"The Impact of Transition and Turbulence Modelling on the SPLEEN High-Speed Low-Pressure Turbine Cascade", ASME GT2025-153288 / J. Turbomach. DOI 10.1115/1.4069487** ([ASME](https://asmedigitalcollection.asme.org/GT/proceedings/GT2025/88865/V010T30A019/1220713)): los cierres lineales son adecuados con flujo adherido incluso transónico, pero **fallan al predecir la mezcla de estela tras burbujas de separación laminar**. Ver también [RANS Prediction of Losses and Transition Onset in a High-Speed LPT Cascade, Energies 16:7348 (2023)](https://doi.org/10.3390/en16217348).

### 2.2 Tratamiento rotor-estator

| Método | Coste | Qué conserva | Cuándo usarlo en Phy-AT |
|---|---|---|---|
| **Mixing plane (stage interface)** | 1 paso por fila, estacionario | flujos promediados circunferencialmente; **no** transporta estelas ni choques | **Estándar para L3** |
| **Frozen rotor** | igual | posición relativa congelada; resultado depende de la posición | ❌ **Prohibirlo para prestaciones**; solo diagnóstico local |
| **Sliding mesh / URANS** | ×10²–10³ | inestacionariedad completa | fuera del lazo |
| **Balance armónico / NLH** | ×20 a ×10³ más barato que URANS | armónicos dominantes | opcional L3+ |

El mixing plane necesita **condiciones de contorno no reflectantes** y suficientes bandas de promediado ([ANSYS Fluent theory](https://ansyshelp.ansys.com/public//Views/Secured/corp/v251/en/flu_th/flu_th_sec_mpm.html); [Implicit and Conservative Mixing-Plane Method, AIAA J., DOI 10.2514/1.J062064](https://doi.org/10.2514/1.J062064)). El método armónico no lineal es **2–3 órdenes de magnitud más rápido que URANS** y el balance armónico 1–2 órdenes ⚠️ ([Harmonic Balance for Multistage Turbomachinery, GT2014](https://asmedigitalcollection.asme.org/GT/proceedings-abstract/GT2014/45615/V02BT39A005/250350)).

### 2.3 ¿Cuántos puntos de η se le pueden creer a un RANS de turbina?

- **En valor absoluto: poco.** Denton (2010) es explícito. En capacidad, Burdett, Hambidge & Povey (2021) declaran que **la precisión CFD en capacidad absoluta es menor que la experimental por un margen apreciable** ([Proc. IMechE Part A, 2021](https://journals.sagepub.com/doi/full/10.1177/0957650920909718)).
- **En diferencias (Δη entre variantes, misma malla/topología/BCs): bastante.** Es el uso "comparativo" que sustenta las ganancias medidas de 0.5–0.9 puntos por endwall contouring (§5.2).
- **RANS vs. URANS**: diferencia en pérdida secundaria ~10 % (de la pérdida secundaria, no de η) ⚠️.
- **Suelo experimental**: en la QinetiQ Turbine Test Facility, U95 del rendimiento = **0.45 % absoluto, 0.25 % relativo** ([Turbine Efficiency Measurement System, QinetiQ TTF](https://www.researchgate.net/publication/245355150_Turbine_Efficiency_Measurement_System_for_the_QinetiQ_Turbine_Test_Facility)). **Ninguna calibración puede declarar precisión por debajo de esto.**

**Cifra recomendada a declarar en Phy-AT**: RANS 3D bien planteado ⇒ **η absoluto ±1–2 puntos**, **Δη entre variantes ±0.3–0.5 puntos**, **capacidad ±1–2 %** sin calibrar. En LPT a bajo Re con separación, degradar un factor ~2.

### 2.4 Mallado

- **Densidad**: de 100–200 mil celdas por canal en los 90 a **0.5–1.0 millones por canal hoy** ⚠️ ([Mesh Resolution Effect on 3D RANS Turbomachinery](https://arxiv.org/pdf/1609.00063)).
- **y⁺ ≈ 1** para resolución de pared (obligatorio con γ-Reθ o flujo de calor); y⁺ ≈ 30 solo con funciones de pared, nunca para transición/térmica.
- **Topología O-H** — TurboGrid la genera automáticamente ([TurboGrid](https://www.ansys.com/products/fluids/ansys-turbogrid)); AutoGrid es el equivalente NUMECA.
- Un solo criterio (y⁺) no basta: forma de celda, celdas en capa límite, razón de expansión y suavidad determinan conjuntamente la precisión.

---

## 3. Escalado y alta fidelidad: LES, DNS, GPU

### 3.1 Códigos que definen el campo

| Código | Origen | Tipo | Notas |
|---|---|---|---|
| **Turbostream 3** | Brandvik & Pullan (Cambridge) | RANS/URANS estructurado, GPU | **3 escalones de turbina con caminos de fuga completos en <10 min con 4 GPUs** ([J. Turbomach. 133(2):021025, 2011](https://asmedigitalcollection.asme.org/turbomachinery/article-abstract/133/2/021025/468994)) |
| **HiPSTAR** | Sandberg (Melbourne) | DNS/LES compresible, GPU | HPT con rugosidad de microescala: **10–20 mil millones de puntos**, semanas en **Frontier**; campañas con 576–1728 GPUs V100 en Summit ([OLCF 2026](https://www.olcf.ornl.gov/2026/01/27/frontier-provides-high-fidelity-insights-into-turbine-aerothermal-performance/)) |
| **MULTALL** | Denton | RANS 3D multi-escalón + Q3D + throughflow | abierto, CPU, minutos ([GitHub mirror](https://github.com/paopaoai11/Multall-open-18.3)) |
| **TRACE** | DLR | RANS/URANS/armónico, transición | LPT 2 escalones rueda completa ([J. Turbomach. 141(5):051012](https://asmedigitalcollection.asme.org/turbomachinery/article-abstract/141/5/051012/368129)) |
| **elsA** | ONERA | multipropósito ([Mechanics & Industry 14:159–174](https://www.mechanics-industry.org/articles/meca/pdf/2013/03/mi130017.pdf)) | |
| **HYDRA** | Rolls-Royce | RANS + **adjunto discreto** | objetivos aero, termo y aeroacústicos |
| **AVBP / TurboAVBP** | CERFACS | LES compresible no estructurada | LES de la etapa MT1 completa |
| **YALES2** | CORIA | LES compresible, VF 4º orden | **WRLES de SPLEEN C1 a Re₂,is=70k** ([IJTPP 10(3):21, 2025](https://doi.org/10.3390/ijtpp10030021)) |
| **Nektar++ / PyFR** | Imperial | alto orden DG/FR, GPU | ILES de T106C |

(Nota: STREAmS no es un código de turbomáquinas — DNS de capa límite compresible.)

### 3.2 Coste: la ley que hay que respetar

- **RANS → LES resuelta en pared a Re ~1e6: ~10⁴** ([Slotnick et al., NASA/CR-2014-218178, CFD Vision 2030](https://komahanb.github.io/files/publications/nasa-cfd-2030.pdf)).
- **WRLES ~Re^2.4**; WMLES ~N_bl² ⚠️.
- Por fila: **RANS ~10⁶ celdas**, **híbrido ~10⁷**, **WRLES ~10⁸–10⁹**, **DNS ~10⁹–10¹⁰**.
- CFD Vision 2030: los métodos dominados por LES para problemas de ingeniería a Re altos **probablemente no serán viables ni en 2030**. Para Phy-AT: **LES nunca entra en el lazo de optimización; solo priors y casos-testigo.**

### 3.3 El punto de inflexión GPU

En turbigen con Turbostream 3, **una simulación RANS ≈ 10 minutos en una NVIDIA A100**, y el bucle de corrección de pérdida/desviación/incidencia ≈ 10 simulaciones ≈ **1 hora de reloj** ([turbigen: Flow solvers](https://turbigen.org/solver.html)). Coste de un "par hi-fi" moderno de L3: **~1 h por punto convergido**, no días. Con eso, **15–40 pares de calibración son un fin de semana de máquina**.

---

## 4. Efectos específicos de turbina que la alta fidelidad captura y las correlaciones no

### 4.1 Interacción cámara de combustión–turbina (hot streaks, swirl residual)

- **Efecto Kerrebrock–Mikolajczak (1970)**: el gas caliente tiene mayor velocidad absoluta → mayor incidencia → **impacta preferentemente en la cara de presión** del rotor.
- **Magnitud**: la distorsión de temperatura incrementa la carga térmica del álabe **hasta un 10–30 %, principalmente en la cara de presión** ([Analysis of Hot Streak Effects on Turbine Rotor Heat Load, J. Turbomach. 119(3):544, 1997](https://asmedigitalcollection.asme.org/turbomachinery/article-abstract/119/3/544/438359)).
- **Proyecto FACTOR (EU)**: simulador de cámara lean-burn 360° + turbina de 1.5 escalones en NG-Turb (DLR Göttingen); clocking cámara-NGV ([Aerospace 8(10):285, 2021](https://www.mdpi.com/2226-4310/8/10/285/htm); [LES of the Combustor Turbine Interface, GT2016](https://asmedigitalcollection.asme.org/GT/proceedings-abstract/GT2016/49798/V05BT17A003/239786)).
- **Swirl residual sobre η**: medido en MT1 ([Effect of Combustor Swirl on Transonic HP Turbine Efficiency, J. Turbomach. 136(1):011002, 2014](https://asmedigitalcollection.asme.org/turbomachinery/article-abstract/136/1/011002/378050)).
- **Consecuencia para Phy-AT**: la BC de entrada de un HPT es un campo 2D (r,θ) de P0, T0 y ángulos con clocking, no un escalar. Ver §8.3.

### 4.2 Purga de cavidades / rim seals

- **0.8 % de η_tt por cada 1 % de gasto de purga** en un escalón altamente cargado ⚠️; con distinto endwall contouring: **1.2 %/% y 0.7 %/%** ([GPPS: rim seal purge + endwall contouring](https://journal.gpps.global/Aerodynamic-influence-of-rim-seal-purge-flow-injection-on-the-main-flow-in-a-1-5,162078,0,2.html)).
- Dos mecanismos: mezcla de la purga con el flujo anular, e **intensificación de los flujos secundarios del rotor**.
- Datos: [Experimental Investigation of Purge Flow Effects on a HP Turbine Stage, J. Turbomach. 137(4):041006](https://asmedigitalcollection.asme.org/turbomachinery/article-abstract/137/4/041006/378534); [Fully Purged HP Turbine Stage, IJTPP 8(3):22](https://www.mdpi.com/2504-186X/8/3/22).

### 4.3 Fugas de punta (no carenado) y de shroud (carenado)

- **No carenado**: el juego de punta es responsable de **~1/3 de las pérdidas aerodinámicas del escalón** ⚠️ ([Huang, MIT](http://dspace.mit.edu/bitstream/handle/1721.1/67067/758647434-MIT.pdf?sequence=2), sobre Denton 1993). RANS captura mal la **rotura del vórtice de fuga**.
- **Carenado**: Rosic & Denton — la mezcla por diferencia de velocidad tangencial entre flujo principal y de fuga crea gran parte de la pérdida del shroud; álabes de enderezamiento en la cavidad de salida dieron **+0.4 % medido al freno** ([Control of Shroud Leakage Loss, J. Turbomach. 130(2):021010, 2008](https://asmedigitalcollection.asme.org/turbomachinery/article-abstract/130/2/021010/451508)). En baja relación de aspecto, pequeños cambios en la cavidad valen **hasta 1 punto de η** ⚠️.

### 4.4 Refrigeración por película y su impacto aerodinámico

- La pérdida aerodinámica es más sensible al soplado en el lado de succión; penalización creciente con tasas de eyección 0–3 % ⚠️.
- **Modelo clásico**: Hartsel, con revisiones modernas ([Research progress of mixing loss model for film cooling, Int. J. Heat Fluid Flow](https://www.sciencedirect.com/science/article/abs/pii/S0142727X24002856)).
- **Truco de coste**: BCs de orificio en la superficie en lugar de mallar plenum+tubos (≈180 % menos malla, ≈300 % menos tiempo) ⚠️.
- **Contabilidad termodinámica**: **Young & Wilcock (2002), "Modelling the Air-Cooled Gas Turbine: Parts 1–2", J. Turbomach. 124(2):207–221**: dividir la expansión por filas, estator y rotor separados; el refrigerante del rotor contribuye al numerador de η. Ver también [Calculating Cooled Turbine Efficiency With Weighted Cooling Flow Distributions, J. Turbomach. 145(6):061007, 2023](https://asmedigitalcollection.asme.org/turbomachinery/article/145/6/061007/1153559).

> ⚠️ **Riesgo de diseño para Phy-AT**: si L0 y el CFD no usan la **misma** definición de rendimiento refrigerado, la calibración afín absorberá una diferencia de definición como si fuera física, y extrapolará mal. Invariante verificado del contrato.

### 4.5 Inestacionariedad estela-álabe y wake recovery en LPT

- Hodson & Howell (2005): la interacción estela-capa límite de succión gobierna la transición en LPT.
- Consecuencia: perfiles **high-lift/ultra-high-lift** (Zweifel alto, menos álabes) funcionan porque la estela periódica suprime la burbuja que existiría en cascada estacionaria. **Un modelo estacionario sobreestima la pérdida de estos perfiles** — sesgo sistemático que Phy-AT debe declarar/calibrar.
- Modelado data-driven del efecto: EARSM entrenados sobre DNS/LES con **gene expression programming** para mezcla de estela en LPT ([ASME GT2018](https://asmedigitalcollection.asme.org/GT/proceedings-abstract/GT2018/51012/V02CT42A009/272402); [Zhao et al., J. Comput. Phys. 2020](https://www.sciencedirect.com/science/article/pii/S002199912030187X)). **Precedente académico del data flywheel, a nivel de cierre.**

### 4.6 Clocking

Efecto real pero **pequeño**: 0.12 % en trabajo y 0.08 % en η en un LPT diseñado para maximizarlo ([Clocking in Low-Pressure Turbines, J. Turbomach. 139(10):101003, 2017](https://asmedigitalcollection.asme.org/turbomachinery/article-abstract/139/10/101003/378793)); 0.086 % en 1.5 etapas ⚠️; +0.44 %/+0.71 % en un caso optimizado de 2 etapas ⚠️.

**Regla para Phy-AT**: 0.08–0.4 % está al nivel del ruido del banco (0.25–0.45 %). **No modelar clocking en L0/L1 y declararlo fuera de alcance.**

### 4.7 Rugosidad y degradación en servicio (frontera 2023–2026)

La rugosidad crece en servicio y aumenta pérdida y flujo de calor. Melbourne + GE Aerospace en Frontier: primer HPT 3D de alta fidelidad representativo de motor con degradación superficial a microescala ([Int. J. Heat Fluid Flow, 2023](https://www.sciencedirect.com/science/article/pii/S0142727X23000334); [OLCF 2026](https://www.olcf.ornl.gov/2026/01/27/frontier-provides-high-fidelity-insights-into-turbine-aerothermal-performance/)). Candidato natural a parámetro de calibración (`k_s`) con pares de máquinas degradadas.

---

## 5. Optimización aerodinámica de turbinas: estado del arte

### 5.1 Métodos adjuntos

- **SU2 (abierto)**: adjunto discreto turbulento multi-escalón con mixing plane conservativo + BCs no reflectantes, AD por sobrecarga de operadores (CoDiPack) ([Multistage Turbomachinery Design Using the Discrete Adjoint Method Within SU2, J. Propulsion and Power, DOI 10.2514/1.B37685](https://arc.aiaa.org/doi/full/10.2514/1.B37685)).
- **HYDRA (RR)**: adjunto multi-objetivo, primer adjunto aeroacústico; Parablading diferenciado algorítmicamente ([Optimization and Engineering, 2019](https://link.springer.com/article/10.1007/s11081-019-09474-x)).
- TRACE (DLR) y elsA (ONERA) tienen ramas adjuntas equivalentes.
- SOTA: **adjunto + caos polinomial adaptativo para optimización robusta** ([Aerospace Science and Technology](https://www.sciencedirect.com/science/article/abs/pii/S1270963823004893)).

### 5.2 Parametrizaciones 3D y ganancias medidas (el techo realista)

| Técnica | Ganancia medida | Fuente |
|---|---|---|
| Endwall no axisimétrico, Trent 500 **HP** | **+0.59 ± 0.25 %** de η de escalón | [J. Turbomach. 125(3):497](https://asmedigitalcollection.asme.org/turbomachinery/article-abstract/125/3/497/459331) |
| Endwall no axisimétrico, Trent 500 **IP** (multi-fila) | **+0.9 ± 0.4 %** | [ASME GT2002-30339](https://asmedigitalcollection.asme.org/GT/proceedings-abstract/GT2002/3610X/119/292264) |
| Endwall, cascada de Durham | −1/3 de la pérdida secundaria | [Hartland/Gregory-Smith](https://www.researchgate.net/publication/30049557_The_benefits_of_turbine_endwall_profiling_in_a_cascade) |
| Endwall, rotor rotativo modelo | **+0.5 a +1.8 %** según incidencia | [Durham repository](https://durham-repository.worktribe.com/output/1360413/) |
| Compound lean + sweep | hasta −14 % de pérdida, +0.93 % η ⚠️ | [AD Technology](https://blog.adtechnology.com/axial-turbine-stacking-best-practices-secondary-flow-suppression) |
| Lean (clásico) | mecanismo: reduce intensidad de flujos secundarios | [Denton & Xu, J. Turbomach. 114(1):184, 1992](https://asmedigitalcollection.asme.org/turbomachinery/article/114/1/184/419987) |

**Lectura para Phy-AT**: el techo de ganancia 3D creíble en turbina es **0.5–1.0 punto de η**. Un optimizador que anuncie +3 puntos está midiendo el error del modelo. **Codificarlo como aviso (guardrail).**

### 5.3 Optimización bayesiana, surrogates y multi-fidelidad

- **Co-kriging AR(1) de Kennedy & O'Hagan (2000)**; formulación recursiva de **Le Gratiet & Garnier** ([arXiv:1210.0686](https://arxiv.org/pdf/1210.0686)); implementación abierta en **SMT (MFCK)** ([docs](https://smt.readthedocs.io/en/v2.9.1/_src_docs/applications/mfck.html)).
- Escala industrial: BO multi-fidelidad multi-objetivo paralela sobre un **compresor axial de 3 escalones con 144 variables** ([Aerosp. Sci. Technol., 2024](https://www.sciencedirect.com/science/article/abs/pii/S1270963824003663)).
- Turbinas: **optimización robusta con transferencia bayesiana sobre álabes ultra-high-lift** ([J. Aerospace Engineering 39(4)](https://ascelibrary.org/doi/10.1061/JAEEEZ.ASENG-6769)) — conceptualmente idéntico a Phy-AT.
- [Comparison of multi-fidelity surrogate models under extreme cost imbalance (2025)](https://link.springer.com/article/10.1186/s40323-025-00316-3).
- Surrogates de pérdida DL: [pérdida de perfil no paramétrica, Energy 2023](https://www.sciencedirect.com/science/article/abs/pii/S0360544223031134); **GNN entrenadas con datos experimentales de pequeña escala para amplio rango de incidencia** ([J. Turbomach. 147(2):021008, 2025](https://asmedigitalcollection.asme.org/turbomachinery/article/147/2/021008/1205175)) — precedente de "pocos pares, mucho prior".
- Revisiones: [AI in turbomachinery aerodynamics, AI Review 2024](https://link.springer.com/article/10.1007/s10462-024-10867-3); [ML Methods in CFD for Turbomachinery, IJTPP 7(2):16](https://www.mdpi.com/2504-186X/7/2/16).

### 5.4 Diseño robusto y UQ

- **Garzon & Darmofal (2003). "Impact of Geometric Variability on Axial Compressor Performance". J. Turbomach. 125(4):692–703** ([ASME](https://asmedigitalcollection.asme.org/turbomachinery/article-abstract/125/4/692/449395)): PCA de mediciones de superficie + Monte Carlo + meanline; **~1 % de η perdido por variabilidad de fabricación representativa**.
- **Montomoli et al., "Uncertainty Quantification in CFD and Aircraft Engines"**, Springer (2015/2019) ([Springer](https://link.springer.com/book/10.1007/978-3-319-92943-9)): capítulos por componente (HP turbine, LP turbine). Tesis: la UQ **indica direcciones de diseño robustas**, no solo barras de error.
- [Adjoint CFD para variaciones de fabricación en un vano de turbina heavy-duty (102 escaneos), arXiv:1901.10352](https://arxiv.org/abs/1901.10352).
- [Mohanamuraly & Müller (2021), adjoint-assisted multilevel multifidelity UQ, IJNME](https://onlinelibrary.wiley.com/doi/10.1002/nme.6617).

---

## 6. Benchmarks y casos test públicos de turbina

| Caso | Tipo | Qué mide | Rango | Dónde | Uso en Phy-AT |
|---|---|---|---|---|---|
| **VKI LS89** (Arts et al.) | NGV transónico 2D, cascada | **transferencia de calor**, velocidad, pérdida | Re₂ 0.5–2.2×10⁶, M₂ 0.7–1.1, Tu 1–6 % | [ASME GT1990](https://asmedigitalcollection.asme.org/GT/proceedings/GT1990/79047/V001T01A106/241577); [AboutFlow](https://www.aboutflow.sems.qmul.ac.uk/events/munich2016/benchmark/testcase2/) | Ancla HPT/aerotérmica; muy sensible a Tu inyectada |
| **SPLEEN C1** (VKI + Safran) | cascada LPT alta velocidad, **con estelas y purga** | pérdidas, presiones, PIV, secundarios inestacionarios | **M₂is 0.70–0.95, Re₂is 65k–120k** | [Zenodo 10.5281/zenodo.8075795](https://zenodo.org/records/8075795); [PIV](https://zenodo.org/records/10253213); [off-design IJTPP 11(1):14](https://doi.org/10.3390/ijtpp11010014) | **El mejor caso público moderno.** Banco de regresión de L2 |
| **T106A/C** | LPT, cascada | presión, pérdida de estela; **referencia DNS** | Re 6×10⁴–2×10⁵ | [Springer](https://link.springer.com/chapter/10.1007/978-3-030-62048-6_17) | Validación de modelos de transición |
| **PAK-B** (P&W/AFRL) | LPT genérico | separación en succión | Re_c 5×10⁴–2×10⁵; Tu 1–4 % | [NTRS 20070010029](https://ntrs.nasa.gov/citations/20070010029) | Curva Yp(Re) |
| **MT1** (QinetiQ → Oxford) | **etapa HP transónica completa** | η, transferencia de calor, **hot streak y swirl** | escala motor | [ORA Oxford](https://ora.ox.ac.uk/objects/uuid:d1df1c3a-9b99-4ea1-95d3-0c826a146cdd) | Único caso público con η de etapa HP + distorsión |
| **Aachen 1.5 etapas** (RWTH) | estator-rotor-estator subsónico | campos entre filas, pérdidas | subsónico | [Tutorial SU2 mixing plane](https://su2code.github.io/tutorials/Aachen_Turbine/) | **Smoke test del pipeline L3** |
| **LISA** (ETH Zurich) | turbina axial 2 escalones **carenada**, torquímetro por escalón | inestacionario, shroud/cavidades | baja velocidad | [ETH LEC](https://lec.ethz.ch/research/turbomachinery_experimental.html) | Shroud/purga/clocking multi-escalón |
| **NASA E³ HPT** (GE, 2 escalones) | HPT + cámara | prestaciones diseño y off-design | escala motor | **[geometría + malla + BCs en data.gov](https://catalog.data.gov/dataset/geometry-grid-and-boundary-condition-data-for-eee-combustor-and-turbine-ef151)** | **Caso semilla "de motor" completo** |
| **FACTOR** (EU FP7) | cámara lean-burn 360° + turbina 1.5 et. | swirl + hot streaks, clocking | representativo motor | [Aerospace 8(10):285](https://www.mdpi.com/2226-4310/8/10/285/htm) | Prior de distorsión de entrada |
| **VKI CT3** | etapa HP transónica, tubo de compresión | flujo de calor y presión inestacionarios | transónico | [VKI CT-3](https://www.vki.ac.be/index.php/research-consulting-mainmenu-107/facilities-other-menu-148/turbomachinery-facilities/90-isentropic-compression-tube-annular-cascade-facility-ct-3) | Espaciado axial y refrigeración |

**Prioridad de adopción**: (1) SPLEEN C1; (2) VKI LS89; (3) Aachen 1.5; (4) NASA E³; (5) MT1; (6) PAK-B/T106.

---

## 7. Cuantificación de incertidumbre y calibración multifidelidad

### 7.1 Los tres marcos que importan

**(a) Calibración bayesiana de Kennedy & O'Hagan (2001)**, JRSS-B 63(3):425–464: `z(x) = ρ·η(x, θ*) + δ(x) + ε` con **δ(x) el término de discrepancia**. Lección crítica: **omitir δ(x) hace que los parámetros calibrados absorban el sesgo estructural y la extrapolación se rompa.**

**(b) Co-kriging AR(1) (Kennedy & O'Hagan, 2000, Biometrika)**: `y_hi = ρ·y_lo + δ(x)`. Recursivo eficiente: Le Gratiet & Garnier ([arXiv:1210.0686](https://arxiv.org/pdf/1210.0686)). **La corrección afín de Phy-AC es el caso degenerado** (ρ=a, δ=b constante, sin UQ). Subir a co-kriging es la evolución natural.

**(c) Multifidelidad para UQ**: **Peherstorfer, Willcox & Gunzburger (2018), SIAM Review 60(3):550–591** ([arXiv:1806.10761](https://arxiv.org/abs/1806.10761)).

### 7.2 Práctica industrial de calibración

1. **Se calibra el CFD contra el banco antes de calibrar nada más** ("calibrated CFD methods" es terminología estándar en capacidad de NGV).
2. **Offsets por familia de máquina** (HPT refrigerada vs LPT no refrigerada), no globales.
3. **Primero capacidad, luego rendimiento** — un error de capacidad re-empareja el escalón y contamina η.
4. **Se declara el envolvente** — fuera del rango muestreado, la corrección se desactiva.

### 7.3 Precisión honesta por fidelidad (síntesis para declarar)

| Fidelidad | η_tt absoluto | Δη comparativo | Capacidad | Fuente del límite |
|---|---|---|---|---|
| Meanline L0 (KO / Craig-Cox) | **±1.5 %** en cargas convencionales | n/a | ±3–5 % | Kacker & Okapuu (1982); hasta 12 % de dispersión entre modelos en turbinas sCO₂ pequeñas ⚠️ |
| Throughflow L1 | ±1–2 pts | ±0.5–1 pt | ±2–3 % | correlaciones + mezcla radial |
| Q3D (MISES-like) | pérdida de perfil ±10–20 % rel. | buena | n/a | fuera de separación gruesa |
| RANS 3D multi-fila | **±1–2 pts** | **±0.3–0.5 pts** | ±1–2 % | Denton (2010); Burdett et al. (2021) |
| RANS en LPT bajo Re | ±2–4 pts | ±0.5–1 pt | — | GT2025-153288 (SPLEEN) |
| URANS/armónico | mejora secundaria ~10 % | ±0.2–0.4 pts | — | ⚠️ |
| LES/DNS | referencia de física, no de η | — | — | Sandberg & Michelassi (2022) |
| **Banco (suelo irreducible)** | **U95 = 0.45 %** | **0.25 %** | ~0.5 % | QinetiQ TTF |

---

## 8. Recomendaciones concretas para Phy-AT

### 8.1 Escalera de fidelidades propuesta

| Nivel | Contenido | Coste objetivo | Qué captura | Qué NO captura |
|---|---|---|---|---|
| **L0** | Meanline stage-stacking. Triángulos, Zweifel, Λ. Pérdidas: **AM revisado + Kacker-Okapuu** default, **Craig-Cox** alternativa, **Coull & Hodson** perfil. Corrección Re por fila. **Refrigeración según Young & Wilcock (2002)**. Guardrail: carta de Smith (1965) | ~0.5–2 ms | η_tt/η_ts, Γ, Λ, α/β, mapa grueso | 3D, endwall, fugas, transición real, hot streaks |
| **L1** | Throughflow SCM axisimétrico (Denton 1978 + mezcla radial) | ~1–5 s | perfiles radiales, matching, capacidad con bloqueo | secundaria 3D real, fuga detallada |
| **L2** ⭐ | (a) motor Q3D de transición (Euler+IBL tipo MISES o surrogate) para el **lapso de Re de LPT**; (b) **corrección afín + residual** con pares (patrón Phy-AC), evolucionable a co-kriging recursivo | ~0.1–10 s | Yp vs Re/M/Tu, transición, high-lift | inestacionariedad, hot streaks, cavidades |
| **L3** | **Lazo CFD externo**: emisión de paquete completo de BCs + retorno de pares, estratificado por (Re, M₂is, ψ, φ, τ/h, purga) | ~10 min (GPU) – 2 h (CPU) por punto | 3D completo | inestacionariedad, LES |
| **L4** (fuera del lazo) | URANS/armónico y LES de referencia; **solo priors** (p.ej. GEP de Sandberg) y casos-testigo | horas–semanas | wake recovery, clocking, rugosidad | — |

### 8.2 Salidas a calibrar (el `KEYS` de HiFiCalibration para turbinas)

```python
KEYS = (
    "Gamma",        # ṁ√T0_in / P0_in — CAPACIDAD. Calibrar SIEMPRE PRIMERO.
    "eta_tt",       # (o eta_ts si la última fila descarga)
    "ER",           # razón de expansión realizada
    "alpha_exit",   # swirl de salida (matching con la fila siguiente)
    "Lambda",       # grado de reacción realizado
    "Yp_row[i]",    # pérdida por fila, si el CFD la desglosa
)
```

**Γ primero**: el aft-loading puede mover el área mínima de control aguas abajo de la garganta geométrica y **alterar la capacidad hasta un 10 %** ([Understanding Capacity Sensitivity of Cooled Transonic NGVs, J. Turbomach. 143(5):051001](https://asmedigitalcollection.asme.org/turbomachinery/article/143/5/051001/1097201)). Calibrar η sin calibrar Γ hace que el ajuste afín absorba un error de matching como pérdida.

**Guardarraíles del ajuste**:
- Rechazar pares con `|a − 1| > 0.15` o `|b|` > 3 puntos de η → error de definición, no física.
- Congelar la corrección fuera de la cáscara convexa de los pares.
- No calibrar por debajo de **0.25 %** de η (ruido de banco).

### 8.3 Condiciones de contorno CFD que Phy-AT debe emitir

El `cfx_boundary_conditions` de Phy-AC (P0, T0, ṁ, RPM, P_out, Tu=0.05, interfaces) **es insuficiente para turbinas**. Paquete propuesto:

**A. Entrada**: `P0_in`, `T0_in` como **perfiles radiales (r, P0, T0)**, no escalares (Denton 2010); `flow_direction` con **swirl residual α(r)** para lean-burn; `hot_streak` opcional (mapa 2D con nº de streaks, amplitud, clocking vs LE del NGV); `turbulence`: **Tu y L_t** + `Tu_target_at_LE` (LS89 es muy sensible; dos usuarios con el mismo Tu_in y distinta L_t obtienen pérdidas distintas).

**B. Salida**: `P_static_out` **con equilibrio radial**; alternativa `mdot_out` para casos ahogados (donde la calibración de Γ es más informativa).

**C. Paredes y térmica** (lo que Phy-AC no tiene): `wall_thermal`: `adiabatic|isothermal|HTC`; si isotermo, **T_wall por fila y superficie** o razón `T0_in/T_wall` (LS89: 1.1 cuasi-adiabático a 1.7 cercano a motor ⚠️ — [RANS Aerothermal Database of LS89, Energies 18:5321](https://doi.org/10.3390/en18195321)); `roughness_ks_m` por superficie.

**D. Refrigeración**: por fila/hilera: `mdot_coolant`, `T0_coolant`, `P0_coolant`, blowing ratio, momentum flux ratio, ángulo, pitch/D. **Método recomendado: BCs de orificio** (no mallar plenum). **Emitir la definición de η usada** (Young & Wilcock 2002) como campo obligatorio `efficiency_definition`.

**E. Fugas y cavidades**: `tip_clearance` por rotor y su variación frío/caliente; `shroud` (geometría de laberinto si carenado); `purge_inlets[]` con ṁ, T0 y **swirl ratio β** — sin él la penalización no es reproducible; sanity check `expected_deta_per_pct_purge = 0.7–1.2` ⚠️.

**F. Interfaces y numérica**: mixing plane con `non_reflecting=True`, ≥32 bandas; **`frozen_rotor_allowed=False`**; `SST` + `gamma-Retheta` **obligatorio si Re_exit < 5e5**; `y_plus_target=1`, ≥25 celdas en capa límite, expansión ≤1.2, ~0.5–1e6 celdas/canal; `gas_model`: cp(T) + **FAR** (productos de combustión, no aire frío γ=1.4).

**G. Contrato de retorno**: `Gamma, eta_tt, eta_ts, ER, alpha_exit(r), Lambda, Yp_row[], mdot_check, T0_out, torque, power, q_wall_avg, y_plus_max, residuals, mesh_count, solver, turbulence_model, transition_model, convergence_flag`. **Sin `y_plus_max`, `mesh_count` y `convergence_flag` el par no se acepta en el flywheel.**

### 8.4 Benchmarks como pares semilla

1. **Aachen 1.5** — smoke test del emisor de BCs (tutorial SU2 público).
2. **SPLEEN C1** — semilla principal de L2 (Yp(Re,M) con estelas y purga).
3. **VKI LS89** — semilla aerotérmica (bloque wall_thermal y Tu/L_t).
4. **PAK-B + T106A/C** — priors del lapso de Reynolds.
5. **NASA E³** — primer caso "de motor" completo (geometría+malla+BC públicas).
6. **MT1** — distorsión de entrada con η medido.
7. **LISA** — shroud/purga multi-escalón con η por escalón.

6 casos × 3–5 puntos ⇒ **~25 pares semilla**, suficiente para arrancar el residual no lineal (umbral 15 de Phy-AC) antes del primer CFD del usuario.

### 8.5 Límites que Phy-AT debe declarar (texto sugerido)

> **Alcance y límites de Phy-AT.**
> 1. Phy-AT **no predice el rendimiento absoluto de una turbina**. Predice tendencias y, dentro del envolvente de sus pares de calibración, valores corregidos (Denton, GT2010-22540: uso comparativo).
> 2. Precisión declarada: **L0 ±1.5–2.5 pts de η_tt y ±3–5 % de capacidad**; **L1 ±1–2 pts**; **L2 calibrado ±0.3–0.5 pts dentro del envolvente**; **L3 ±1–2 pts absolutos, ±0.3–0.5 comparativos**.
> 3. **Suelo de precisión: 0.25 % (relativo) / 0.45 % (absoluto)** — U95 de banco. No se reportan mejoras por debajo como significativas.
> 4. **Fuera de alcance L0/L1/L2**: clocking (0.08–0.4 %, bajo el ruido); wake recovery (⇒ L0/L1 **sobreestiman** la pérdida de perfiles high-lift); migración de hot streaks / Kerrebrock-Mikolajczak; rotura del vórtice de fuga; modos de cavidad de rim seal.
> 5. **Solo por correlación (candidatos prioritarios a calibración)**: purga (0.7–1.2 pts/% ⚠️), fuga de punta (~1/3 de la pérdida), fuga de shroud, penalización de film cooling, lapso de Re en LPT.
> 6. **Techo de ganancia creíble**: 0.5–1.0 punto de η por diseño 3D avanzado. Aviso automático si el optimizador propone más.
> 7. **Invariante de definición**: L0, L1 y el CFD usan la **misma** definición de η refrigerado (Young & Wilcock 2002). Se rechazan pares con definiciones distintas.
> 8. **Variabilidad de fabricación**: no modelada por defecto (~1 pt según Garzon & Darmofal 2003).

### 8.6 Evolución del bloque de calibración (HiFiCalibration)

1. **Ponderación por incertidumbre**: `w=1/σ` del par (CFD ≈ 1 pt, banco ≈ 0.3 pts).
2. **Regularización hacia la identidad**: prior `a~N(1, 0.05²)`, `b~N(0, 0.5²)` — con 2–3 pares el OLS puede dar pendientes absurdas.
3. **Estratificación por régimen**: (HPT refrigerada / LPT no refrigerada) y bandas de Re.
4. **Co-kriging recursivo con banda de incertidumbre** (Le Gratiet) cuando `n_pairs ≥ ~15`.
5. **δ(x) de Kennedy-O'Hagan** — es lo que ya hace el residual del deep ensemble; documentarlo y citarlo como tal.

---

## Fuentes

(ver hipervínculos in-line; principales: Denton GT2010-22540; Sandberg & Michelassi ARFM 54; Hodson & Howell ARFM 37; Tyacke et al. PAS 110; Lee, Dawes & Coull JGPPS 5; Denton MULTALL J.Turbomach 139(12); Brandvik & Pullan J.Turbomach 133(2); Slotnick et al. CFD Vision 2030; GT2025-153288 SPLEEN; Young & Wilcock J.Turbomach 124; Kennedy & O'Hagan 2000/2001; Peherstorfer, Willcox & Gunzburger SIAM Rev 60; Garzon & Darmofal J.Turbomach 125(4); Montomoli et al. Springer; Burdett et al. Proc IMechE A 2021; QinetiQ TTF; casos test: SPLEEN C1 Zenodo, LS89, T106, PAK-B, MT1, Aachen, LISA, NASA E³ data.gov, FACTOR, CT3.)
