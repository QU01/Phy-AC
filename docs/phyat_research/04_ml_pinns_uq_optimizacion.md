# Informe de investigación — ML científico para turbinas axiales: bases para **Phy-AT (Quasar)**

**Fecha:** 19 agosto 2026 · **Alcance:** literatura 2020–2026, con anclajes clásicos donde son load-bearing
**Contexto de destino:** replicar la arquitectura Phy-AC (prior físico L0/L1 → deep ensemble residual con puerta de calidad → NSGA-II con dominancia de Deb + LCB + k-means → calibración afín L2 → data flywheel), en NumPy puro, con presupuesto de 150–500 evaluaciones físicas.

> **Nota metodológica y de honestidad.** El proxy de red de la sesión de investigación bloquea el acceso directo (WebFetch) a `arxiv.org`, `asmedigitalcollection.asme.org`, `link.springer.com`, `mdpi.com`, `zenodo.org`, `sciencedirect.com`, `ntrs.nasa.gov` y `semanticscholar.org`. La investigación se hizo con búsqueda web (que sí devuelve títulos, URLs, DOIs y resúmenes indexados) más una única lectura directa exitosa (GitHub/TurboFlow). **Consecuencia práctica:** los datos bibliográficos (autores/volumen/página) provienen de índices, no de la lectura del PDF en la mayoría de casos. Se marca con ⚠️ las citas cuyos metadatos no se pudieron verificar contra la fuente primaria. Antes de que estas referencias entren en `Quasar_PhyAT_Science.md` conviene resolver cada DOI una vez con red abierta.

---

## 1. PINNs en turbomáquinas (2020–2026)

### 1.1 Qué se ha publicado realmente

**Flujo en cascada (el caso mejor documentado).**
La referencia más citada del dominio es Li, Montomoli y Sharma, *"Investigation of Compressor Cascade Flow Using Physics-Informed Neural Networks with Adaptive Learning Strategy"*, **AIAA Journal 62(4):1400, 2024** ([doi:10.2514/1.J063562](https://arc.aiaa.org/doi/10.2514/1.J063562), preprint [arXiv:2308.04501](https://arxiv.org/abs/2308.04501)). Resuelven el problema **directo e inverso** en una cascada de compresor con pesos adaptativos por término de pérdida y learning rate dinámico para mitigar el desbalance de gradientes; el resultado destacado es la **reconstrucción del campo a partir de vectores de velocidad parciales y presión cerca de pared**, y robustez frente a ruido aleatorio en los datos etiquetados. El propio abstract enmarca el PINN como *"an additional and promising option alongside current dominant CFD methods"* — es decir, complemento, no sustituto.

**Asimilación de datos en flujos transicionales/separados (turbina).**
Hanrahan, Kozul y Sandberg, *"Data Assimilation of Transitional and Separated Turbomachinery Flows With Physics-Informed Neural Networks"*, **J. Turbomach. 147(11):111011, 2025** ([doi:10.1115/1.4068396](https://doi.org/10.1115/1.4068396)). Punto clave: el PINN resuelve las **RANS no cerradas** — no requiere modelo de turbulencia — usando datos dispersos para cerrar el sistema. Esto es exactamente el nicho donde el PINN gana: **problemas inversos donde RANS falla por modelado (no-equilibrio, separación masiva, unsteadiness coherente)**. Precede el trabajo de conferencia GT2023-87103 del mismo grupo (Melbourne), *"Predicting Transitional and Turbulent Flow Around a Turbine Blade With a Physics-Informed Neural Network"*. ⚠️

**Cuántos datos hacen falta (la pregunta incómoda).**
*"The Effect of Training Data on Predicting Turbulent Flow Through a Linear Cascade Using Physics-Informed Neural Networks"*, **J. Turbomach. 148(3):031014** ([enlace ASME](https://asmedigitalcollection.asme.org/turbomachinery/article-abstract/148/3/031014/1222105/The-Effect-of-Training-Data-on-Predicting)). El propio planteamiento del paper admite que *"the specific characteristics of the data — such as quantity and location — required for accurate predictions remain largely uncertain"*. Traducción: **no hay receta**; el PINN de cascada es un método de asimilación cuyo coste de datos no está caracterizado a priori. ⚠️

**Transferencia de calor en álabes refrigerados.**
- PINN para el **problema inverso de transferencia de calor en cavidades rotantes de turbina de gas**, J. Turbomach. 147(7):071010, 2025 ([ASME](https://asmedigitalcollection.asme.org/turbomachinery/article-abstract/147/7/071010/1208627/A-Physics-Informed-Neural-Network-for-Solving-the)). ⚠️
- PINN con transformada de Fourier + atención para **superposición de film cooling multi-fila**, Physics of Fluids 37(6):065174, 2025 ([AIP](https://pubs.aip.org/aip/pof/article-abstract/37/6/065174/3350978/Physics-informed-neural-network-for-predicting)). ⚠️
- *Extended multiphysics-informed neural network for conjugate heat transfer problems*, Int. J. Heat Mass Transfer, 2025 ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0017931025004375)), validado sobre un álabe de turbina refrigerado internamente 2D. ⚠️
- PINN para estimar coeficientes convectivos en **jet impingement** a partir de temperatura dispersa y ruidosa ([arXiv:2507.09356](https://arxiv.org/pdf/2507.09356)).

**Reconstrucción de campos desde experimento.**
*Trade-off between reconstruction accuracy and physical validity in modeling turbomachinery PIV data by Physics-Informed CNN* ([arXiv:2403.00183](https://arxiv.org/pdf/2403.00183)) — el título mismo es el hallazgo: existe un **trade-off explícito** entre ajustar el PIV y satisfacer las ecuaciones.

**Estructuras / aeromecánica** (fuera del alcance aerodinámico de Phy-AT pero relevante para la capa 1s): PINNs para ROM de **blisks con mistuning**, J. Eng. Gas Turbines Power 148(9):091009 ([ASME](https://asmedigitalcollection.asme.org/gasturbinespower/article-abstract/148/9/091009/1231762/Physics-Informed-Neural-Networks-for-Reduced-Order)). ⚠️

**Revisiones del campo.**
- Zou, Xu, Chen, Yao y Fu, *"Application of artificial intelligence in turbomachinery aerodynamics: progresses and challenges"*, **Artificial Intelligence Review 57(8):222, 2024** ([doi:10.1007/s10462-024-10867-3](https://link.springer.com/article/10.1007/s10462-024-10867-3), [dblp](https://dblp.org/rec/journals/air/ZouXCYF24.html)). Es la revisión de referencia del dominio; cubre diseño, validación y O&M.
- *Physics-Informed Neural Networks for Industrial Gas Turbines: Recent Trends, Advancements and Challenges*, [arXiv:2506.19503](https://arxiv.org/pdf/2506.19503) (2025). Su conclusión operativa es explícita: *"traditional CFD and FEM remain the workhorses for detailed turbine analysis; PINNs are emerging as powerful complementary tools particularly for inverse problems and scenarios with limited data"*.
- *Physics-Informed Machine Learning for Intelligent Gas Turbine Digital Twins: A Review*, Energies 18(20):5523, 2025 ([doi:10.3390/en18205523](https://doi.org/10.3390/en18205523)). Taxonomía útil de 4 categorías: modelos termodinámicos aumentados por ANN, arquitecturas operacionales integradas en física, **redes con restricciones físicas (PcNN) como surrogates de CFD**, y enfoques generativos / de descubrimiento de modelos. ⚠️
- *Meta-PINNs* ([arXiv:2603.07740](https://arxiv.org/abs/2603.07740), 2026): meta-aprendizaje para evitar el **reentrenamiento completo bajo cambio de condiciones de contorno** — evaluado en cilindro y en pasaje de cascada de compresor a varios ángulos de ataque. Es el reconocimiento explícito de la limitación #1 de los PINNs para diseño paramétrico.

### 1.2 Limitaciones prácticas documentadas (esto es lo que decide el veredicto)

1. **Modos de fallo de optimización, no de expresividad.** Krishnapriyan et al. (NeurIPS 2021) caracterizaron los *failure modes* de los PINNs: en PDEs rígidas el cuello de botella **no es la capacidad de la red sino la optimización**, con el error concentrándose donde hay gradientes fuertes. Wang, Teng y Perdikaris (SIAM J. Sci. Comput. 2021) atribuyeron fallos al **flujo de gradiente desbalanceado entre términos de pérdida**; Wang et al. (CMAME 2021) al **sesgo espectral** (la red aprende primero las bajas frecuencias). Confirmado y resumido en trabajos recientes ([Dual Cone Gradient Descent, NeurIPS 2024](https://proceedings.neurips.cc/paper_files/paper/2024/file/b2b781badeeb49896c4b324c466ec442-Paper-Conference.pdf); [R3 sampling, PMLR 202](https://proceedings.mlr.press/v202/daw23a/daw23a.pdf)).
2. **Turbulencia.** Sin datos, los PINNs no cierran RANS: se documenta que PINNs **fallan** en un backward-facing step turbulento 2D con k-ε/k-ω salvo que se les inyecte DNS, y que **incluso con asimilación parcial persisten desviaciones significativas** ([PT-PINNs, arXiv:2503.17704](https://arxiv.org/pdf/2503.17704)).
3. **Fenómenos dinámicos.** *Predictive Limitations of PINNs in Vortex Shedding* ([arXiv:2306.00230](https://arxiv.org/pdf/2306.00230)): los PINNs *data-free* no predicen vortex shedding, y los *data-driven* lo exhiben **solo mientras hay datos**, revirtiendo a la solución estacionaria al cortarlos. Extrapolable a estelas y a inestabilidades de estela de rotor.
4. **Choques.** El sesgo espectral y la concentración de error en gradientes fuertes es exactamente lo que hay en el **cuello sónico y la onda de choque en el borde de fuga de una turbina transónica** — el régimen central de una HPT. Las mitigaciones (descomposición de dominio, adaptativos, TSONN, operadores) *"incur substantial computational costs and lack validation for 3D turbulent flows"* (PT-PINNs, ibid.).
5. **Coste vs CFD.** Entrenar un PINN por caso cuesta del orden de una solución CFD (o más), y **no amortiza sobre un barrido paramétrico** salvo con meta-aprendizaje o parametrización explícita — de ahí Meta-PINNs. Una crítica reciente y dura: *"Fundamental flaws of physics-informed neural networks and explainability methods in engineering systems"*, Eng. Appl. AI, 2025 ([ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0360835225008502)). ⚠️

### 1.3 Veredicto honesto sobre PINNs en una escalera de fidelidades de diseño

**Dónde aportan:**
- **Problemas inversos**: reconstruir un campo desde medidas dispersas (PIV, presiones de pared, termopares), estimar HTC/flujo térmico en un álabe refrigerado, inferir condiciones de contorno faltantes. Aquí el PINN gana a CFD porque CFD *no puede resolver el problema planteado*.
- **Asimilación / cierre**: RANS sin modelo de turbulencia con datos escasos (Hanrahan et al. 2025).
- **Post-proceso de validación experimental**, es decir, en el flujo de datos que alimenta la **calibración L2** de Phy-AT: convertir una medida puntual en un campo consistente.

**Dónde NO aportan (para Phy-AT específicamente):**
- **Como sustituto del prior L0/L1.** El meanline stage-stacking de Phy-AC corre en ~0.5 ms/punto y es diferenciable en forma cerrada. Un PINN, tras entrenar horas, sería más lento en inferencia por diseño evaluado, menos fiable en choque y sin garantías de conservación. **No hay ningún trabajo publicado que demuestre un PINN batiendo a un meanline calibrado en coste-precisión para diseño preliminar.**
- **Como surrogate del objetivo dentro del bucle de optimización 15-D.** Con 150–500 evaluaciones, entrenar PINNs paramétricos es inviable; el trabajo de Meta-PINNs existe precisamente porque los PINNs vanilla exigen reentrenamiento por condición.
- **Como "physics loss" añadida al ensemble residual.** Tentador, pero el residuo verdad−L0 **no obedece ninguna PDE conocida**: es un término de discrepancia de modelo. Añadir un residual de Euler/RANS al MLP del ensemble sería física decorativa, no informada. Lo correcto es lo que Phy-AC ya hace: **meter la física por el embedding de features, no por la loss**.

**Recomendación Phy-AT: no implementar PINNs en capas 2–4.** Reservarlos como herramienta *offline* de la capa de datos/validación (L2), y documentarlo así en el `Science.md` para que la decisión sea explícita y defendible en revisión.

---

## 2. Operadores neuronales, GNNs y surrogates de campo

### 2.1 Estado del arte aplicado a turbomáquinas

**FNO y variantes con atención.** *An attention-enhanced Fourier neural operator model for predicting flow fields in turbomachinery cascades* (A-FNO), **Physics of Fluids 37(3):036121, 2025** ([AIP](https://pubs.aip.org/aip/pof/article/37/3/036121/3339141/An-attention-enhanced-Fourier-neural-operator)) — self-attention tipo Galerkin sobre FNO, motivado explícitamente por que los modelos previos fallan **en las regiones de gradiente alto, que son justo donde se generan las pérdidas aerodinámicas**. Complementario: *Prediction of compressor blade cascade flow field based on Fourier neural operator*, Aerospace Science and Technology, 2025 ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S1270963825002792)).

**Transformer neural operator (TNO).** *A panoramic aerodynamic performance prediction method for turbomachinery cascades using transformer-enhanced neural operator*, Chinese Journal of Aeronautics, 2025 ([ScienceDirect](https://www.sciencedirect.com/science/article/pii/S1000936125000792)) — reporta errores menores que FNO y DeepONet en cascadas.

**GNN / MeshGraphNets.** Dos papers ASME muy pertinentes:
- *Predicting Time-Averaged Unsteady Flows in Turbomachinery via Graph Neural Networks*, **J. Turbomach. 148(1):011003** ([ASME](https://asmedigitalcollection.asme.org/turbomachinery/article/148/1/011003/1219581/Predicting-Time-Averaged-Unsteady-Flows-in), [doi:10.1115/1.4069140](https://doi.org/10.1115/1.4069140)). Idea elegante: la GNN aprende **la diferencia entre RANS estacionario y URANS promediado**, operando directamente sobre la malla numérica. **Esto es residual-learning sobre un prior de menor fidelidad — el mismo principio de Phy-AC, aplicado a campos.** ⚠️
- *Prediction of Steady and Unsteady Flow Quantities Using Multiscale Graph Neural Networks*, **J. Turbomach. 147(7):071015, 2025** ([ASME](https://asmedigitalcollection.asme.org/turbomachinery/article/147/7/071015/1209226/Prediction-of-Steady-and-Unsteady-Flow-Quantities)) — GNN multiescala para superar la distancia de información en mallas grandes. ⚠️

**CNN industrial: C(NN)FD.** Bruni, Maleki y Krishnababu (Univ. Lincoln / **Siemens Energy Industrial Turbomachinery**): [arXiv:2306.05889](https://arxiv.org/abs/2306.05889) (framework) y [arXiv:2503.14369](https://arxiv.org/pdf/2503.14369) (*Deep Learning Modelling of Multi-Stage Axial Compressors Aerodynamics*), más [arXiv:2310.04264](https://arxiv.org/html/2310.04264) publicado en **Data-Centric Engineering (Cambridge)**. Predicen en tiempo real el impacto de **variaciones de fabricación y montaje (tip clearance)** sobre campo y rendimiento de compresores axiales multietapa, con precisión comparable al benchmark CFD. Es el ejemplo mejor documentado de surrogate de campo llevado a uso industrial en turbomáquinas.

**Rendimiento de cascadas / álabes con DL.**
- *Application of Deep Learning for Fan Rotor Blade Performance Prediction in Turbomachinery*, J. Turbomach. 147(11):111002, 2025 ([ASME](https://asmedigitalcollection.asme.org/turbomachinery/article/147/11/111002/1215227/Application-of-Deep-Learning-for-Fan-Rotor-Blade)). ⚠️
- *A non-parametric high-resolution prediction method for turbine blade profile loss based on deep learning*, Energy, 2023 ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0360544223031134)) — CNN/ResNet + transfer learning sobre perfiles de turbina; se reportan bibliotecas de **63.450 álabes 2D** con 400 puntos por perfil y sus distribuciones de presión estática.
- *Performance prediction and design optimization of turbine blade profile with deep learning method*, Energy, 2022 ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0360544222012543)).
- *Machine-Learning Models for Loss and Deviation Angle of Compressor Cascades*, **AIAA SciTech 2024-1765** ([AIAA](https://arc.aiaa.org/doi/abs/10.2514/6.2024-1765)) — DNN de desviación y coeficiente de pérdida en función de incidencia, calado, solidez, Mach y Reynolds. **Este es el molde exacto de lo que Phy-AT podría hacer para corregir su modelo de pérdidas de turbina.**

### 2.2 Datasets públicos: qué existe realmente para turbina axial

La respuesta honesta: **no existe un "AirfRANS de turbinas"**, pero hay tres cosas aprovechables y una emergente.

| Recurso | Contenido | Utilidad para Phy-AT |
|---|---|---|
| **SPLEEN C1** (VKI + Safran, H2020 Clean Sky 2) — [Zenodo doi:10.5281/zenodo.7264761](https://zenodo.org/records/10253213), registros [10253213](https://zenodo.org/records/10253213) y [13712403](https://zenodo.org/records/13712403) | Base **experimental** abierta de una cascada LPT de alta velocidad representativa de rotor de LPT moderna: 23 álabes, span 165 mm, M₂ 0.7–0.95, Re 70k–120k, PIV + estelas no estacionarias + purge flow | **Alta.** Es la mejor fuente pública de "verdad" para **anclar/calibrar L2** de un modelo de pérdidas de LPT y para validar el prior en régimen transónico-bajo Re |
| **PLAID / PLAID-datasets** — [arXiv:2505.02974](https://arxiv.org/html/2505.02974v3), [HuggingFace](https://huggingface.co/PLAID-datasets), [docs](https://plaid-lib.readthedocs.io/en/latest/plaid_benchmarks.html) | Estándar de datos + 6 datasets. **Rotor37**: RANS 3D compresible (elsA) del rotor 37 en conducto, con variabilidad de malla y **velocidad de rotación**. **2D_profile**: RANS 2D transónico sobre perfiles tipo álabe con grandes deformaciones | **Media-alta.** Rotor37 es compresor, no turbina, pero es el único benchmark serio de turbomáquina 3D con protocolo ML. Útil como banco de pruebas de metodología |
| **AirfRANS** (Bonnet et al., NeurIPS 2022) | 1.000 perfiles, RANS incompresible subsónico, mallas muy refinadas, escasez deliberada de datos | **Baja-media.** Perfil aislado incompresible ≠ cascada de turbina; sirve para benchmarking metodológico. Ver [arXiv:2504.15993](https://arxiv.org/pdf/2504.15993) para comparativa de arquitecturas |
| **DATED (Zenodo)** — centrífugo | 22M muestras generadas por meanline | No hay análogo axial de turbina equivalente publicado |

**Hallazgo crítico y contraintuitivo del benchmark PLAID:** *"FNO suffers on datasets featuring unstructured meshes with pronounced anisotropies, due to the loss of accuracy introduced by projections to and from regular grids (e.g., **Rotor37** and 2D_profile)"* ([PLAID Benchmarks](https://plaid-lib.readthedocs.io/en/latest/plaid_benchmarks.html)). Es decir: **el FNO, el operador neuronal de moda, es estructuralmente inadecuado para mallas de turbomáquina**, que son anisótropas y no estructuradas por construcción.

Y el contrapunto que más importa para el diseño de Phy-AT: **MMGP** (Casenave et al., **NeurIPS 2023**, [arXiv:2305.12871](https://arxiv.org/pdf/2305.12871), [proceedings](https://proceedings.neurips.cc/paper_files/paper/2023/file/89379d5fc6eb34ff98488202fb52b9d0-Paper-Conference.pdf)) — *Mesh Morphing Gaussian Process* — **no usa GNNs**: morfa mallas a un soporte común, reduce dimensión y ajusta **procesos gaussianos**. Resultado: competitivo o superior a GCNN y MeshGraphNets en Rotor37, **entrenable en CPU**, con **incertidumbre predictiva nativa** y sin necesidad de parametrización de forma. Los autores señalan explícitamente que las GNNs *"depend on extensive datasets and are limited in providing built-in predictive uncertainties"*. Ganó el 1er puesto en la competición NeurIPS asociada ([NeurIPS 2024](https://neurips.cc/virtual/2024/108612)).

**Lectura para Phy-AT:** en el régimen de datos escasos (cientos de puntos) y sin GPU, **los métodos gaussianos/de baja dimensión efectiva baten a los operadores neuronales**, incluso en tareas de campo. Esto refuerza —no debilita— la apuesta de Phy-AC por un modelo pequeño sobre un embedding físico.

### 2.3 Veredicto

Los operadores neuronales y GNNs son **tecnología de capa 1.5** (sustituir/acelerar CFD 3D), no de capas 2–4 de Phy-AT. Requieren O(10³–10⁴) simulaciones de entrenamiento; Phy-AT tiene un presupuesto de 150–500 evaluaciones *totales*. **No adoptar.** El único elemento transferible es **conceptual**: la GNN de J. Turbomach. 148(1):011003 aprende *la diferencia entre dos fidelidades*, validando el residual-learning de Phy-AC en un contexto independiente.

---

## 3. Surrogates multifidelidad

### 3.1 Marco clásico y su aplicación a turbomáquinas

El marco canónico es **Kennedy & O'Hagan (2000), Biometrika 87(1):1–13**, con la descomposición autorregresiva AR(1):

$$f_{hi}(x) = \rho \, f_{lo}(x) + \delta(x)$$

donde ρ es un factor de escala y δ el término de discrepancia, ambos GPs. Su forma práctica para ingeniería es **co-kriging** à la Forrester, Sóbester & Keane (2007), *"Multi-fidelity optimization via surrogate modelling"*, **Proc. R. Soc. A** ([enlace](https://royalsocietypublishing.org/doi/10.1098/rspa.2007.1900)). La ventaja operativa citada consistentemente es que la formulación **provee estimaciones de error de predicción utilizables para buscar puntos de infill**.

Aplicaciones documentadas en turbomáquinas y afines:
- *Multi-Fidelity Surrogate Models for Predicting the Aerodynamic Performance of Turbine Airfoils* / *Aerodynamic Optimization of Turbine Airfoils Using Multi-fidelity Surrogate Models* ([Springer chapter](https://link.springer.com/chapter/10.1007/978-3-319-97773-7_50)). ⚠️
- *Analysis of dataset selection for multi-fidelity surrogates for a turbine problem*, **Struct. Multidiscip. Optim.**, 2018 ([Springer](https://link.springer.com/article/10.1007/s00158-018-2001-8)) — trata la pregunta que Phy-AT tendrá: *cuántos puntos hi-fi y dónde*. ⚠️
- *Comparison of multi-fidelity surrogate models for multi-objective aerodynamic optimization in turbomachinery **under extreme cost imbalance***, Adv. Model. Simul. Eng. Sci., 2025 ([doi:10.1186/s40323-025-00316-3](https://link.springer.com/article/10.1186/s40323-025-00316-3)). El título describe literalmente el régimen de Phy-AT (L0 en ms vs CFD en horas). ⚠️
- *An efficient parallel multi-fidelity multi-objective Bayesian optimization method and application to **3-stage axial compressor with 144 variables***, Aerospace Sci. Tech., 2024 ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S1270963824003663)). ⚠️
- *Multi-fidelity graph neural network for flow field data fusion of turbomachinery*, Energy, 2023 ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0360544223027998)) — LF sin turbulencia, HF con SST k-ω. ⚠️
- *Interdisciplinary design optimization of compressor blades combining low- and high-fidelity models*, Struct. Multidiscip. Optim., 2023 ([doi:10.1007/s00158-023-03516-w](https://link.springer.com/article/10.1007/s00158-023-03516-w)). ⚠️

### 3.2 Multifidelidad con deep learning y transferencia

- **Redes compuestas de Meng & Karniadakis** (JCP 2020): dos subredes, una lineal y una no lineal, para aprender la correlación entre fidelidades — el ancestro de la familia "concatenación de fidelidades". ⚠️ (verificar cita exacta)
- *Multi-fidelity Residual Neural Processes for scalable surrogate modeling*, **ICML 2024** ([ACM](https://dl.acm.org/doi/10.5555/3692070.3693626)) — modela explícitamente **el residuo** entre la agregación de fidelidades bajas y la verdad de la fidelidad alta.
- *Residual multi-fidelity neural network computing*, BIT Numerical Mathematics, 2025 ([Springer](https://link.springer.com/article/10.1007/s10543-025-01058-9)). ⚠️
- Evidencia cuantitativa de ahorro: en un caso de flujo subsuperficial, 2.500 corridas LF + 200 HF dieron **~90% de reducción de coste de simulación** con precisión "casi igual" a un surrogate entrenado solo con HF ([Use of multifidelity training data and transfer learning…, JCP 2022](https://www.sciencedirect.com/science/article/abs/pii/S0021999122008634)).
- El principio general, tal como lo resume la literatura de *delta learning*: *"broad trends are represented by the low-fidelity model, while systematic discrepancies are captured by the residual model using a smaller set of high-fidelity samples"*.

### 3.3 MFBO: BOCA, MF-MES y el estado 2024–2026

- **BOCA** — Kandasamy et al., *Multi-fidelity Bayesian Optimisation with Continuous Approximations*, ICML 2017 ([arXiv:1703.06240](https://arxiv.org/pdf/1703.06240)).
- **MF-MES** — Takeno et al., *Multi-fidelity Bayesian Optimization with Max-value Entropy Search and its parallelization*, ICML 2020 ([PMLR v119](https://proceedings.mlr.press/v119/takeno20a.html), [arXiv:1901.08275](https://arxiv.org/abs/1901.08275)).
- **Revisión**: *Multi-fidelity Bayesian Optimization: A Review*, [arXiv:2311.13050](https://arxiv.org/abs/2311.13050).
- **Buenas prácticas y advertencia**: *Best Practices for Multi-Fidelity Bayesian Optimization in Materials and Molecular Research*, [arXiv:2410.00544](https://arxiv.org/abs/2410.00544) (2024) — da **guías sobre *cuándo* usar MFBO**, reconociendo implícitamente que a veces no compensa.
- **MFBO con restricciones**: [arXiv:2510.10984](https://arxiv.org/pdf/2510.10984) y *Constrained multi-fidelity BO with automatic stop condition* ([arXiv:2503.01126](https://arxiv.org/pdf/2503.01126)).

### 3.4 ¿Qué aporta frente al residual-learning + calibración afín de Phy-AC? (respuesta directa)

Hay que ser preciso sobre qué es cada cosa:

| Componente Phy-AC | Equivalente formal en el marco KOH |
|---|---|
| Ensemble aprende **verdad − L0** | Es exactamente el término de discrepancia **δ(x)**, con ρ ≡ 1 fijado a priori |
| **Calibración afín L2** con pares hi-fi | Es el término **ρ** (más un offset), ajustado *globalmente* y no *como función de x* |

Es decir: **Phy-AC ya es un modelo multifidelidad AR(1), con ρ constante y δ estimado por un deep ensemble en lugar de un GP.** Lo que le falta frente a co-kriging completo es:

1. **ρ(x) dependiente del diseño** — que la escala de la corrección varíe por región del espacio (p. ej. la corrección de η es distinta cerca del choque que lejos).
2. **Propagación consistente de incertidumbre entre niveles** — el co-kriging da σ que incluye la incertidumbre de haber visto pocos puntos hi-fi; la calibración afín de Phy-AC no la propaga.
3. **Selección de fidelidad como decisión de la adquisición** (MFBO) — decidir *si el siguiente punto se evalúa con L0, L1 o CFD*, ponderando por coste.

**¿Cuándo compensa?** Sólo cuando (a) hay **≥3 niveles de fidelidad con costes separados por 1–2 órdenes de magnitud**, (b) el nivel caro se ejecuta **decenas de veces** dentro del bucle, y (c) la correlación LF-HF es alta pero **no afín**. Phy-AT, en su versión inicial, tendrá L0 (ms) + L1 (through-flow, segundos) + un puñado de anclas hi-fi *fuera* del bucle. En ese régimen, **MFBO no compensa**: el overhead de mantener GPs acoplados con ρ e hiperparámetros, con <20 puntos hi-fi, tiene más varianza de estimación que beneficio. Es un caso claro de complejidad prematura.

**Lo que sí compensa ya:** generalizar la calibración afín de **ρ escalar** a **ρ + término lineal en 2–3 features físicas** (p. ej. Mach de salida y coeficiente de carga), es decir, una regresión de discrepancia de bajísima dimensión con regularización fuerte. Es la mejora KOH más barata y de menor riesgo.

---

## 4. UQ calibrada

### 4.1 El campo de juego

- **Deep ensembles** — Lakshminarayanan, Pritzel & Blundell, NeurIPS 2017 ([arXiv:1612.01474](https://arxiv.org/pdf/1612.01474)). La base de Phy-AC. El paper original sugiere M=5.
- **MC-dropout** — Gal & Ghahramani, ICML 2016. Consenso empírico posterior: *"deep ensembles generally outperform MC dropout due to more decorrelated inference models"* y *"significantly outperform MC-dropout in terms of calibration"* (múltiples estudios comparativos; ver [arXiv:2303.16210](https://arxiv.org/pdf/2303.16210) y [NASA NTRS 20230017659](https://ntrs.nasa.gov/api/citations/20230017659/downloads/Unc_Quan_NASA_Final_revised.pdf)). Ovadia et al. (NeurIPS 2019) establecieron que los ensembles son los más robustos bajo *dataset shift* ⚠️ (cita de memoria, verificar).
- **GPs** — *"in low-dimensional datasets, GP is often considered the gold standard for uncertainty quantification"*, y *"Gaussian processes as surrogate models are hard to beat on smaller datasets and optimization budgets, but they scale poorly with amount of data, cannot easily capture non-stationarities and are rather slow at prediction time"* ([Trieste docs, Secondmind](https://secondmind-labs.github.io/trieste/1.0.0/notebooks/deep_ensembles.html)).

### 4.2 ¿Cuántos miembros? (revisión de la elección K=5 de Phy-AC)

Matiz importante y bien documentado: la regla "5 basta" de Lakshminarayanan et al. **no se transfiere limpiamente a regresión**. Estudios de regresión multi-salida encuentran que *"using M = 5 as suggested by Lakshminarayanan et al. (2017) does not guarantee sufficient UQ quality in regression tasks"*, y que la calidad de incertidumbre satura alrededor de **8–12 miembros** ([Towards Reliable Uncertainty Quantification via Deep Ensembles in Multi-output Regression Task, arXiv:2303.16210](https://arxiv.org/pdf/2303.16210); versión revista en Eng. Appl. AI, 2024, [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0952197624000290)). Con MLPs NumPy de (96,96,64) el coste marginal de pasar de 5 a 8–10 miembros es trivial (segundos), así que **es una mejora casi gratis**.

### 4.3 Conformal prediction para regresión en ingeniería (2022–2026) — el hallazgo más accionable

**Fundamentos:** Vovk, Gammerman & Shafer (2005); Papadopoulos et al. (inductive/split CP); Lei, G'Sell, Rinaldo, Tibshirani & Wasserman (2018, JASA) para regresión split-conformal; **Romano, Patterson & Candès (2019, NeurIPS), *Conformalized Quantile Regression* (CQR)**; Angelopoulos & Bates (2021), *A Gentle Introduction to Conformal Prediction*. ⚠️ (clásicos, citados de memoria — verificar).

**La aplicación directa a surrogates de simulación:** Gopakumar et al., *"Uncertainty Quantification of Surrogate Models using Conformal Prediction"*, [arXiv:2408.09881](https://arxiv.org/abs/2408.09881), publicado en **Machine Learning: Science and Technology (IOP)**, [doi:10.1088/2632-2153/ae2e7b](https://iopscience.iop.org/article/10.1088/2632-2153/ae2e7b). Los puntos que importan para Phy-AT, citados casi literalmente del abstract:
- proporciona **cobertura marginal estadísticamente garantizada**, de forma **agnóstica al modelo** y con **coste computacional prácticamente nulo**;
- resuelve la queja de que los métodos previos *"fail to provide statistical guarantees over error bars […] and require computationally expensive ensemble training, extensive sampling, or architectural modifications"*;
- validado sobre **MLP, U-Net, FNO, ViT y GNN**, en PDEs, MHD, meteorología y diagnóstico de fusión;
- la calibración toma **segundos a minutos en hardware estándar**.

Aplicaciones más recientes al dominio: *Multi-Granularity Conformal Prediction for Reliable Neural-Operator Automotive Aerodynamic Surrogates* ([arXiv:2607.17297](https://arxiv.org/html/2607.17297)); *Stratified Conformal Prediction for Neural Fluid Surrogates* ([EngrXiv 7361](https://engrxiv.org/preprint/view/7361)); revisión general *Conformal Prediction: A Data Perspective*, **ACM Computing Surveys** ([doi:10.1145/3736575](https://dl.acm.org/doi/10.1145/3736575)).

**El caveat que hay que documentar y no barrer bajo la alfombra.** CP garantiza cobertura bajo **intercambiabilidad**. El aprendizaje activo la rompe por construcción:
- *"When systems have autonomy to collect their own data, such as in black-box optimization and active learning, their actions induce sequential feedback-loop shifts in the data distribution"* — **feedback covariate shift**: Fannjiang, Bates, Angelopoulos, Listgarten & Jordan, *Conformal prediction under feedback covariate shift for biomolecular design*, **PNAS 2022** ([arXiv:2202.03613](https://arxiv.org/pdf/2202.03613)).
- Tibshirani, Barber, Candès & Ramdas (2019), *Conformal prediction under covariate shift* ([PDF Candès](https://candes.su.domains/publications/downloads/WeightedCP.pdf)) — la reparación es **reponderación por likelihood ratio**, que *"becomes unstable when training and test distributions exhibit limited support overlap"*.
- Barber, Candès, Ramdas & Tibshirani (2023), *Conformal prediction beyond exchangeability*, **Annals of Statistics 51(2)** ([Project Euclid](https://projecteuclid.org/journals/annals-of-statistics/volume-51/issue-2/Conformal-prediction-beyond-exchangeability/10.1214/23-AOS2276.pdf)) — cotas de degradación de cobertura cuando la intercambiabilidad falla.

**CP dentro de optimización bayesiana** (existe y funciona):
- Stanton, Maddox & Wilson, *Bayesian Optimization with Conformal Prediction Sets*, **AISTATS 2023** ([PMLR v206](https://proceedings.mlr.press/v206/stanton23a.html), [arXiv:2210.12496](https://arxiv.org/abs/2210.12496)) — dirige las consultas hacia regiones donde las predicciones tienen validez garantizada.
- Deshpande, Marx & Kuleshov, *Online Calibrated and Conformal Prediction Improves Bayesian Optimization* ([PMC11482741](https://pmc.ncbi.nlm.nih.gov/articles/PMC11482741/)) — la calibración *online* evita el problema de intercambiabilidad usando garantías adversariales en lugar de i.i.d.
- *Robust Bayesian Optimization via Localized Online Conformal Prediction* ([arXiv:2411.17387](https://arxiv.org/pdf/2411.17387)).

### 4.4 Recalibración paramétrica (alternativa/complemento barato)

**Kuleshov, Fenner & Ermon (2018), ICML, *Accurate Uncertainties for Deep Learning Using Calibrated Regression***: aprende por **regresión isotónica** un mapeo monótono entre niveles de cuantil esperados y observados, logrando uniformidad del PIT. Es el baseline estándar para recalibrar regresores ([slides Stanford](https://ai.stanford.edu/~kuleshov/papers/uai2018-slides.pdf)). Ventaja para Phy-AC/AT: **la regresión isotónica es ~30 líneas de NumPy (PAVA)** y no necesita dependencias.

Complementos: *Stratification of uncertainties recalibrated by isotonic regression* ([arXiv:2306.05180](https://arxiv.org/pdf/2306.05180)); *Distribution-Free Model-Agnostic Regression Calibration via Nonparametric Methods* ([arXiv:2305.12283](https://arxiv.org/pdf/2305.12283)); *Model-Free Local Recalibration of Neural Networks* ([arXiv:2403.05756](https://arxiv.org/pdf/2403.05756)).

### 4.5 GP vs ensembles en ~15-D con ~500 puntos — recomendación honesta

Los hechos relevantes:
- 500 puntos en 15-D es **régimen de datos escasos** (≈ 33 puntos/dimensión, muy por debajo de las reglas empíricas de 6·n_d a 24·n_d **por cada** región de interés).
- Un GP anisótropo (ARD) con 500 puntos es computacionalmente trivial (Cholesky de 500×500) y **es donde el GP es "gold standard"**. Los 15 hiperparámetros de longitud de escala son estimables con 500 puntos, aunque con varianza notable.
- Pero: la superficie de respuesta de Phy-AC/AT es **discontinua en el choque y escalonada por n_stages entero** — exactamente las **no estacionariedades** que el GP estándar (kernel estacionario) *"cannot easily capture"*.
- Y: Phy-AC **no aprende la función, aprende el residuo sobre L0**. El residuo es más suave y de menor amplitud que la función — lo cual favorece a **ambos** métodos, pero elimina el argumento de que "el GP es imprescindible porque captura mejor la estructura".
- Y el constraint duro: **NumPy puro**. Un GP con ARD + optimización de hiperparámetros por L-BFGS en NumPy sin autodiff es implementable pero introduce fragilidad numérica (jitter, mal condicionamiento) que un ensemble de MLPs no tiene.

**Recomendación:** **mantener el deep ensemble**, y arreglar su punto débil (calibración) con las dos herramientas baratas de arriba (conformal split + recalibración isotónica) en lugar de cambiar de familia de modelo. Un GP sería una alternativa defendible pero **no una mejora clara**, y costaría el invariante de simplicidad del proyecto. Si algún día se quiere probar GP, el punto de entrada correcto es la evidencia de MMGP (§2.2): GP sobre representación reducida gana en datos escasos.

---

## 5. Optimización

### 5.1 ¿Sigue NSGA-II siendo razonable en 2026 para 15-D multiobjetivo con restricciones?

**Respuesta: sí, en el papel exacto que juega en Phy-AC — como buscador *sobre el surrogate/prior barato*, no como consumidor del presupuesto de evaluaciones caras.** Hay que separar dos preguntas que la literatura suele mezclar.

**(a) NSGA-II como optimizador de presupuesto caro.** Aquí la evidencia le es desfavorable:
- *"Bayesian Optimization Algorithms have shown convergence improvements of between 5.9% and 31.9% over NSGA-II on problems with a budget of 500 design evaluations"* ⚠️ (métrica reportada en índice; verificar fuente primaria).
- MORBO (Daulton et al., ICML 2022, [arXiv:2109.10964](https://arxiv.org/pdf/2109.10964)): en el problema Mazda (222 dimensiones, 54 restricciones black-box), *"while NSGA-II made progress from the initial feasible solution, it is not competitive with MORBO"*.
- Benchmark en diseño estructural: BO *"performed considerably better in terms of rate-of-improvement, final solution quality, and variance across repeated runs"* ([Multi-objective constrained Bayesian optimization for structural design, SMO 2020](https://link.springer.com/article/10.1007/s00158-020-02720-2)). ⚠️

**(b) NSGA-II como buscador interno sobre un modelo barato.** Aquí sigue siendo la elección estándar y sensata: *"NSGA-II […] combines fast non-dominated sorting with a crowding-distance diversity preserver and **has become the default benchmark** for multi-objective evolutionary algorithms"* (Deb, Pratap, Agarwal & Meyarivan, **IEEE Trans. Evol. Comput. 6(2):182–197, 2002**). Y sobre manejo de restricciones: *"NSGA-II with CDP [constraint domination principle] provided better uniform spread"* frente al método de penalización — respalda directamente la elección de Phy-AC de usar **dominancia restringida de Deb sobre las g exactas del prior**.

**El matiz decisivo que salva a Phy-AC:** el paper de MORBO admite que *"in some expensive design problems, state-of-the-art methods such as qNEHVI do not outperform NSGA-II"*. Y en Phy-AC, NSGA-II **no gasta evaluaciones caras**: corre 96×60 ≈ 5.760 evaluaciones del prior barato + surrogate, y el LCB + k-means es quien decide las 14 evaluaciones caras. **Esa arquitectura ya es un SAEA**, y el papel de la adquisición LCB es el papel de la BO. La crítica "NSGA-II es peor que BO" no aplica al uso que Phy-AC le da.

### 5.2 Alternativas y lo que aportarían

| Método | Referencia | ¿Aporta a Phy-AT? |
|---|---|---|
| **qNEHVI** | Daulton, Balandat & Bakshy, NeurIPS 2021 — *"qEHVI is the current state-of-the-art for batch multi-objective optimization"* ⚠️ | **Conceptualmente sí, en implementación no.** Requiere GPs + BoTorch/PyTorch. Pero su *idea* — maximizar mejora de **hipervolumen** en lote, no LCB por objetivo — es transferible a NumPy |
| **TuRBO / SCBO** | Eriksson et al., NeurIPS 2019; Eriksson & Poloczek, AISTATS 2021 (SCBO = versión con restricciones escalable). Avances 2024–2026: [REI](https://arxiv.org/pdf/2412.11456), [AdaScale-TuRBO](https://arxiv.org/abs/2604.22967), [MG-TuRBO](https://arxiv.org/pdf/2604.08569) | **Región de confianza: idea barata y valiosa.** Restringir la búsqueda a una caja alrededor del mejor punto factible, encogiéndola tras fracasos, es ~40 líneas y mitiga el fallo típico de surrogates en alta dimensión |
| **CMA-ES con restricciones** | Hansen | No aporta: monoobjetivo por naturaleza; escalarizar destruiría el frente Pareto que Phy-AC entrega |
| **SAEA** | Revisión: Expert Systems with Applications, 2023 ([doi:10.1016/j.eswa.2022.119495](https://dl.acm.org/doi/10.1016/j.eswa.2022.119495)). Turbomáquinas: *An efficient surrogate-assisted differential evolution algorithm for turbomachinery cascades optimization with **more than 100 variables***, Aerospace Sci. Tech., 2023 ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S1270963823005710)) | **Phy-AC ya es un SAEA.** Advertencia del campo: *"the over-reliance on the accuracy of the surrogate model causes optimization performance of most SAEAs to decrease drastically with increasing dimensionality"* — **la puerta de calidad de Phy-AC es precisamente la mitigación de esto, y es una decisión de diseño acertada y poco común** |
| **Híbrido EA+BO** | **EGBO**: *Evolution-guided Bayesian optimization for constrained multi-objective optimization*, **npj Computational Materials, 2024** ([doi:10.1038/s41524-024-01274-x](https://www.nature.com/articles/s41524-024-01274-x)) — presión de selección evolutiva **en paralelo** con qNEHVI, *"achieving better coverage of the PF and limiting sampling in the infeasible space"* | **Es la validación externa del patrón de Phy-AC**: EA + adquisición bayesiana juntos, no uno u otro |
| Aplicación turbomáquina | *Integrated Surrogate Model-Based Approach for Aerodynamic Design Optimization of Three-Stage Axial Compressor*, Energies 18(17):4514, 2025 ([doi:10.3390/en18174514](https://doi.org/10.3390/en18174514)) — SAEA + reducción de dimensionalidad ⚠️ | Confirma que el patrón surrogate+EA es el estándar industrial en compresores multietapa |

### 5.3 Diseño generativo (GANs / VAEs / difusión) para álabes de turbina

Explosión real en 2024–2026, y es **la parte más "hype-sensible"** del informe.

**Difusión — lo más relevante:**
- *Generative Inverse Aerodynamic Design of a Single-Stage Turbine Using Conditional Denoising Diffusion Probabilistic Model*, **J. Turbomach. 148(4):041001** ([ASME](https://asmedigitalcollection.asme.org/turbomachinery/article/148/4/041001/1222148/Generative-Inverse-Aerodynamic-Design-of-a-Single)). Toma prestaciones aerodinámicas + indicadores clave como entrada y **genera directamente los parámetros de diseño**. Reporta *"maximum relative error of less than 1.5% compared to target performance metrics"* ⚠️ — cifra a verificar, y en todo caso medida contra el propio modelo de entrenamiento, no contra CFD independiente.
- *A New Paradigm for 3D Turbomachinery Design: Generative Diffusion Model Based Framework with **Direct Geometry Encoding***, [arXiv:2607.27093](https://arxiv.org/abs/2607.27093) (2026) — entrenado directamente sobre geometrías 3D **sin parametrización**.
- *Data-Driven Inverse Design of Turbine Blade Passages*, Energies 19(12):2796, 2026 ([MDPI](https://www.mdpi.com/1996-1073/19/12/2796)).
- *An inverse aerodynamic design framework for compressor blades based on generative model*, Phys. Fluids 37(7):076108, 2025 ([AIP](https://pubs.aip.org/aip/pof/article-abstract/37/7/076108/3351690/An-inverse-aerodynamic-design-framework-for)).
- *Dflow-SUR: Enhancing Generative Aerodynamic Inverse Design using Differentiation Throughout Flow Matching* ([arXiv:2512.08336](https://arxiv.org/pdf/2512.08336)).

**VAE/GAN como reducción de dimensionalidad:**
- *Knowledge transfer accelerated turbine blade optimization via a sample-weighted variational autoencoder*, Aerospace Sci. Tech., 2024 ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S1270963824001317)).
- **Bézier-GAN** (Chen et al.): selecciona 3 códigos latentes de 13 parámetros capturando la variabilidad principal del espacio de perfiles ⚠️.
- DFFD-VAEGAN ([GPPS-TC-2025 paper 179](https://gpps.global/wp-content/uploads/2025/09/GPPS-TC-2025_paper_179.pdf)) — con el diagnóstico honesto del campo: *"GANs often collapse into mode failure while VAEs struggle to balance generation quality with latent space continuity"*.

**El contrapeso crítico:** *Benchmarking Generative AI Against Bayesian Optimization for Constrained Multi-Objective Inverse Design* ([arXiv:2511.00070](https://arxiv.org/pdf/2511.00070)). **Este es el paper que hay que leer antes de dejarse llevar por la ola generativa.**

**Veredicto:** los modelos generativos requieren **miles de geometrías etiquetadas** (63.450 álabes en el estudio de Energy 2023). Phy-AT no las tendrá hasta que el data flywheel haya girado mucho. Además, **generan geometría, y el problema de Phy-AT en capas 2–4 es de 15 variables meanline, no de forma libre**. Su lugar natural es la **capa 5 (geometría)** y en un futuro lejano.

**Nota tangencial pero relevante para el posicionamiento de Quasar:** ya existe *TurboAgent: An LLM-Driven Autonomous Multi-Agent Framework for Turbomachinery Aerodynamic Design* ([arXiv:2604.06747](https://arxiv.org/pdf/2604.06747), 2026), y hay un cuerpo creciente de **DRL para diseño de álabes**: *Turbine blade optimization considering smoothness of the Mach number using deep reinforcement learning* (Information Sciences, 2023, [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0020025523006515)) y *Automation of the CFD-based design process for turbomachinery blades using deep reinforcement learning* (Computers & Fluids, 2026, [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0045793026002057)), con agentes especializados para malla, convergencia CFD y diseño multi-condición. **DRL no es viable para Phy-AT** (necesita 10⁴–10⁶ interacciones; el presupuesto es 500).

---

## 6. Meanline / loss models de turbina aumentados por ML

Esta sección es, junto con la §4.3, **la que más valor específico tiene para Phy-AT**, porque toca la capa 1 (el prior) y por tanto determina la calidad de todo lo demás.

### 6.1 El anclaje clásico que Phy-AT debe implementar en L0

La cadena canónica de modelos de pérdida de turbina axial (que la literatura ML usa siempre como baseline): **Ainley & Mathieson (1951)** → revisado por **Dunham & Came (1970)** → **Kacker & Okapuu (1982)** ("modelo KO"), con la alternativa **Craig & Cox (1971)** y la corrección de pérdida secundaria de **Benner, Sjolander & Moustapha (2006)**, que *"took the influence area of the secondary flow into account to improve the KO model"*. Correlaciones de perfil y desviación mejoradas también en la línea de *Improved Profile Loss and Deviation Correlations for Axial-Turbine Blade Rows* ([ResearchGate](https://www.researchgate.net/publication/267499762_Improved_Profile_Loss_and_Deviation_Correlations_for_Axial-Turbine_Blade_Rows)). Y como marco de diseño: el **diagrama de Smith (1965)** (η vs ψ, φ) y el **criterio de Zweifel (1945)** para solidez óptima.

**Implementación open-source de referencia (MIT):** **TurboFlow** — Anderson, Agromayor, Haglind & Nord, **JOSS 10(111):7588, 2025** ([JOSS](https://joss.theoj.org/papers/10.21105/joss.07588), [GitHub](https://github.com/turbo-sim/turboflow), [docs](https://turbo-sim.github.io/turboflow/)). Verificado por lectura directa del repo: Python, **licencia MIT**, formulación **equation-oriented** compatible con optimización basada en gradiente, propiedades de gas real vía CoolProp, y **submodelos intercambiables de pérdida, desviación y *choking***. Papers asociados: *Equation-Oriented Meanline Method for Axial Turbine Performance Prediction Under **Choking Conditions***, J. Turbomach. 147(4):041002, 2025 ([ASME](https://asmedigitalcollection.asme.org/turbomachinery/article-abstract/147/4/041002/1206960/Equation-Oriented-Meanline-Method-for-Axial)) y *Equation-Oriented and Black-Box Design Optimization of Axial Turbines Using Gradient-Based and Gradient-Free Algorithms*, J. Turbomach. 148(1):011011 ([ASME](https://asmedigitalcollection.asme.org/turbomachinery/article/148/1/011011/1221042/Equation-Oriented-and-Black-Box-Design)).

> **Lección de arquitectura para Phy-AT:** el tratamiento del **choking** en turbinas es un problema de primer orden (a diferencia del compresor, la turbina *opera* rutinariamente con la garganta crítica), y hay un paper ASME dedicado solo a formularlo de modo continuo y diferenciable. Phy-AC ya aprendió esta lección en compresor (penalización continua de choque plegada en g, verificada a través de la frontera). **En turbina el requisito es más duro: el choque no es un fallo, es un modo de operación normal, y g debe distinguir "gargantado por diseño" de "gargantado por error".**

### 6.2 ML aumentando modelos de pérdida (la evidencia)

**El paper de referencia — y el más filosóficamente cercano a Phy-AT:**
**Senior, A.C. & Miller, R.J. (Whittle Lab, Cambridge), *"A Data-Centric Approach to Loss Mechanisms"*, J. Turbomach. 146(4):041007, 2024** (paper TURBO-23-1233, online dic. 2023; versión de conferencia GT2023-87080; [ASME](https://asmedigitalcollection.asme.org/turbomachinery/article/146/4/041007/1171656/A-Data-Centric-Approach-to-Loss-Mechanisms), [repositorio Apollo](https://www.repository.cam.ac.uk/items/71f06fd3-59d5-4bb6-8dcf-6525b7b8a090)).
Tesis del paper: descomponer la pérdida total en **modelos físicos de bajo orden** es la forma potente de construir modelos de pérdida, pero *"in complex flows it is often not clear how to break a flow down physically without making assumptions […] which often leads to loss models of low accuracy that only work in a limited part of the design space"*. Su propuesta: **usar ML para *aumentar al diseñador* en el descubrimiento de la descomposición física correcta**, no para reemplazarla. Esto es, textualmente, la filosofía de Phy-AC/Phy-AT.

**Otras aportaciones concretas:**
- *Knowledge enhanced modeling of low-pressure turbine profile loss by combining physical-based and data-driven methods*, **Energy, 2025** ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S036054422500756X)). Modelo híbrido de pérdida de perfil de LPT que **considera intensidad de turbulencia de corriente libre y estelas periódicas**, distingue los modos separación-readherencia y separación-sin-readherencia, y predice pérdida en diseño y fuera de diseño **más las posiciones de transición y readherencia**. Reivindica explícitamente *"accuracy, generalization and interpretability"*.
- *Development of helium turbine loss model based on knowledge transfer with Neural Network and its application on aerodynamic design*, **Energy, 2024** ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0360544224011009), preprint [arXiv:2309.06709](https://arxiv.org/pdf/2309.06709)) — **transfer learning entre fluidos de trabajo** para el modelo de pérdidas. Muy relevante si Phy-AT quiere cubrir turbinas de gas y de vapor/ORC con un solo prior.
- *A Spanwise Loss Model of Turbine Cascade with Tip Clearance Based on Machine Learning* ([Springer, 2024](https://link.springer.com/chapter/10.1007/978-981-97-3998-1_90)) — el efecto del *tip leakage* sobre la distribución **radial** de pérdida, aprendido con ML. Directamente aplicable a la capa L1 (through-flow) de Phy-AT.
- *Machine-Learning Models for Loss and Deviation Angle of Compressor Cascades*, **AIAA SciTech 2024-1765** ([AIAA](https://arc.aiaa.org/doi/abs/10.2514/6.2024-1765)) — DNN entrenada con CFD sobre (incidencia, calado, solidez, Mach, Reynolds). Es la plantilla de features exacta para un modelo de pérdidas de cascada de turbina.
- *A Hybrid Data-Driven Adaptive Correction Model for **Axial Compressor Meanline** Performance Prediction*, JMSE 14(9):825 ([doi:10.3390/jmse14090825](https://doi.org/10.3390/jmse14090825)) — meanline + **factores de ganancia** que corrigen desviación y pérdida de presión total mediante corrección adaptativa data-driven. **Esta es exactamente la calibración L2 de Phy-AC, publicada de forma independiente.** ⚠️
- Contexto amplio: *RANS Turbulence Model Development using CFD-Driven Machine Learning* ([arXiv:1902.09075](https://arxiv.org/pdf/1902.09075)) y *Accelerating CFD-driven training of transition and turbulence models for turbine flows* (Computers & Fluids, 2025, [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0045793025003871)).

### 6.3 Lectura para Phy-AT

El campo ha convergido a una conclusión que valida la arquitectura Phy-AC y que Phy-AT debe heredar literalmente: **la combinación ganadora es física de bajo orden interpretable + corrección ML de la discrepancia, con las features físicas correctas como entrada.** No modelos ML puros, no física pura. Senior & Miller (2024) y el paper de LPT de Energy (2025) lo dicen con esas palabras.

---

## 7. Aprendizaje activo: estrategias de adquisición para presupuestos pequeños con restricciones

### 7.1 Manejo de restricciones en la adquisición

- Base clásica: **EIc / cEI** — *"Expected Constrained Improvement is a popular acquisition function for constrained BO, where the expected improvement is weighted by the probability that a point is a feasible design"*. Orígenes: Schonlau et al. (1998); Gardner et al. (ICML 2014); Gelbart, Snoek & Adams (UAI 2014, *BO with unknown constraints*). ⚠️
- Revisión: *Constrained Bayesian Optimization: A Review* ([UHasselt](https://documentserver.uhasselt.be/bitstream/1942/45255/1/Constrained%20Bayesian%20Optimization_%20A%20Review.pdf)).
- **La crítica de 2025 que importa:** *"Existing approaches modify EI by incorporating feasibility probabilities, **requiring an initial feasible point**, and often **restricting exploration to the feasible region**"* — de ahí la propuesta de un criterio *feasibility-infeasibility weighted improvement* que explora deliberadamente a ambos lados de la frontera: **J. Comput. Graph. Stat.**, [doi:10.1080/10618600.2025.2611111](https://www.tandfonline.com/doi/abs/10.1080/10618600.2025.2611111).
- Multiobjetivo con restricciones: *Constrained Multi-objective Bayesian Optimization through Optimistic Constraints Estimation* ([arXiv:2411.03641](https://arxiv.org/pdf/2411.03641)); *Constrained BO with Adaptive Active Learning of Unknown Constraints* ([arXiv:2310.08751](https://arxiv.org/pdf/2310.08751)); *Constrained BO with merit functions* ([arXiv:2403.13140](https://arxiv.org/pdf/2403.13140)).

**Ventaja estructural de Phy-AC/Phy-AT que conviene enunciar explícitamente:** en Phy-AC **las restricciones g NO son black-box** — se calculan exactas y baratas desde el prior. Toda la literatura de "restricciones desconocidas modeladas con un GP auxiliar" es **inaplicable, y eso es una ventaja, no una carencia**. Deb-CDP sobre las g exactas es estrictamente mejor que estimar la factibilidad con un surrogate. Documentarlo evita que un revisor lo lea como una omisión.

### 7.2 Adquisición por lotes

- **k-means++ como aproximación greedy a un DPP**: *"K-means++ can be viewed as a computationally efficient, greedy approximation to sampling from a DPP where the similarity between items is inversely related to the distance between their gradient embeddings"* — es el argumento de **BADGE** (Ash et al., ICLR 2020) ⚠️, y **la justificación teórica retroactiva del k-means de Phy-AC**. Vale la pena citarlo en el `Science.md`: la de-duplicación por k-means no es un truco ad hoc, es un DPP barato.
- DPP explícito: Bıyık et al., *Batch Active Learning Using Determinantal Point Processes* ([arXiv:1906.07975](https://arxiv.org/abs/1906.07975)) — más caro computacionalmente, mejora marginal.
- *Scalable Batch Acquisition for Deep Bayesian Active Learning* ([arXiv:2301.05490](https://arxiv.org/pdf/2301.05490)).
- Alternativa clásica: **local penalization** (González et al., AISTATS 2016) para lotes con GPs. ⚠️

### 7.3 Qué mejoraría sobre el LCB + k-means de Phy-AC

Diagnóstico honesto del esquema actual: LCB por objetivo + explotación del Pareto predicho + exploración de σ alta + k-means es **razonable, robusto y barato**, pero tiene tres debilidades identificables:

1. **LCB por objetivo no es lo mismo que mejora del frente.** Dos puntos con LCB excelente pueden estar pegados en el espacio objetivo y aportar hipervolumen casi nulo. La corrección canónica es **EHVI/qNEHVI** (Daulton et al. 2021). Versión implementable en NumPy puro: **calcular la contribución de hipervolumen esperada por Monte Carlo** — muestrear del ensemble (μ, σ) por miembro, calcular ΔHV respecto al frente actual, promediar. Con 2 objetivos el hipervolumen 2-D es un cálculo O(n log n) trivial.
2. **No hay presión de factibilidad explícita en el lote.** Con ~20% de región factible, un lote de 14 puede desperdiciar varios en infactibles. Mitigación barata: ponderar la adquisición por un **margen de factibilidad suave** sobre las g exactas — y, siguiendo la crítica de JCGS 2025, **no** excluir los infactibles marginales, sino sesgar hacia la frontera (donde suelen vivir los óptimos restringidos).
3. **Exploración global sin región de confianza.** En 15-D con un ensemble entrenado con 320+ puntos, la exploración de σ alta tiende a los rincones del hipercubo, donde el prior también es menos fiable. La corrección barata es **TuRBO-lite**: una caja de confianza centrada en el mejor factible, que se encoge tras rondas sin mejora y se expande tras éxitos (Eriksson et al. 2019; SCBO para la versión con restricciones).

---

# Recomendaciones concretas para Phy-AT (capas 2–4)

## A. Conservar de Phy-AC tal cual — y por qué (con evidencia)

| Elemento | Conservar | Justificación con literatura |
|---|---|---|
| **Residual-learning (verdad − L0)** | **Sí, sin cambios** | Es AR(1)/KOH con ρ=1 (Kennedy & O'Hagan 2000; Forrester et al. 2007). Validado independientemente en turbomáquinas por la GNN que aprende URANS−RANS (J. Turbomach. 148(1):011003) y por Multi-fidelity Residual Neural Processes (ICML 2024). Su virtud clave: **el surrogate degrada hacia la física, no hacia el ruido** |
| **Embedding físico en lugar de física en la loss** | **Sí** | Los PINNs meten la física en la loss y pagan sesgo espectral, patologías de gradiente y fallo en choque (Krishnapriyan 2021; Wang 2021). El residuo verdad−L0 no obedece ninguna PDE: meterle un residual de Euler sería física decorativa |
| **Deep ensemble de MLPs pequeños en NumPy** | **Sí** | Deep ensembles > MC-dropout en calibración (consenso 2019–2025). GP sería alternativa, no mejora, dada la no estacionariedad por choque y n_stages entero — y MMGP (NeurIPS 2023) muestra que la ventaja del GP viene de la **reducción de dimensión**, que aquí ya la da el embedding físico |
| **Puerta de calidad antes de guiar la búsqueda** | **Sí — y presumir de ella** | La revisión de SAEAs identifica *"over-reliance on the accuracy of the surrogate model"* como **la** causa de degradación en alta dimensión ([ESWA 2023](https://dl.acm.org/doi/10.1016/j.eswa.2022.119495)). La puerta es una mitigación explícita y poco frecuente en la literatura. **Es un diferenciador defendible de Quasar** |
| **NSGA-II + dominancia restringida de Deb sobre g exactas** | **Sí** | Deb et al. (2002) sigue siendo el benchmark por defecto; CDP > penalización en spread uniforme. Las críticas "BO bate a NSGA-II" aplican a NSGA-II **consumiendo presupuesto caro**, no a NSGA-II sobre el prior barato. EGBO (npj Comput. Mater. 2024) es evidencia externa de que EA+adquisición bayesiana juntos es la combinación correcta |
| **g exactas y no modeladas** | **Sí** | Toda la literatura de constrained-BO existe porque las restricciones suelen ser black-box. Aquí no lo son. Es una ventaja estructural — **documentarla explícitamente** |
| **k-means para de-duplicación del lote** | **Sí** | Es una aproximación greedy barata a un DPP (argumento de BADGE). El DPP explícito (Bıyık et al. 2019) es más caro con mejora marginal |
| **Calibración afín L2 con pares hi-fi** | **Sí, como base** | Es el término ρ de KOH. Independientemente reinventada como "gain factors" en meanline de compresor (JMSE 14(9):825). Ver mejora M4 |
| **Data flywheel (dataset.csv por run)** | **Sí** | Es el activo estratégico. **No existe un dataset público de turbina axial comparable a AirfRANS**; el corpus propio es la ventaja competitiva real |

## B. Mejoras puntuales a adoptar, con evidencia y coste

Ordenadas por **(valor / riesgo)** decreciente.

**M1 — Conformal prediction split sobre el ensemble (ALTA prioridad, ~50 líneas NumPy).**
Reservar un *calibration set* (10–15% del dataset), calcular puntuaciones de no-conformidad normalizadas `s_i = |y_i − μ(x_i)| / σ(x_i)`, y tomar el cuantil ⌈(n+1)(1−α)⌉/n. Da un intervalo `μ ± q_α·σ` con **cobertura marginal garantizada**, agnóstico al modelo, con **coste computacional prácticamente nulo** (Gopakumar et al., *Mach. Learn.: Sci. Technol.*, [doi:10.1088/2632-2153/ae2e7b](https://iopscience.iop.org/article/10.1088/2632-2153/ae2e7b); [arXiv:2408.09881](https://arxiv.org/abs/2408.09881)).
*Impacto sobre Phy-AT:* la **puerta de calidad deja de ser heurística**. Hoy exige "cobertura ±2σ ∈ [0.80, 1]" — un umbral empírico. Con conformal, el factor `q_α` se *deriva* de los datos y la cobertura es una garantía, no una esperanza. Y el σ que entra en el LCB queda escalado correctamente.
*Caveat obligatorio a documentar:* la garantía es marginal y bajo intercambiabilidad. **El aprendizaje activo la viola** (feedback covariate shift; Fannjiang et al., PNAS 2022; Barber et al., Ann. Statist. 51(2), 2023). Mitigación práctica y suficiente: **recalibrar `q_α` en cada ronda** sobre los datos acumulados (Deshpande, Marx & Kuleshov; Stanton et al., AISTATS 2023), y **reportar la cobertura empírica observada por ronda** en el run log. Llamarlo "cobertura conformal por ronda, no garantía asintótica bajo shift".

**M2 — Ampliar el ensemble a K = 8–10 miembros (prioridad ALTA, coste ~cero).**
La regla M=5 de Lakshminarayanan et al. (2017) es de clasificación; en **regresión** se documenta que M=5 *"does not guarantee sufficient UQ quality"*, con saturación en 8–12 ([arXiv:2303.16210](https://arxiv.org/pdf/2303.16210); Eng. Appl. AI 2024). Con MLPs (96,96,64) en NumPy son segundos. Medir el efecto sobre la cobertura conformal antes/después y dejar el número que la evidencia soporte — **no aumentarlo por fe**.

**M3 — Recalibración isotónica como red de seguridad (prioridad MEDIA, ~30 líneas PAVA).**
Kuleshov, Fenner & Ermon (ICML 2018): mapeo monótono entre cuantiles esperados y observados. Complementa a M1 (conformal da cobertura marginal a un nivel; isotónica da **cobertura correcta a todos los niveles**, útil para el diagrama de fiabilidad del informe de run). Es la métrica que hace auditables las barras de error de Phy-AT.

**M4 — Calibración L2 afín → afín con dependencia física de bajo orden (prioridad MEDIA).**
Sustituir `y_hi ≈ a·y_L0 + b` por `y_hi ≈ (a₀ + a₁·z₁ + a₂·z₂)·y_L0 + b`, con z₁, z₂ = 2 features físicas normalizadas (sugeridas: **Mach absoluto de salida del estator** y **coeficiente de carga de etapa**), con regularización ridge fuerte y validación *leave-one-out* sobre los pares hi-fi. Es el paso mínimo hacia ρ(x) de Kennedy-O'Hagan sin pagar co-kriging completo, y tiene precedente directo en turbomáquinas (*gain factors* de JMSE 14(9):825). **Regla de puerta: solo activar el término lineal si LOO mejora; si no, colapsar a la afín actual.**

**M5 — Features físicas de turbina para el embedding (prioridad ALTA — es *el* trabajo de dominio de Phy-AT).**
El embedding de Phy-AC (12 features) está diseñado para compresor: difusión (DF), margen de surge de Koch, Mach relativo del primer rotor. **Ninguna de esas tres tiene sentido en turbina.** Conjunto de reemplazo sugerido por la literatura (AIAA 2024-1765; Kacker-Okapuu; Benner et al.; Zweifel 1945; Smith 1965; Ainley-Mathieson 1951):

| # | Feature | Motivación / fuente |
|---|---|---|
| 1–2 | **Zweifel del estator y del rotor**, Z = 2·(s/c)·cos²α₂·(tanα₁+tanα₂) | Criterio canónico de carga de álabe de turbina (Zweifel 1945); fija la solidez óptima |
| 3 | **Mach absoluto máximo de salida de estator** (M₂) | Discrimina subsónico / transónico / gargantado; la no linealidad dominante en la pérdida de perfil |
| 4 | **Mach relativo máximo de salida de rotor** (M₃,rel) | Ídem para el rotor; controla la pérdida de choque de borde de fuga |
| 5 | **Margen de garganta / relación de área crítica** (ṁ/ṁ_choke por fila, mínimo) | Reemplaza el "min Koch SM". En turbina el choque es modo normal → **feature continua**, no penalización binaria (cf. J. Turbomach. 147(4):041002) |
| 6 | **Reacción mínima en el cubo** (min R_hub) | La reacción negativa en cubo es el fallo clásico de turbinas de bajo HTR |
| 7 | **Ángulo de giro máximo (deflexión) por fila** | Variable primaria de Ainley-Mathieson y Soderberg |
| 8 | **Reynolds de cuerda mínimo** (Re_c) | Correlacionado con pérdida de perfil en LPT (SPLEEN opera a Re 70k–120k); crítico si Phy-AT cubre LPT |
| 9 | **τ/h — holgura de punta sobre altura de álabe** (última etapa) | Impulsor del tip leakage; confirmado como variable de primer orden por el spanwise loss model ML |
| 10 | **AN²** (área de anillo × RPM²) | Restricción estructural canónica: esfuerzo centrífugo de raíz ∝ AN² |
| 11 | **t_TE/o — bloqueo de borde de fuga sobre garganta** | Término explícito en Kacker-Okapuu; controla la pérdida de borde de fuga |
| 12 | **Swirl de salida de la última etapa** (α_exit) | Penaliza energía cinética residual; determina el difusor de salida |
| 13–14 | **Coordenadas del diagrama de Smith (φ, ψ)** — media y pendiente | Smith (1965) es *literalmente* la superficie de η(φ,ψ) que el surrogate debe aprender |
| 15 | **log(TR) y η_L0** del prior | Análogo directo de log PR_L0 y η_L0 en Phy-AC |

*(Si Phy-AT ha de cubrir turbinas refrigeradas, añadir la **fracción de caudal de refrigeración** y la **temperatura metálica normalizada** — pero solo si L0 las modela; una feature que el prior no consume es ruido.)*

**M6 — Adquisición: ΔHipervolumen esperado por Monte Carlo, en lugar de LCB por objetivo (prioridad MEDIA, riesgo BAJO).**
La contribución de hipervolumen es la métrica correcta de "cuánto mejora el frente" (Daulton et al., qNEHVI, NeurIPS 2021), y en 2 objetivos es O(n log n) en NumPy puro. Implementarlo **junto al LCB actual, con un flag y un A/B sobre los benchmarks de regresión** — no como reemplazo ciego. Evidencia de que hay que verificarlo: *"in some expensive design problems, qNEHVI does not outperform NSGA-II"* (MORBO, [arXiv:2109.10964](https://arxiv.org/pdf/2109.10964)).

**M7 — Región de confianza tipo TuRBO alrededor del mejor factible (prioridad MEDIA-BAJA, ~40 líneas).**
Caja que se encoge tras rondas sin mejora y se expande tras éxitos (Eriksson et al., NeurIPS 2019; con restricciones, SCBO, AISTATS 2021). Mitiga la deriva de la exploración de σ alta hacia rincones del hipercubo en 15-D. **Adoptar solo si el A/B lo justifica**; añade dos hiperparámetros y es la mejora con más superficie de sintonía.

**M8 — Adquisición consciente de la frontera de factibilidad (prioridad MEDIA).**
Con ~20% de región factible, sesgar deliberadamente parte del lote hacia la **frontera g≈0**, no hacia el interior factible ni hacia el infactible profundo. Justificación: la crítica de 2025 a cEI, que *"restricts exploration to the feasible region"* cuando los óptimos restringidos viven en la frontera ([JCGS, doi:10.1080/10618600.2025.2611111](https://www.tandfonline.com/doi/abs/10.1080/10618600.2025.2611111)). Barato aquí porque g es exacta.

**M9 — Anclar la validación a SPLEEN C1 (prioridad ALTA para credibilidad, coste = trabajo de datos).**
Es el mejor test case público experimental de cascada LPT de alta velocidad ([Zenodo 10.5281/zenodo.7264761](https://zenodo.org/records/10253213)). Usar sus condiciones (M₂ 0.7–0.95, Re 70k–120k) como **anclas hi-fi de la calibración L2** y como caso de regresión del prior. Da a Phy-AT algo que Phy-AC no tiene: **validación contra medida experimental abierta y citable**.

**M10 — Modelo de pérdidas de turbina "híbrido explícito" en L0 (prioridad ALTA, es capa 1 pero condiciona 2–4).**
Implementar Kacker-Okapuu con la corrección secundaria de Benner et al. como base, y dejar **enganches (hooks) para factores de corrección aprendidos** por bloque de pérdida (perfil / secundaria / TE / tip leakage), no un único factor global. Justificación: Senior & Miller (J. Turbomach. 146(4):041007, 2024); el modelo híbrido de LPT de Energy 2025 lo confirma. **Referencia de implementación:** TurboFlow (MIT, JOSS 2025) para el tratamiento equation-oriented del choking.

## C. Qué NO adoptar (todavía) — y la razón

| No adoptar | Razón |
|---|---|
| **PINNs en cualquier capa 2–4** | Coste de entrenamiento ≥ CFD por caso; sesgo espectral y concentración de error en gradientes fuertes = el choque de una HPT transónica; requieren reentrenamiento por condición de contorno; no hay evidencia de PINN batiendo a un meanline calibrado en coste-precisión. El residuo verdad−L0 **no obedece ninguna PDE**. *Sí* como herramienta offline de reconstrucción de campos experimentales para alimentar L2 |
| **FNO / DeepONet / neural operators** | Requieren O(10³–10⁴) simulaciones; presupuesto de Phy-AT: 150–500 evaluaciones **totales**. PLAID documenta que **FNO se degrada en mallas anisótropas no estructuradas — Rotor37 incluido** |
| **GNNs / MeshGraphNets** | Mismo problema de datos, dependencia de PyTorch/GPU (rompe el invariante NumPy puro), y *"limitations in providing built-in predictive uncertainties"* (MMGP, NeurIPS 2023) — la UQ es el corazón de la puerta de calidad |
| **Sustituir el ensemble por un GP / co-kriging completo** | Con 15-D, 500 puntos, superficie no estacionaria (choque, n_stages entero) y NumPy sin autodiff, ganancia incierta y coste en fragilidad numérica cierto. La ventaja del GP viene de la **reducción de dimensión**, que el embedding físico ya provee |
| **MFBO / selección de fidelidad en la adquisición (BOCA, MF-MES)** | Requiere ≥3 niveles con costes separados 1–2 órdenes **y decenas de evaluaciones del nivel caro dentro del bucle**. Con <20 puntos hi-fi, estimar ρ e hiperparámetros acoplados tiene más varianza que beneficio. Revisar cuando el flywheel dé ≥50 pares hi-fi; empezar por [arXiv:2410.00544](https://arxiv.org/abs/2410.00544) |
| **qNEHVI/BoTorch, TuRBO como frameworks externos** | Dependencia de PyTorch/GPyTorch. Las **ideas** (hipervolumen esperado, región de confianza) sí — reimplementadas mínimamente en NumPy (M6, M7). Los frameworks no |
| **Modelos generativos (difusión, GAN, VAE) para geometría** | Necesitan miles de geometrías etiquetadas; el problema de capas 2–4 es de **15 variables meanline**, no de forma libre. *"GANs often collapse into mode failure while VAEs struggle…"*. Ver [arXiv:2511.00070](https://arxiv.org/pdf/2511.00070) |
| **Deep RL para diseño** | 10⁴–10⁶ interacciones frente a un presupuesto de 500 |
| **DPP explícito para el lote** | k-means++ **ya es** la aproximación greedy barata a un DPP (BADGE); el DPP exacto cuesta más con mejora marginal |
| **Modelar las restricciones g con un surrogate** | Serían peores que las g exactas del prior — ventaja estructural de la arquitectura |

## D. Riesgos específicos de turbina que Phy-AT hereda mal de Phy-AC

1. **El choque no es un fallo.** En compresor, el choque es una frontera de operación; en turbina, la garganta crítica es normal. La penalización continua de choque de Phy-AC debe reinterpretarse como **feature continua de margen de garganta** + restricción solo cuando la garganta gargantada **invalida el punto de operación pedido**. Referencia de formulación continua: J. Turbomach. 147(4):041002 (2025).
2. **Reacción negativa en cubo.** El fallo característico de turbinas de bajo HTR con alta carga. Debe ser una **g exacta y continua**, no un chequeo posterior.
3. **La escalera Smith(φ,ψ,R) de turbina es distinta de la de compresor.** El rango útil de ψ en turbina llega a 2–2.5 (vs 0.22–0.45 en compresor de Phy-AC) y φ a 0.4–1.0. **Copiar los rangos de DESIGN_VARS de Phy-AC produciría un espacio de diseño físicamente absurdo.** Es lo primero que hay que reescribir.
4. **Bajo Reynolds y transición.** Si Phy-AT cubre LPT, la pérdida de perfil está dominada por burbuja de separación y transición inducida por estelas — no por difusión. El modelo híbrido de Energy 2025 y el caso SPLEEN son las referencias correctas; ignorar Re en el embedding sería un error de primer orden.

---

## Anexo — Referencias clave por prioridad de lectura

**Leer primero (deciden decisiones de arquitectura):**
1. Senior & Miller (2024), *A Data-Centric Approach to Loss Mechanisms*, J. Turbomach. 146(4):041007 — la filosofía híbrida física+ML, del Whittle Lab.
2. Gopakumar et al. (2025), *Uncertainty Quantification of Surrogate Models using Conformal Prediction*, Mach. Learn.: Sci. Technol., [doi:10.1088/2632-2153/ae2e7b](https://iopscience.iop.org/article/10.1088/2632-2153/ae2e7b) — la mejora M1.
3. Casenave et al. (2023), *MMGP*, NeurIPS, [arXiv:2305.12871](https://arxiv.org/pdf/2305.12871) — por qué los métodos gaussianos/reducidos baten a los operadores neuronales en datos escasos.
4. Zou et al. (2024), *Application of AI in turbomachinery aerodynamics*, AI Review 57:222 — mapa del campo.
5. Anderson et al. (2025), *TurboFlow*, JOSS 10(111):7588 + J. Turbomach. 147(4):041002 — meanline de turbina equation-oriented con choking, MIT.

**Leer para las mejoras:**
6. Kuleshov, Fenner & Ermon (2018), ICML — recalibración isotónica (M3).
7. Fannjiang et al. (2022), PNAS + Barber et al. (2023), Ann. Statist. 51(2) — los límites honestos de conformal bajo aprendizaje activo.
8. Daulton, Balandat & Bakshy (2021), qNEHVI, NeurIPS + Daulton et al. (2022), MORBO, [arXiv:2109.10964](https://arxiv.org/pdf/2109.10964) — hipervolumen y su matiz frente a NSGA-II (M6).
9. Eriksson et al. (2019) TuRBO + Eriksson & Poloczek (2021) SCBO (M7).
10. Base de datos **SPLEEN C1**, [Zenodo doi:10.5281/zenodo.7264761](https://zenodo.org/records/10253213) (M9).

**Leer para saber qué NO hacer:**
11. Krishnapriyan et al. (2021) NeurIPS + Wang, Teng & Perdikaris (2021) SIAM JSC — modos de fallo de PINNs.
12. [PLAID Benchmarks](https://plaid-lib.readthedocs.io/en/latest/plaid_benchmarks.html) — FNO se degrada en mallas anisótropas (Rotor37).
13. *Best Practices for Multi-Fidelity Bayesian Optimization*, [arXiv:2410.00544](https://arxiv.org/abs/2410.00544) — cuándo MFBO **no** compensa.
14. *Benchmarking Generative AI Against Bayesian Optimization for Constrained Multi-Objective Inverse Design*, [arXiv:2511.00070](https://arxiv.org/pdf/2511.00070) — el contrapeso al hype generativo.
