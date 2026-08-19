# Informe: XLB, JAX-Fluids y NVIDIA Warp como motor de alta fidelidad para Phy-AT

**Fecha**: 19 ago 2026 · **Método**: clonado y lectura directa del código fuente de los tres repositorios + búsqueda web.
**Marcadores de fiabilidad**: `[V-code]` verificado leyendo el código fuente clonado · `[V-doc]` verificado leyendo README/docs oficiales del repo · `[S]` procede de resúmenes de búsqueda (dominios académicos bloqueados por el proxy en la sesión de investigación; github/raw.githubusercontent accesibles).

---

## 0. Resumen ejecutivo

**Ninguna de las tres librerías puede correr hoy un pasaje de turbina axial transónica refrigerada con y+≈1.** No es cuestión de esfuerzo de integración: les faltan piezas físicas de primer orden, no de conveniencia.

| Librería | Qué es realmente | Bloqueante duro para turbina |
|---|---|---|
| **XLB** | LBM **incompresible isotermo** (Ma≲0.3), rejilla cartesiana + IB | No tiene modelo compresible ni térmico. "Supersonic flows" está en la **wishlist**, no en el roadmap activo `[V-doc]` |
| **JAX-Fluids** | FV compresible de alto orden, **cartesiano puro**, ILES | Sin malla body-fitted, sin marco rotatorio, sin periodicidad de paso, sin NSCBC, sin cp(T), **sin modelo de pared** `[V-code]` |
| **Warp** | **No es un solver**: compilador Python→CUDA con AD | Escribir un RANS de turbomáquina desde cero ≈ 3–5 años-persona; ningún precedente publicado `[V-code, V-doc]` |

**Recomendación de una línea**: mantener el plan actual (**emitir BCs → CFD externo/MULTALL → pares escalares**), y capturar el valor del "ángulo diferenciable" donde sí es barato y sí funciona: **haciendo diferenciables los modelos L0/L1 de Phy-AT (NumPy→JAX), no el CFD**. Vigilar JAX-Fluids como generador de datos 2D/cuasi-3D de física canónica, y **SU2** (que ya tiene mixing-plane + NRBC + adjunto discreto validado en turbinas) como la alternativa realista si algún día se quiere gradiente a través de la alta fidelidad.

---

## 1. XLB — Accelerated Lattice Boltzmann (Autodesk)

**Repo**: https://github.com/Autodesk/XLB · **Licencia**: Apache-2.0 `[V-code]` · Último commit leído: `9470e54` (29 may 2026, "Differentiable LBM example") → proyecto vivo. **Paper**: Ataei & Salehipour, CPC 300:109187 (2024) `[S]`; [arXiv:2311.16080](https://arxiv.org/abs/2311.16080).

### 1.1 Inventario del código `[V-code]`

| Componente | Contenido real |
|---|---|
| Colisiones | BGK, KBC, Smagorinsky-LES-BGK, forced. **Sin cumulante, MRT ni RR/HRR** |
| Equilibrios | cuadrático (Hermite 2º orden) — **solo el isotermo estándar** |
| Velocity sets | D2Q9, D3Q19, D3Q27. **Sin D3Q39 ni retículas de alto orden** (las necesarias para compresible) |
| Steppers | `IncompressibleNavierStokesStepper` (literalmente ese nombre), multires, IBM |
| BCs | equilibrium, bounce-back (full/halfway/interpolado), Zou-He, regularized, extrapolation-outflow, hybrid |
| Precisión | **FP64 solo en el backend JAX**; con Warp/Neon lanza `ValueError` |
| Backends | JAX (multi-GPU shardmap), Warp (single-GPU), Neon (multi-GPU + multiresolución) |

### 1.2 El problema duro: LBM es incompresible e isotermo — y XLB no lo resuelve

El grep de `thermal|compressible|temperature|energy` sobre `xlb/` devuelve **cero** implementaciones de física térmica o compresible `[V-code]`. El README lo confirma `[V-doc]`: "Fluid-Thermal" y "Adjoint-based Shape Optimization" en **Work in Progress**; **"Supersonic Flows"** en la **Wishlist** ("Contributions welcome. Please submit PRs"). Para M 0.7–1.3 y γ variable: descalificatorio de partida.

### 1.3 Lo que la literatura tuvo que construir (y XLB no tiene)

El LBM transónico industrial **existe**, pero es otro método `[S]`: **HRR térmico** (LB para masa+momento + FV acoplado para la energía — Feng, Boivin, Jacob, Sagaut, JCP 2019); **LS89 transónico con LBM** (2024) resuelto con **ProLB** (comercial, CERFACS/Safran), no con XLB; geometrías rotativas vía **overset/chimera** (2021), que XLB no tiene; **PowerFLOW** sí hace transónico (hasta Mach 2) pero con conmutación **D3Q19↔D3Q39** — XLB no tiene D3Q39 `[V-code]`. El camino XLB→turbina no es "configurar": es **reescribir el núcleo físico de la librería**.

### 1.4 La malla: el argumento que mata el caso incluso ignorando la compresibilidad

Para una tobera tipo LS89 (c=50 mm, Re=10⁶): y+=1 → Δy₁ ≈ 1.28 µm.

| Estrategia | Celdas | Veredicto |
|---|---|---|
| Cartesiano uniforme a 1.28 µm | ≈ 10¹⁴ (~10 PB) | Absurdo |
| Multiresolución (Neon), banda fina de 1 mm | ≈ 3×10¹² (~300 TB) | Absurdo |
| Con **modelo de pared** y+≈50 | ≈ 10⁸ | El régimen de ProLB/PowerFLOW |
| **Body-fitted O-H, RANS, y+≈1** | 10⁶–5×10⁶ | Lo que hacen MULTALL/Turbostream |

**El modelo de pared es la condición de existencia del LBM cartesiano** (baja el coste 4 órdenes) — **y XLB no tiene ninguno** `[V-code]` (su canal turbulento es Re_τ=180, DNS de juguete). Además Δt ∝ Δx en LBM: refinar en pared multiplica también los pasos.

### 1.5 Rotor-estátor y periodicidad

`grep -i "coriolis|centrifugal|rotating frame|MRF"` → **cero resultados** `[V-code]`. No hay BC periódica rotacional (solo periodicidad cartesiana del streaming). Lo que sí hay: `wind_turbine_ibm.py` (IBM con marcadores lagrangianos móviles, Ma≈0.01, flujo abierto) y `rotating_sphere_3d.py` — **irrelevantes** para un pasaje confinado transónico.

### 1.6 Diferenciabilidad — el hallazgo más importante

Del docstring de `examples/cfd/differentiable_lbm.py`, **literal** `[V-code]`:

> *"XLB's Warp stepper kernel doesn't have adjoint implementations, so gradients are zero when using wp.Tape (verified by test_stepper_autodiff.py). JAX uses source transformation which works through the stepper."*

1. La diferenciabilidad **solo existe en el backend JAX** (el lento). 2. El backend Warp devuelve **gradientes cero**. 3. Neon es fork de Warp → tampoco. **Rendimiento y diferenciabilidad son mutuamente excluyentes en XLB hoy.** El ejemplo diferenciable es 128×128, 2D, 50 pasos, optimiza condiciones iniciales para dibujar una letra — **no es optimización de forma**. Ningún caso publicado de optimización de forma con XLB `[S]`.

### 1.7 Veredicto XLB

> **DESCARTAR** para L3 de turbina axial: solver incompresible isotermo, compresible en wishlist, sin modelo de pared, sin marco rotatorio ni periodicidad de paso, y con la diferenciabilidad rota en el backend rápido. Único uso residual: visualización cualitativa externa a bajo Ma — no justifica la dependencia. Revisar solo si el CHANGELOG añade stepper térmico/compresible con retícula de alto orden **y** wall model.

---

## 2. JAX-Fluids (Bezgin, Buhendwa, Adams — TUM)

**Repo**: https://github.com/tumaer/JAXFLUIDS · **Licencia**: **MIT** `[V-code]` · Último commit `af2b7cf` (5 ago 2026) → muy activo. **Papers**: CPC 284:108527 (2023) y CPC 308:109433 (2025, "JAX-Fluids 2.0").

### 2.1 Qué hace bien `[V-code/V-doc]`

FV compresible explícito (Euler/RK2/RK3); **WENO-3/5/7, WENO-CU6, TENO-5/6/8, WENO-3NN**; Riemann HLL/HLLC/HLLC-LM/AUSM+/Rusanov/Roe/CATUM; **level-set** sharp-interface para sólidos inmersos y dos fases; **positivity-preserving**; `jax.pmap` hasta **512 A100 / 2048 TPU-v3**; FP64; e infraestructura ML de primera (`jaxfluids_nn`, RL env, `feed_forward` con **`jax.checkpoint` en dos niveles** — la maquinaria exacta del learning-in-the-loop).

### 2.2 Auditoría contra los requisitos de una turbina `[V-code]`

| Requisito | ¿Existe? | Evidencia |
|---|---|---|
| Malla body-fitted / curvilínea / ALE | ❌ Solo cartesiano con estiramiento tensor-producto por eje (PIECEWISE, CHANNEL, BOUNDARY_LAYER) | `domain/mesh_creation/` |
| Periodicidad rotacional de paso | ❌ PERIODIC es cartesiana pura (sin rotación del vector en el halo) | `halos/outer/` |
| Marco rotatorio (Coriolis+centrífuga) | ❌ Solo gravedad y términos axisimétricos | `solvers/source_term_solver.py` |
| Malla deslizante / mixing plane | ❌ (los "moving solids" del ej. 13 no son eso) | |
| BCs no reflectantes (NSCBC/Giles) | ❌ Solo sponge layer como mitigación | `halos/outer/__init__.py` |
| Gas caliente cp(T), γ(T), R(FAR) | ❌ `IdealGas`: γ y cp **constantes**; `get_specific_heat_capacity(T)` ignora T | `materials/single_materials/ideal_gas.py` L29-48 |
| RANS | ❌ | |
| LES explícito / **modelo de pared** / transición | ❌ Solo ILES (ALDM, TENO) | |
| Refrigeración/film | Parcial: `MASSTRANSFERWALL` solo en paredes planas del dominio, no en el level-set | |

### 2.3 El coste

Level-set sharp-interface sin función de pared → **wall-resolved obligatorio en cartesiano** → los mismos 10¹²–10¹⁴ celdas de §1.4. Y es **explícito**: Δ=1.28 µm → Δt≈2.6 ns → ~2×10⁶ pasos para estadísticas. El estiramiento tensor-producto no salva un perfil combado (refinar y para el intradós refina todo el pasaje).

### 2.4 Lo que SÍ podría hacer

Los ejemplos existentes (`cylinder`, `double_mach`, `bowshock`, `NACA`, `diamond_airfoil`, TGV, HIT) marcan la frontera — **ni un caso de cascada ni de turbomáquina** `[V-code]`. Usos acotados: (1) cascada 2D a Re 10⁴–10⁵ ILES sin gas caliente — ejercicio académico, no un par (Γ, η_tt); (2) **estudios canónicos de física de pérdidas** (SBLI, burbuja de separación, mezcla de estela) como priors cualitativos; (3) banco de aprendizaje de esquemas numéricos (para lo que fue diseñado).

### 2.5 Diferenciabilidad y coste del backprop

Genuinamente buena: `lax.scan` + `jax.checkpoint` en dos niveles `[V-code]`; gradientes end-to-end demostrados (ML-ILES; control de flujo por AD — [arXiv:2410.23415](https://arxiv.org/abs/2410.23415)) `[S]`. Memoria: 10⁶ celdas × 2000 pasos sin remat = 80 GB (imposible); con checkpointing √N ≈ 1.8 GB + ~2× tiempo. **El backprop es viable en 2D/3D pequeño — justo la escala donde JAX-Fluids puede correr una cascada.** Coherente, pero pequeño.

### 2.6 Veredicto JAX-Fluids

> **VIGILAR**, con rol experimental acotado — **nunca L3**. El mejor código diferenciable compresible abierto (MIT, alto orden, positivity, checkpointing, HPC probado), pero su restricción cartesiana es **arquitectónica**. Añadir marco rotatorio (~1 sem), periodicidad de paso (~1 mes), NSCBC (~2 meses), cp(T) (~1 sem) y wall model (~6 meses) ≈ **1 año-persona y seguiría sin poder mallar un álabe combado a y+≈1**. Rol: *banco de física canónica 2D* para priors cualitativos y experimentos de aprendizaje de cierres. **Ninguna salida suya entra como par hi-fi en L2.**

---

## 3. NVIDIA Warp

**Repo**: https://github.com/NVIDIA/warp · **Licencia**: **Apache-2.0** `[V-code]` · v1.17.0.dev4, commit del 19 ago 2026.

### 3.1 Qué es realmente

**No es un solver CFD**: transpilador Python→CUDA/CPU con AD por transformación de fuente (`wp.Tape`), BVH/HashGrid/NanoVDB, `warp.fem` (grid/quad/tri/tet/hex/nanogrid — sí mallas no estructuradas), `warp.sparse`, `warp.optim`, interop JAX/PyTorch. **`warp.sim` ya no existe en el árbol** (externalizado a Newton) `[V-code]` — la API de alto nivel se mueve.

### 3.2 Limitaciones de AD que importan (de `docs/user_guide/differentiability.rst`, leído) `[V-code]`

- **`*=` y `/=` in-place NO son diferenciables** ("incorrect results in the backward pass"); solo `+=`/`-=`.
- Componentes de vector/matriz asignables con `=` **una sola vez**.
- **Los bucles dinámicos no se reproducen en el backward pass** (~250 líneas de workarounds en el manual).

Un residuo RANS real está lleno de bucles sobre caras, limitadores min/max y updates multiplicativos: **el adjunto es un campo de minas**.

### 3.3 ¿CFD serio en Warp?

`PUBLICATIONS.md` (~130 entradas 2024–2026) `[V-code]`: ~90 % robótica/gráficos/MPM. CFD compresible: **una** entrada (Chamarthi, [arXiv:2604.02757](https://arxiv.org/abs/2604.02757)). **Ningún solver de turbomáquinas, ninguno RANS body-fitted multibloque.** Los ejemplos de fluidos del repo son demos (Stam 2D, SPH, NS 2D incompresible FEM).

### 3.4 ¿Qué costaría el RANS de turbomáquina en Warp?

La lista mínima para igualar a MULTALL: malla multibloque O-H, flujos+limitadores, SA/SST con pared, marco rotatorio con rotalpía, periodicidad + **mixing plane conservativo** + Giles, gas imperfecto con búsqueda de choking, multigrid/paso local, fugas/cavidades/refrigeración, y V&V completa. Es **exactamente lo que Turbostream 3 hizo a mano en CUDA** (Brandvik & Pullan 2011) con un equipo del Whittle Lab y ~15 años de continuidad. Warp ahorraría ~20 % (gestión de kernels), **nada** del 80 % que es física numérica. Y el "adjunto gratis" tiene letra pequeña: el adjunto de un RANS estacionario no se obtiene diferenciando el pseudo-tiempo — se resuelve el sistema adjunto sobre el estado convergido; Warp solo da la derivada del kernel.

### 3.5 Lock-in y madurez

CUDA-only en GPU (hay CPU, sin ROCm/Metal) — y **el usuario de Phy-AC corre en laptop CPU**. Churn de API real (warp.sim fuera, deprecaciones) — un pasivo para la reproducibilidad bit-exacta que Phy-AC se autoimpone.

### 3.6 Veredicto Warp

> **DESCARTAR como base de solver** (3–5 años-persona para replicar lo que MULTALL da gratis y Turbostream vende). **VIGILAR como utilidad táctica** en dos nichos: aceleración GPU de geometría (BVH/SDF/voxelización) para 5a/5c, y FEM estructural diferenciable — ambos problemas que Phy-AC ya resuelve en NumPy en milisegundos, así que hoy no hay problema que resolver.

---

## 4. Ecosistema cercano — una línea cada uno

| Proyecto | Método | GPU | ¿Turbomáquinas? | Licencia | Nota |
|---|---|---|---|---|---|
| **[MULTALL](https://github.com/paopaoai11/Multall-open-18.3)** | RANS 3D multietapa + Q3D + throughflow | ❌ CPU | ✅ | dominio público (F77) | **La referencia: minutos en CPU** |
| **Turbostream 3/4** | RANS/URANS GPU | ✅ | ✅ | comercial | 3 etapas con fugas <10 min/4 GPU |
| **[turbigen](https://gitlab.developers.cam.ac.uk/jb753/turbigen)** | diseño+exploración | vía TS | ✅ | GPL-3 | el solver que llama es comercial |
| **[SU2](https://su2code.github.io/)** | RANS/URANS FV no estructurado | parcial | ✅ **mixing plane + NRBC Giles + adjunto discreto (CoDiPack) validado en turbinas** | LGPL-2.1 | ⭐ **el competidor serio si se quiere gradiente hi-fi** |
| **[PyFR](https://pyfr.org/)** | FR alto orden | ✅ (2048 GH200) | ✅ ILES/DNS T106C, MTU-T161 | BSD-3 | LES de cascada; sin AD ni mixing plane |
| **[PhysicsNeMo](https://github.com/NVIDIA/physicsnemo)** (ex-Modulus) | PINNs/operadores | ✅ | solo aero externa | Apache-2.0 | surrogate, no solver — mismo veredicto que el informe 04 |
| **[JAX-CFD](https://github.com/google/jax-cfd)** | incompresible diferenciable | ✅ | ❌ | Apache-2.0 | el paper fundacional (Kochkov PNAS 2021), inaplicable |
| **[PhiFlow](https://github.com/tum-pbs/PhiFlow)** | PDE diferenciable incompresible | ✅ | ❌ | MIT | pedagógico |
| **[Lettuce](https://github.com/lettucecfd/lettuce)** | LBM+PyTorch dif. | ✅ | ❌ | MIT | mismos límites que XLB |
| **[FluidX3D](https://github.com/ProjectPhysX/FluidX3D)** | LBM OpenCL (8799 MLUPS/A100) | ✅ | ❌ | **no comercial** | ⛔ la licencia lo descalifica |
| waLBerla / OpenLB / Palabos | LBM HPC | ✅ | ❌ | GPL/varias | sin transónico |
| **[MFC](https://mflowcode.github.io/)** | FV multifase compresible | ✅ (43k GPUs) | ❌ | MIT | dominio equivocado |
| **[STREAmS-2](https://github.com/STREAmS-CFD/STREAmS-2)** | DNS compresible canónico | ✅ | ❌ | GPL-3 | canal/BL/SBLI |
| **[Trixi.jl](https://github.com/trixi-framework/Trixi.jl)** | DG adaptativo Julia | parcial | ❌ | MIT | AD nativa, inmaduro en GPU |
| **[Nektar++](https://www.nektar.info/)** | DG/CG espectral | parcial | ✅ LES cascada | MIT | alternativa a PyFR |
| **JAX-FVM** ([arXiv:2607.07385](https://arxiv.org/html/2607.07385), 2026) | **FV entropy-stable NO estructurado, JAX, diferenciable** | ✅ | ❌ aún | ? | ⭐⭐ **rompe la barrera "diferenciable ⇒ cartesiano"** — lo más prometedor del horizonte |
| **DiFVM** ([arXiv:2603.15920](https://arxiv.org/pdf/2603.15920), 2026) | FV dif. poliédrico, compatible OpenFOAM | ✅ | ❌ | ? | ⭐ la otra vía; preprint sin código maduro |

**Respuesta directa a "¿qué aportaría reescribir en JAX/Warp?": exactamente una cosa — el gradiente. Nada más.** En física, cobertura, validación y tiempo-a-resultado, MULTALL en CPU gana a cualquier cosa escrita en JAX/Warp en los próximos 5 años. Y si se quiere gradiente sobre alta fidelidad, **SU2 ya lo tiene validado en turbinas** — cuesta días, no años.

---

## 5. El ángulo diferenciable: por qué SÍ interesa, y por qué hoy no pasa por el CFD

### 5.1 La promesa (real y documentada)

Kochkov et al. (PNAS 2021, 8–10× resolución efectiva aprendiendo correcciones dentro del solver); ML-ILES de JAX-Fluids; *Differentiable Turbulence* ([arXiv:2307.03683](https://arxiv.org/pdf/2307.03683)) `[S]`. Traducido: calibrar los **coeficientes de los modelos de pérdidas** minimizando ‖y_L1(θ_loss) − y_hifi‖² por gradiente exacto.

### 5.2 El malentendido a deshacer

**Para calibrar los modelos baratos por gradiente NO hace falta que el solver caro sea diferenciable.** El gradiente necesario es ∂y_L1/∂θ_loss — la derivada del **modelo barato** respecto a **sus propios coeficientes**; y_hifi es una constante del problema. Diferenciar el CFD solo hace falta para optimización de **forma** con adjunto o para cierres **incrustados en el propio CFD** — ninguno es el cuello de botella de Phy-AT.

### 5.3 Comparación con el lazo actual (leída la `HiFiCalibration` real de `physics_core.py` `[V-code]`)

| Aspecto | Afín actual (a·y+b) | Gradiente sobre L0/L1 diferenciable | Adjunto a través del CFD |
|---|---|---|---|
| Parámetros | 2 por salida, global | 10–30 coeficientes físicos | ∞ (forma) |
| Pares necesarios | ≥2 | 15–40 | 1/iteración |
| Extrapolación | ❌ recta global | ✅ la física del modelo la lleva | n/a |
| Implementación | hecha | **2–4 semanas** (meanline en jax.numpy) | años, o SU2 |
| Hardware | laptop CPU ✅ | **laptop CPU ✅** | GPU ❌ |
| Sobreajuste | bajo | medio (regularizar hacia valores publicados) | alto |

**El salto de calidad disponible — de la recta de 2 parámetros a la calibración física regularizada de 10–30 coeficientes — no requiere GPU, ni XLB, ni JAX-Fluids, ni Warp**: requiere escribir los modelos de Phy-AT con `jax.numpy` (casi una sustitución de import, con `.at[].set()` en vez de asignación in-place). Encaja con el KOH del informe 04: y_hifi = ρ·y_L1(θ, θ_loss) + δ(θ).

### 5.4 Coste realista

El usuario de Phy-AC corre en **laptop CPU** — XLB/Warp sin GPU son inútiles; JAX-Fluids en CPU hace 2D pequeño. Una A100 por horas es barata, pero 10¹² celdas no caben en ninguna cantidad razonable de A100s. Ninguna de las tres pasa el listón de "cero recalibración por máquina".

---

## 6. Recomendación para Phy-AT

### 6.1 Tabla comparativa

| | **XLB** | **JAX-Fluids** | **Warp** | *(ref)* MULTALL | *(ref)* SU2 | *(ref)* PyFR |
|---|---|---|---|---|---|---|
| Método | LBM incompresible | FV compresible WENO/TENO ILES | framework kernels | RANS 3D estructurado | RANS FV no estructurado | FR alto orden |
| Compresible transónico | ❌ (wishlist) | ✅ | n/a | ✅ | ✅ | ✅ |
| Body-fitted | ❌ (IB cartesiano) | ❌ (cartesiano+level-set) | fem ✅ | ✅ O-H | ✅ | ✅ |
| Rotor-estátor | ❌ | ❌ | n/a | ✅ mixing plane | ✅ mixing plane + Giles | ⚠️ |
| Periodicidad de paso | ❌ | ❌ | n/a | ✅ | ✅ | ✅ |
| cp(T) | ❌ | ❌ (constantes) | n/a | ✅ | ✅ | ⚠️ |
| Turbulencia/pared/transición | Smagorinsky/❌/❌ | ILES/❌/❌ | n/a | ML+pared | SA, SST, γ-Reθ | ILES |
| AD/adjunto | solo backend JAX (Warp=grad 0) | ✅ end-to-end + checkpointing | ✅ con 3 limitaciones serias | ❌ | ✅ **adjunto discreto** | ❌ |
| Multi-GPU / FP64 | ✅/solo JAX | ✅ (512 A100)/✅ | ✅/✅ | ❌/✅ | ⚠️/✅ | ✅/✅ |
| Licencia | Apache-2.0 | MIT | Apache-2.0 | abierta | LGPL-2.1 | BSD-3 |
| **Madurez para turbina** | **0/10** | **1/10** | **0/10** | **8/10** | **8/10** | **6/10** |

### 6.2 Veredictos

- 🔴 **XLB — DESCARTAR** (rol: ninguno). Ecuaciones equivocadas, malla equivocada, sin BCs de turbina, diferenciabilidad rota en el backend rápido. Re-evaluar solo con stepper compresible + wall model en el CHANGELOG.
- 🟡 **JAX-Fluids — VIGILAR** (rol: banco de física canónica 2D + laboratorio de cierres; **nunca L3, nunca pares de calibración**). Vigilar también JAX-FVM y DiFVM (AD sobre malla no estructurada — el verdadero horizonte).
- 🔴/🟡 **Warp — DESCARTAR como base de solver / VIGILAR como utilidad** (geometría GPU en 5a/5c; FEM estructural diferenciable — hoy sin problema que resolver).
- 🟢 **Adoptar: JAX aplicado a los modelos BARATOS.** Portar el meanline de Phy-AT (y opcionalmente el SCM) a `jax.numpy`, exponer los coeficientes de pérdidas como parámetros y calibrarlos por gradiente con regularización hacia los valores publicados. Coste: 2–4 semanas. Hardware: laptop CPU. Beneficio: sustituye la recta de 2 parámetros por calibración física estructurada que extrapola.

### 6.3 Contra el plan actual

El plan **"emitir BCs → CFD externo → pares escalares"** es correcto y se mantiene: coste marginal cero, es la práctica industrial, el contrato escalar es robusto/auditable, y corre en laptop. Sus debilidades se atacan **sin** las tres librerías:

| Debilidad | Remedio |
|---|---|
| La afín de 2 parámetros no extrapola | **JAX sobre L0/L1** (§6.2) |
| Depender del usuario para cada par es lento | **Integrar MULTALL como L3 interno de referencia** (Meangen/Stagen desde el contrato → parser de salida → decenas de pares desatendidos en CPU). **El mayor retorno por esfuerzo de este informe** |
| Los pares escalares no informan del reparto spanwise | ampliar el contrato con **perfiles radiales** (α_exit(r), p₀(r)) que MULTALL/CFX ya producen |
| No hay gradiente de forma | si algún día hace falta: **SU2**, no un solver propio |

### 6.4 Plan de acción (prioridad decreciente)

1. **[Ya]** Meanline L0 de Phy-AT en `jax.numpy` **desde el primer día** (coste ~0 en diseño inicial; carísimo como retrofit).
2. **[Pronto]** Empaquetar **MULTALL** como L3 interno opcional.
3. **[Mantener]** El paquete de BCs para CFD externo como fidelidad primaria.
4. **[Vigilar, sin invertir]** JAX-Fluids; JAX-FVM y DiFVM; SU2+GPU.
5. **[No hacer]** No integrar XLB. No escribir un solver en Warp. No perseguir PINNs/operadores como sustituto de física.

---

## 7. Marketing vs realidad del código

| El marketing dice | La realidad del código |
|---|---|
| XLB: "fully differentiable… state-of-the-art performance" | Diferenciable **o** rápido, nunca a la vez (backend Warp: gradientes **cero**, verificado en su propio test) |
| XLB: "Fluid-Thermal Simulation Capabilities" | En **Work in Progress**; "Supersonic Flows" en la **Wishlist** |
| XLB: el GIF del aerogenerador | IBM a Ma≈0.01 en flujo abierto — cero relación con un pasaje transónico |
| JAX-Fluids: "fully-differentiable CFD solver for 3D compressible flows" | Cierto, **y** cartesiano tensor-producto, sin RANS, sin wall model, sin cp(T), sin BCs características |
| Warp: "differentiable physics simulation" | Para CFD: `*=` rompe el gradiente, componentes de vector una sola asignación, bucles dinámicos no se reproducen en backward |
| Warp: ecosistema | 130 publicaciones, ~90 % robótica/gráficos, **cero turbomáquinas** |
| PhysicsNeMo: "physics-ML for CFD" | Surrogates que necesitan datos que Phy-AT no tiene; no es un solver |

---

## Fuentes

**Verificadas por lectura de código clonado**: [Autodesk/XLB](https://github.com/Autodesk/XLB) (commit 9470e54) · [tumaer/JAXFLUIDS](https://github.com/tumaer/JAXFLUIDS) (commit af2b7cf) · [NVIDIA/warp](https://github.com/NVIDIA/warp) (v1.17.0.dev4) · `physics_core.py` de Phy-AC (`HiFiCalibration`).

**De resúmenes de búsqueda**: XLB CPC 300:109187 / [arXiv:2311.16080](https://arxiv.org/abs/2311.16080) · JAX-Fluids [arXiv:2203.13760](https://arxiv.org/abs/2203.13760) / [2402.05193](https://arxiv.org/abs/2402.05193) · HRR térmico (Feng et al. JCP 2019) · LBM del LS89 transónico (ProLB, 2024) · HRR+overset rotativo (2021) · etapa de compresor HP con LBM (J. Turbomach. 146(11):111007) · PowerFLOW transónico (blogs SIMULIA) · SU2 adjunto en turbomáquinas ([J. Prop. Power](https://arc.aiaa.org/doi/10.2514/1.B37685), [arXiv:2405.06056](https://arxiv.org/pdf/2405.06056)) · MULTALL · turbigen · PyFR (T106C, MTU-T161) · Kochkov PNAS 2021 · [arXiv:2307.03683](https://arxiv.org/pdf/2307.03683) · [arXiv:2410.23415](https://arxiv.org/pdf/2410.23415) · JAX-FVM [arXiv:2607.07385](https://arxiv.org/html/2607.07385) · DiFVM [arXiv:2603.15920](https://arxiv.org/pdf/2603.15920) · JAX-CFD · PhiFlow · Lettuce · FluidX3D · MFC · STREAmS-2 · Trixi.jl · Nektar++ · PhysicsNeMo (DoMINO, Transolver).
