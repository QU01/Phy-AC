# Métodos Through-Flow / Streamline Curvature y Equilibrio Radial para TURBINAS AXIALES
## Investigación de base para la fidelidad L1 de Phy-AT (Quasar)

---

### Nota metodológica previa (honestidad sobre las fuentes)

La sesión de investigación corre tras un proxy de egreso con lista blanca. **Bloqueados** (403): `asmedigitalcollection.asme.org`, `ntrs.nasa.gov`, `link.springer.com`, `sciencedirect.com`, `arxiv.org`, `researchgate.net`, `mdpi.com`, `euroturbo.eu`, `pcaeng.co.uk`, `reports.aerade.cranfield.ac.uk`, `joss.theoj.org`, `ora.ox.ac.uk`, `dtic.mil`, `archive.org`, `iopscience`, `api.github.com`. **Accesibles**: motor de búsqueda (200 consultas, agotadas) y `raw.githubusercontent.com`.

Marcadores de fiabilidad:
- **[V-CÓDIGO]** — verificado leyendo código fuente abierto descargado en la sesión (TurboFlow, DTU/NTNU). Fórmulas exactas, transcritas del fuente.
- **[V-ABS]** — verificado contra el resumen/abstract indexado del artículo.
- **[REF]** — referencia bibliográfica completa verificada desde el `bibliography.bib` de TurboFlow; el *contenido* es conocimiento de dominio estándar.
- **[DOM]** — conocimiento de dominio no re-verificado contra fuente primaria en la sesión. "Hipótesis fuerte a confirmar" antes de codificar.

---

## 1. Fundamentos SCM para turbinas

### 1.1 Genealogía del método

| Año | Autor | Aportación | Estado |
|---|---|---|---|
| 1952 | Wu, C.-H., *A General Theory of Three-Dimensional Flow…*, NACA TN 2604 | Superficies **S1** (álabe-a-álabe) y **S2** (meridional); todo through-flow moderno es una S2 promediada circunferencialmente | [V-ABS] — [NTRS 19930083325](https://ntrs.nasa.gov/citations/19930083325) |
| 1966 | Smith, L. H. Jr., *The Radial-Equilibrium Equation of Turbomachinery*, J. Eng. Power 88(1):1–12 | La ecuación DENTRO de la fila, con interpretación física de cada término | [V-ABS] |
| 1967 | Novak, R. A., *Streamline Curvature Computing Procedures for Fluid-Flow Problems*, J. Eng. Power 89(2):478–490, DOI 10.1115/1.3616716 | La formulación numérica que implementan todos los SCM industriales: cuasi-ortogonales, curvatura retrasada, continuidad por integración radial | [V-ABS] |
| 1969/70 | Wilkinson, D. H., *Stability, Convergence and Accuracy of 2-D SCM Using Quasi-Orthogonals*, IMechE Paper 35 | Factor de amortiguamiento óptimo en función de la relación de aspecto de malla y el Mach; diferencias polinómicas mejores que splines para la 2ª derivada | [V-ABS] |
| 1970 | Frost, D. H., ARC R&M 3687 | Primer SCM público con **estaciones internas dentro de la fila** | [V-ABS] |
| **1978** | **Denton, J. D., *Throughflow Calculations for Transonic Axial Flow Turbines*, J. Eng. Power 100(2):212–218, DOI 10.1115/1.3446336** | **La referencia central para turbinas.** "El choking ocurre en la garganta entre dos álabes; una buena estimación del gasto de bloqueo debe usar estimaciones precisas de las áreas de garganta y de su posición" | [V-ABS] |
| 1993 | Denton, *Loss Mechanisms in Turbomachines*, J. Turbomach. 115(4):621–656 | Marco entrópico de pérdidas; modelo simple de fuga de punta | [REF] |
| 2010 | Casey, M. & Robinson, C., *A New Streamline Curvature Throughflow Method for Radial Turbomachinery*, J. Turbomach. 132(3):031021 | Paredes/LE/TE curvos, estaciones internas, mezcla spanwise, **redistribución del flujo en el span por choking**, flujo relativo supersónico | [V-ABS] |
| 2013 | Tiwari, Stein & Lin, *Dual-Solution and Choked Flow Treatment in a SCM Throughflow Solver*, J. Turbomach. 135(4):041004 | **La doble solución sub/supersónica** de la continuidad y cómo guiar al solver a la rama correcta | [V-ABS] |
| 2013 | Petrović & Wiedermann, *Through-Flow Analysis of Air-Cooled Gas Turbines*, J. Turbomach. 135(6):061019 | Through-flow (función de corriente + FEM) para turbinas multietapa **refrigeradas**: distribución radial de pérdidas, mezcla spanwise, film cooling, inyección de TE, flujos de disco | [V-ABS] |

### 1.2 La ecuación de equilibrio radial en turbinas

La ecuación de Phy-AC se traslada literalmente al marco absoluto (estátor):

$$C_m\frac{\partial C_m}{\partial r} = \frac{\partial h_0}{\partial r} - T\frac{\partial s}{\partial r} - \frac{C_u}{r}\frac{\partial (r C_u)}{\partial r} - \frac{C_m^{2}\cos\gamma}{r_c}$$

Para el **rotor** hay que sustituir h₀ por la **rotalpía** (conservada en el marco relativo):

$$I = h_0^{\text{rel}} - \tfrac12 U^2 = h + \tfrac12 W^2 - \tfrac12 (\Omega r)^2 \quad\Rightarrow\quad \frac{DI}{Dt}=0 \ \text{(adiabático, sin refrigeración)}$$

$$C_m\frac{\partial C_m}{\partial r} = \frac{\partial I}{\partial r} - T\frac{\partial s}{\partial r} - \frac{W_u}{r}\frac{\partial (r C_u)}{\partial r} - \frac{C_m^{2}\cos\gamma}{r_c} + F_r \qquad\text{[DOM, forma estándar]}$$

**El término de fuerza del álabe F_r es lo que Phy-AC omite y Phy-AT no puede omitir** (dentro de fila, con lean): la fuerza promediada es normal a la superficie del álabe:

$$\frac{F_r}{F_\theta} = -\tan\lambda \qquad \rho\,C_m\,\frac{\partial (r C_u)}{\partial m} = r\,F_\theta$$

(λ = ángulo de lean tangencial). Es el mecanismo por el que el **compound lean** entra en un through-flow sin cálculo álabe-a-álabe auxiliar. Formulación explícita para turbinas: Wang & Xu, *The efficient modeling of blade lean effects within the turbomachinery throughflow method*, Int. J. Heat & Fluid Flow 10(3) [V-ABS]; Denton & Xu (1999), *The Exploitation of 3D Flow in Turbomachinery Design*, VKI LS 1999-02 [V-ABS].

### 1.3 Diferencias clave TURBINA vs COMPRESOR (y su consecuencia numérica)

**(a) Flujo acelerado → el lazo exterior converge mejor.** El cierre β₂(C_m) del compresor es rígido; en turbina el ángulo de salida está fijado casi geométricamente por la garganta, así que el acoplamiento C_m↔C_u es blando. *Predicción: `SCM_RELAX_R` puede subir de 0.35 a ~0.5–0.6 y `BLADE_RELAX` de 0.35 a ~0.6.* [DOM]

**(b) Pero aparece la doble solución sub/supersónica.** La continuidad ρ(C_m)·C_m·A tiene **dos raíces** y ρC_m tiene un máximo en M_m=1. La bisección de `_solve_station_cm` de Phy-AC **no está bien planteada cerca del bloqueo**: sin raíz por encima del gasto máximo, dos por debajo. Es el problema que Tiwari, Stein & Lin (2013) resuelven [V-ABS]. **Riesgo nº 1 de portar `scm_core.py` sin tocarlo.**

**(c) El caudal lo fija la primera tobera.** Bloqueada la primera tobera, ṁ√T₀/p₀ es constante. Ley de la elipse de Stodola (rama no bloqueada):

$$\frac{\dot m \sqrt{T_{0,in}}}{p_{0,in}} = K\sqrt{1-\left(\frac{p_{out}}{p_{0,in}}\right)^{2}}$$

(refinamiento: Cooke 1985, DOI 10.1115/1.3239778 [REF]). Esto invierte la lógica L0↔L1: **L1 puede DESCUBRIR que el gasto de L0 es imposible** (excede el máximo bloqueado de la garganta que L0 dimensionó). No es divergencia: es información, y hay que reportarla como tal.

**(d) Mucho más giro por fila** (tobera 65–75°, rotor 90–110° vs 20–40° en compresor). El término −(C_u/r)∂(rC_u)/∂r es dominante → el perfil C_m(r) está más determinado por la ley de vórtice; pero el gradiente radial de C_u dentro de la fila es enorme → **las estaciones internas de fila importan más** (Frost 1970, Casey & Robinson 2010).

**(e) Anillos con flare fuerte → la curvatura meridional deja de ser pequeña.** En una LPT el anillo se abre >25–30° de semi-ángulo en el shroud. `CURV_MAX = 8.0` de Phy-AC (que allí protege contra ruido) **en una turbina de BP recortaría física real** — sustituir por suavizado de la línea (filtro sobre r''), no clamp sobre 1/r_c.

**(f) Gradientes radiales de T impuestos, no calculados.** La entrada de la primera tobera es el perfil de salida de la cámara (OTDF/RTDF, pico 1.05–1.15× la media); la refrigeración inyecta masa fría a radios específicos. **∂h₀/∂r y T∂s/∂r son términos de primer orden impuestos desde fuera.** Un SCM de turbina que no acepta perfil radial de T₀ de entrada no sirve.

**(g) Gas caliente**: γ≈1.30–1.33, cp≈1150–1250 J/kg·K, R modificado por dosado y por el aire de refrigeración. Sustituir cp(T) por **cp(T, FAR)** y R(FAR); gasto y composición cambian estación a estación.

**(h) Mach absoluto alto en el marco ABSOLUTO**: la salida de tobera está a M₂ = 0.85–1.3 rutinariamente; con vórtice libre, el M₂ del cubo puede ser 0.3 mayor que el de la punta. La pérdida de onda/post-expansión se evalúa **por línea de corriente, en el marco de cada fila**.

---

## 2. Leyes de vórtice en turbinas

### 2.1 Vórtice libre (FVD)

r·C_u = const, h₀ = const(r), C_m = const(r). Solución exacta de la ecuación simplificada. Se usa más que en compresores (trabajo uniforme en el span es deseable). Desventaja documentada [V-ABS]: mucha torsión de raíz a punta; y con U∝r la reacción cae hacia el cubo:

$$R(r) = 1 - \frac{1-R_m}{(r/r_m)^2}$$

Con HTR ≲ 0.6 **la reacción en el cubo se hace negativa** → separación y pérdida secundaria masiva. El límite duro del vórtice libre. [DOM]

### 2.2 Vórtice general / exponencial

C_u = a·rⁿ (n=−1 libre). Forma cerrada:

$$C_m^2(r) = C_{m}^2(r_m) - \frac{n+1}{n}\,a^2\left(r^{2n}-r_m^{2n}\right) \quad (n\neq 0)$$

**Es la función `physics_core.vortex_cx` de Phy-AC (test T21): se reutiliza sin cambios** — agnóstica de si la máquina comprime o expande.

### 2.3 Ángulo de tobera constante (CNA)

α₂ = const(r) ⟹ C_u ∝ r^(−cos²α₂). Motivación de fabricación (NGV sin torsión) [V-ABS, Adiwidodo et al., J. Phys. Conf. Ser. 1005:012026, 2018]. Con α₂≈70°, C_u ∝ r^(−0.117): reacción mucho más plana en el span. El compromiso clásico de turbinas pequeñas.

### 2.4 Vórtice controlado (CVD) y controlado por presión (PCVD)

Se prescribe Δh₀(r) y/o C_m(r). PCVD [V-ABS]: controla la velocidad axial y la presión radial en el hueco estátor-rotor para **reducir pérdidas secundarias** (*Using a pressure controlled vortex design method…*, Chinese J. Aeronautics 2013). Comparativa de tres leyes: Springer LNME 2020 [V-ABS].

### 2.5 Compound lean a nivel through-flow

Álabe inclinado en sentidos opuestos en cubo y punta ("C"): la componente radial de la fuerza empuja el flujo hacia el midspan en ambos extremos, suprimiendo el vórtice de pasaje [V-ABS]. Entra por el término F_r = −F_θ·tanλ(r) de §1.2. Ver Harrison, PhD Cambridge ([repository](https://www.repository.cam.ac.uk/items/124d7c25-babe-4d40-a92b-7c2d590f4c44)).

**Recomendación Phy-AT L1**: no implementar compound lean en v1, pero **escribir el término F_r con tanλ(r) ≡ 0** por defecto y documentado — añadirlo después es rellenar un cero, no refactorizar.

---

## 3. Cierre de pérdidas y ángulos en el span

### 3.1 El ángulo de salida NO se cierra con Carter — se cierra con la garganta

**La diferencia estructural más importante entre un SCM de compresor y uno de turbina.**

$$\boxed{\;\cos\alpha_2 = \frac{o}{s}\;}\qquad\text{(regla de la garganta, α desde el eje; exacta en } M_2 = 1\text{)}$$

**Advertencia de convención**: algunos textos miden desde la tangencial (sin α₂ = o/s). Fijar la convención una vez y auditar cada fórmula importada.

**Forma robusta con flare (verificada en código)**: cociente de ÁREAS — el **gauging angle**:

$$\beta_g = \arccos\!\left(\frac{A_{\text{throat}}}{A_{\text{out}}}\right)$$

**[V-CÓDIGO]** — transcrito de `turboflow/axial_turbine/deviation_model.py`:
```python
gauging_angle = math.arccosd(geometry["A_throat"] / geometry["A_out"])
```
Absorbe automáticamente el flare del anillo y el bloqueo del TE.

**Desviación subsónica** — dos correlaciones verificadas en código:

*Ainley & Mathieson (1951)* [V-CÓDIGO], tres tramos con rodillas en M=0.5 y M_crit; aviso del propio código: inexacto si β_g > 70°, sin sentido si β_g > 72°.

*Aungier (2006)* [V-CÓDIGO], **C¹** (interpolante quíntico 1−10X³+15X⁴−6X⁵):
$$\delta_0 = \sin^{-1}\!\left[\frac{A_t}{A_o}\left(1+\left(1-\frac{A_t}{A_o}\right)\left(\frac{\beta_g}{90}\right)^{2}\right)\right]$$

**Para Phy-AT, usar Aungier por defecto** (las rodillas de AM producen chattering en el lazo exterior).

### 3.2 Corrección supersónica: expansión post-garganta

Para M₂ > 1, por continuidad (no correlación):

$$\boxed{\;\cos\alpha_2 = \frac{\dot m^{*}}{\rho_2\,C_2\,A_{\text{out}}}\;}$$

**[V-CÓDIGO]** (`choking_criterion.py`). Como ρC tiene máximo en M=1, para M₂>1: cosα₂ > o/s → **α₂ disminuye al crecer M₂** (desviación supersónica) [V-ABS].

### 3.3 Reparto radial de las pérdidas: el corazón del L1

#### 3.3.1 El esquema de Benner et al. (2006) — el correcto para el span

Benner, Sjolander & Moustapha (2006), J. Turbomach. 128(2):273–280 (Parte I) y 281–291 (Parte II) [REF; V-CÓDIGO].

Idea clave: **el vórtice de pasaje "come" una fracción del span** donde ya no hay flujo de perfil 2D; la **profundidad de penetración** Z_TE/H se correlaciona:

$$\frac{Z_{TE}}{H} = \frac{0.10\,F_t^{0.79}}{\sqrt{CR}\;(H/c)^{0.55}} + 32.70\left(\frac{\delta^{*}}{H}\right)^{2},\quad CR=\frac{\cos\beta_{in}}{\cos\beta_{out}},\quad \frac{\delta^{*}}{H}=\left(\frac{\delta^{*}}{H}\right)_{ref}\!\left(\frac{Re_{in}}{3\times10^{5}}\right)^{-1/7}$$

Uso [V-CÓDIGO]:
```python
Y_p  *= 1 - ZTE      # el perfil solo actúa fuera de la zona del vórtice
Y_te *= 1 - ZTE
```

Secundaria [V-CÓDIGO], dos ramas por relación de aspecto:

$$AR \le 2:\ Y_s = \frac{0.038 + 0.41\tanh(1.20\,\delta^{*}/H)}{\sqrt{\cos\xi}\;CR\;(H/c)^{0.55}\left(\frac{\cos\beta_{out}}{\cos\xi}\right)^{0.55}}\qquad AR > 2:\ Y_s = \frac{0.052 + 0.56\tanh(1.20\,\delta^{*}/H)}{\sqrt{\cos\xi}\;CR\;(H/c)\left(\frac{\cos\beta_{out}}{\cos\xi}\right)^{0.55}}$$

**Interpretación para Phy-AT: la banda de pared no es `WALL_BAND_FRAC = 0.30`. Es Z_TE/H, y sale de la carga, el Reynolds y la relación de aspecto de esa fila** (~0.10 en rotor largo de BP a >0.35 en tobera corta de HP). Un salto real de fidelidad.

#### 3.3.2 Pérdida de perfil y choque: Kacker & Okapuu (1982)

El término de choque de KO **ya es función del radio**:

$$Y_{shock} = 0.75\,[\max(0,\ f_{hub}\,M_{1,rel} - 0.4)]^{1.75}\;\frac{r_h}{r_t}\;\frac{p_{01,is}-p_1}{p_{02}-p_2}$$

con f_hub tabulado en r_h/r_t (clamp inferior en 0.5) [V-CÓDIGO]. **KO admite que el choque nace en el CUBO** y lo corrige al midspan porque es un meanline. **Un SCM no necesita el factor: evalúa el choque por línea de corriente donde de verdad ocurre.** También: K_p = 1−K₂(1−K₁), f_Ma = 1+60(M_out−1)² solo para M_out > 1, **por línea de corriente** (cubo a M₂=1.15, punta a 0.85 es normal).

#### 3.3.3 Fuga de punta: con y sin banda (shroud)

Modelo KO/Dunham-Came [V-CÓDIGO]: Y_cl = B·Z·(c/H)·(t_cl/H)^0.78, B = 0 estátor, 0.37 rotor CON banda; ~0.47 sin banda [DOM — verificar].

Física distinta [V-ABS]: sin banda la fuga es un chorro circunferencial (vórtice de punta, mezcla en el 15–30 % exterior + **defecto de trabajo local**); **con banda la fuga puentea la fila entera** por encima del shroud → **modelar como fuente/sumidero de masa**, no como Δs en banda. Refs: Yaras & Sjolander (1992) J. Turbomach. 114(1):204–210 [REF]; Denton (1993); Pacciani et al. (vórtices de Lamb-Oseen + fuentes/sumideros) [V-ABS].

#### 3.3.4 Familias de pérdidas de turbina (mapa)

| Familia | Referencia | Carácter |
|---|---|---|
| Ainley & Mathieson (1951) | ARC R&M 2974 (datos: R&M 2891) | base histórica |
| Dunham & Came (1970) | J. Eng. Power 92(3):252–256 | +Reynolds, +AR |
| **Kacker & Okapuu (1982)** | J. Eng. Power 104(1):111–119 | estándar industrial; choque y compresibilidad |
| Moustapha, Kacker & Tremblay (1990) | J. Turbomach. 112(2):267–276 | **incidencia off-design** |
| **Benner, Sjolander & Moustapha (2006)** | J. Turbomach. 128(2):273–291 | **desglose + penetración; reparto spanwise nativo** |
| Craig & Cox (1970) | Proc. IMechE 185(1):407–424 | familia alternativa (vapor) |
| Traupel (2001) | *Thermische Turbomaschinen*, Springer | escuela alemana/suiza |
| Zhu & Sjolander (2005) | GT2005-69077 | perfil + desviación actualizados |
| Coull | *An Improved Correlation for Turbine Endwall Loss*, J. Turbomach. 148(2):021008 | pared, actualización reciente |
| Revisión AM 2022 | *A Reliable Update of the AM Profile and Secondary Correlations*, IJTPP (open access) | revisión moderna |

### 3.4 Cómo entra la refrigeración en un through-flow

**Nivel 0 — Contabilidad conservativa (obligatorio).** Por fila k, ṁ_c a T₀c, p₀c, en radio prescrito:

$$\dot m_2 = \dot m_1 + \dot m_c;\quad \text{ESTÁTOR: } h_{02} = \frac{\dot m_1 h_{01} + \dot m_c h_{0c}}{\dot m_2};\quad \text{ROTOR: } I_2 = \frac{\dot m_1 I_1 + \dot m_c I_c}{\dot m_2},\ I_c = h_c + \tfrac12 W_c^2 - \tfrac12(\Omega r_c)^2$$

**La forma rotálpica no es cosmética**: el refrigerante entra por el disco a bajo radio y sale a alto radio — absorbe **trabajo de bombeo centrífugo** ½Ω²(r_out²−r_in²).

**Nivel 1 — Entropía de mezcla (Hartsel / Young & Wilcock).** Hartsel (1972), AIAA 72-11: mezcla **a presión estática constante** [V-ABS]. Young & Wilcock (2002), J. Turbomach. 124: separación de la generación de entropía en **término térmico** (igualar temperaturas) y **término cinético** (igualar velocidades) [V-ABS]:

$$\Delta s_{mix} \approx \varepsilon\,c_p\!\left[\ln\frac{T_g}{T_c}-\left(1-\frac{T_c}{T_g}\right)\right] + \varepsilon\,\frac{|\vec C_g - \vec C_c|^{2}}{2\,T_g} + \dots \quad\text{[DOM, forma]}$$

**El término cinético depende del VECTOR de velocidad** → una inyección a contracorriente es mucho peor que una alineada. **Un SCM lo captura porque conoce C⃗ en cada línea; un meanline no.** Segunda gran razón para tener L1 en una turbina.

**Nivel 2 — Lazo de refrigeración (modelo m*)**: LUAX-T (Lund) — Genrup et al. (2005) GT2005-68716; Sammak, Thern & Genrup (2013) GT2013-95469 [REF]. Arquitectura de **tres lazos anidados: refrigeración / entropía / geometría** [V-ABS] — la diferencia estructural con Phy-AC (dos lazos). Estado del arte: Petrović & Wiedermann (2013); Li, Gu & Song (2015) DOI 10.1177/0957650915594294; Ba et al. (2018) DOI 10.1177/0957650917731629 [V-ABS].

---

## 4. Choking y flujo transónico en through-flow de turbinas

### 4.1 Por qué es el problema central

$$\frac{\partial(\rho C_m)}{\partial C_m} = \rho\left(1 - M_m^{2}\right)$$

cambia de signo en M_m=1: dos soluciones para ṁ < ṁ_max, ninguna por encima, y derivada nula en el límite (Newton/secante mal condicionados).

### 4.2 Las estrategias de los códigos reales

**(A) Garganta como plano virtual + limitación de gasto por tubo de corriente** (Casey & Robinson 2010 [V-ABS]): evaluar el bloqueo por streamline en la garganta y **redistribuir el flujo en el span** — exactamente lo que un meanline no puede hacer. Argumento de venta del L1 de turbina.

**(B) Mach crítico corregido por pérdidas (<1)** [V-CÓDIGO, `get_mach_crit`]:

$$M_{crit}^{2}=\frac{2}{\gamma-1}\left[\frac{4\alpha-2}{(2\alpha+\eta-3)+\sqrt{(1+\eta)^{2}+4\alpha(1+\alpha-3\eta)}}-1\right],\qquad \alpha=\frac{\gamma}{\gamma-1}$$

con η el rendimiento de la fila hasta la garganta (η=1 → M_crit=1 exacto; η=0.95, γ=1.33 → M_crit≈0.95–0.97).

**(C) Optimización explícita del gasto en garganta (KKT/Lagrange)** — formulación equation-oriented de TurboFlow [V-CÓDIGO]; Anderson et al. (2024), J. Turbomach. [REF]. TurboFlow expone los tres criterios como opciones — **lección de diseño: el criterio de bloqueo es incertidumbre modelística explícita**.

**(D) Guiar la rama sub/supersónica** — Tiwari, Stein & Lin (2013) [V-ABS].

### 4.3 La primera tobera fija el caudal

Ley de Stodola (rama no bloqueada); Cooke (1985) [REF]; para pocas etapas, la "semi-elipse" (*Enhancement to the Traditional Ellipse Law…*, J. Eng. Gas Turbines Power 139(11):112603 [V-ABS]).

### 4.4 Salidas supersónicas (M₂>1)

Ángulo por continuidad (§3.2) + penalización post-garganta. Datos: Graham & Kost (1979), ASME 79-GT-37 [REF]; Martelli & Boretti (1987) [REF]. **Carga límite** (limit loading): Chen (2018), NASA — por encima, más relación de expansión no produce trabajo adicional → modo de fallo honesto [REF].

---

## 5. Códigos y sistemas existentes: qué resuelve cada uno y qué lección aporta

### 5.1 Abiertos / públicos

**MULTALL-OPEN (Denton)** — [sitio](https://sites.google.com/view/multall-turbomachinery-design). MEANGEN (meanline) → STAGEN (geometría) → MULTALL (NS 3D multietapa + Q3D + throughflow axisimétrico) [V-ABS]. Dominio público. GT2017-63993 / J. Turbomach. 139(12):121001 [REF]. **Lecciones**: separación estricta meanline→geometría→solver con ficheros de intercambio (la escalera de Quasar); publicado "para individuos o pequeñas empresas sin sistema in-house"; solver deliberadamente barato — el valor está en el ciclo de diseño.

**T-AXI / T-C_DES (Turner, U. Cincinnati)** — [web](https://www.downingdesign.net/gtsl/codes/taxi/). Axisimétrico de diseño, **turbinas y compresores**, construido sobre **MTFLOW de Drela**: Newton totalmente acoplado con la posición de las líneas como incógnita — **no sufre la doble raíz** porque no bisecciona la continuidad [V-ABS]. Ref: Turner et al. (2011), J. Turbomach. 133(3):031017. **Lección: si el lazo tipo Novak resulta frágil en transónico, la alternativa madura es Newton acoplado, no más relajación.**

**turbigen (Brind, Whittle Lab)** — [turbigen.org](https://turbigen.org/). **Sí soporta turbinas** (`meanline/axial_turbine.py`) [V-ABS]. Abierto. **Lección incómoda**: turbigen usa CFD 3D rápido en GPU para **saltarse meanline empírico, throughflow y Q3D**. La respuesta de Phy-AT debe ser cuantitativa: el L1 existe para que el ensemble residual tenga un residual que aprender a ~3 s/máquina y sin dependencias, que es un requisito distinto.

**TurboFlow (DTU/NTNU)** — [GitHub](https://github.com/turbo-sim/turboflow), JOSS 10(111):7588 (2025). Meanline de turbinas equation-oriented, submodelos intercambiables de pérdidas/desviación/choking, gas real. **Fuente verificada línea a línea de §3–4.** Lecciones: consistencia diseño/análisis; criterios de choking como opciones explícitas; **los clamps marcados como TODO en el propio fuente** (disciplina `SCMDiverged`); configuración YAML.

**AXOD/AXOD2 (Glassman, NASA)** — [NTRS 19950004441]; TD2-2; Flagg (1967) GE R66FPD258 [REF]. Off-design meanline, dominio público. **SP-290 (Glassman) sigue siendo la referencia pedagógica libre.**

**OTAC (NASA GRC)** — meanline/streamline sobre NPSS, con **equilibrio radial simple** (sin curvatura — Phy-AC ya es más completo). **Hendricks (2016), AIAA 2016-0119**: NASA descubrió la limitación de flujo bloqueado en turbinas *después* de construir OTAC [V-ABS/REF]. **Lección accionable: el tratamiento del choking es requisito de la v1, no extensión posterior.**

**LUAX-T (Lund)** — tres lazos anidados (refrigeración/entropía/geometría); AM modificado + modelo m* de refrigeración [V-ABS]. Tesis de Dahlquist (2008) — comparativa de modelos de pérdidas [REF].

### 5.2 Comerciales / industriales

| Sistema | Lección |
|---|---|
| **PCA SC90T / SC90C / Vista TF** ([web](https://www.pcaeng.co.uk/vista-software)) | **Códigos GEMELOS turbina/compresor desde 1990, no uno parametrizado** → `scm_core.py` y el SCM de Phy-AT deben ser módulos hermanos con utilidades compartidas, no un `if machine_type`. Mismo código en modo análisis y modo cuasi-diseño (prescribiendo rC_u) |
| **AxSTREAM (SoftInWay)** | Escalera 1D→2D streamline→through-flow explícita; **hasta 49 estaciones radiales** para álabes largos (vs 9 de Phy-AC — rehacer el estudio de independencia con flare); "ángulos reales del álabe en cada posición del span" = la filosofía de scm_core |
| **Concepts NREC AXIAL/AxCent/Agile** | El traspaso 1D→3D como artefacto de primera clase — el problema de deriva de `PR_WINDOW` de Phy-AC es un fallo de ese contrato |
| **SOCRATES (RR/Cranfield)** | El SLC sigue siendo "la piedra angular del diseño por razones económicas y prácticas"; Pachidis et al. (2010), *dynamic convergence control*, DOI 10.1243/09544100JAERO565 — relevante a los factores de relajación |

### 5.3 Notas

- **TD3/turbo-design**: Phy-AC lo usó como L1 hasta la fase 11 y lo descartó (ODE de equilibrio radial divergente con >1 streamline). **No repetir el experimento en turbina.** Sí es referencia de modelos de pérdidas de turbina (TD2, AM, KO, Craig-Cox, Traupel implementados).
- **zTurbo (VKI)**: no localizado con el presupuesto disponible. Hueco.
- Dubitsky, Wiedermann et al. (2003), *The Reduced Order Through-Flow Modeling of Axial Turbomachinery* [REF].

---

## 6. Métodos through-flow modernos

### 6.1 Through-flow CFD (Euler/NS axisimétrico con fuerzas de álabe)

Origen: campo de fuerzas de Marble (1964), usado por Denton (1978), Hirsch & Warzee (1976) [V-ABS].
- **Simon & Léonard** (2005, 2007, 2009): NS axisimétrico con fuerzas de álabe y bloqueo; *On the Role of the Deterministic and Circumferential Stresses in Throughflow Calculations*, J. Turbomach. 131(3):031019 — **cuantifica lo que el SCM promediado pierde** [V-ABS].
- **Pacciani et al.** (2016, Proc. IMechE A, DOI 10.1177/0957650915607091; IJTPP 2(3):11 open access): Euler axisimétrico con esquema de TRAF, superficie S2 adaptativa, secundarios y fuga como vórtices de Lamb-Oseen + fuentes/sumideros; aplicado a la **HPT transónica PW E³ y a una LPT** [V-ABS].
- Euler-2D para turbinas multietapa: ETC2021-661; *Time-marching throughflow model…*, J. Braz. Soc. Mech. Sci. Eng. (2023) [V-ABS].

**Ventaja del time-marching**: el bloqueo y el supersónico salen solos (sin doble raíz). **Coste**: dos órdenes de magnitud más que un SCM → descartado para L1 (~3 s/máquina), pero destino natural de un futuro L1.5.

### 6.2 Through-flow con deep learning embebido (2020–2026)

- Revisión: *AI in turbomachinery aerodynamics*, AI Review (2024), DOI 10.1007/s10462-024-10867-3 [V-ABS].
- Transfer learning de modelos de pérdidas (helio ← gas): Energy (2023/24) [V-ABS] — **el patrón de la capa 2 de Quasar publicado independientemente**.
- **Pérdida spanwise con ML**: *A Spanwise Loss Model of Turbine Cascade with Tip Clearance Based on ML*, Springer LNME (2024) [V-ABS] — predice el perfil radial de pérdida, directamente aplicable a L1.
- **Desviación secundaria física para throughflow**: *A physically based semi-empirical model for secondary loss in throughflow analysis…*, Phys. Fluids 37(4):046113 (2025) — magnitud Y **distribución en el span**, descompuesta en vórtice de esquina y capa límite de pared [V-ABS]. **Estado del arte 2025 para §3.3.**

### 6.3 Through-flow como generador de datos multifidelidad

Usos: (1) anclaje — el residual L1−L0 es lo que aprende la capa 2; (2) barridos masivos; (3) calibración afín L2. Contra-tesis de turbigen presente: **la respuesta debe ser cuantitativa** (medir que el residual L1−L0 tiene estructura aprendible).

---

## 7. Verificación

### 7.1 Soluciones analíticas

- **(V1) Vórtice libre exacto**: C_m(r) estrictamente uniforme (error < 1e−5). Atrapa signos en −(C_u/r)∂(rC_u)/∂r.
- **(V2) Vórtice general** C_u=arⁿ: forma cerrada de §2.2 — **reutilizar el test T21 de Phy-AC sin cambios**; añadir el caso de ángulo constante n=−cos²α₂ con α₂=70°.
- **(V3) Actuator disc**: perturbación de C_x decae como e^(−π|z|/h), mitad del cambio en el plano del disco [DOM; familia de soluciones en Bessel confirmada V-ABS]. **EL test del término de curvatura meridional** — imprescindible antes de tocar un anillo con flare.
- **(V4) Conservación exacta**: rotalpía por rotor adiabático; h₀ por estátor; con refrigeración, balance de masa y energía cerrado a máquina redonda; la mezcla spanwise conserva la media másica.
- **(V5) Independencia de discretización**: N_sl ∈ {5,7,9,11,15,21} — rehacer con flare fuerte.
- **(V6) Tobera bloqueada**: sin pérdidas, el gasto debe reproducir
$$\dot m^{*} = A_{throat}\,p_0\sqrt{\frac{\gamma}{RT_0}}\left(\frac{2}{\gamma+1}\right)^{\frac{\gamma+1}{2(\gamma-1)}}$$
y el ángulo arccos(o/s) exacto. Un test de una línea que atrapa la mayoría de errores de §4.

### 7.2 Validación — bancos públicos

| Caso | Referencia | Qué aporta |
|---|---|---|
| Turbina 1 etapa turbofán pequeño | **Kofskey & Nusbaum (1972), NASA TN D-6967** — [NTRS 19720024422](https://ntrs.nasa.gov/citations/19720024422) | Geometría completa + prestaciones en aire frío. Ancla L1 |
| Turbina turborreactor bajo coste | **Kofskey, Roelke & Haas (1974), NASA TN D-7625** — [NTRS 19740018139](https://ntrs.nasa.gov/citations/19740018139) | Diseño + aire frío |
| Ídem, **estátor abierto** | Kofskey, Nusbaum & Haas (1974), NASA E-7776 — [NTRS 19740019165](https://ntrs.nasa.gov/citations/19740019165) | **Barrido de área de garganta** ⇒ verifica choking y Stodola |
| Cascada off-design | Tremblay, Sjolander & Moustapha (1990), 90-GT-314 | Pérdidas por incidencia |
| Tobera transónica baja AR | Moustapha, Carscallen & McGeachy (1993), J. Turbomach. 115(3):400–408 | Peor caso para el reparto spanwise |
| Cascadas transónicas alto giro | Graham & Kost (1979), 79-GT-37 | Choque/capa límite |
| **PW E³ HPT** | Pacciani et al. | Comparación con through-flow CFD publicado |

**Nota honesta**: **no existe un "benchmark de through-flow de turbina" normalizado con datos radialmente resueltos + geometría + resultados de referencia de otros códigos.** Los anclajes L1 serán los rigs de Kofskey más lo reproducible de la literatura; declararlo en la documentación.

---

## 8. Recomendaciones concretas para Phy-AT L1

### 8.1 Arquitectura: `scm_turbine.py`, módulo HERMANO de `scm_core.py`

**Dos módulos, no uno parametrizado** (lección PCA SC90T/SC90C). Compartir `scm_common.py`: `_ddr`, `_curvature`, `_span_mix`, `_reposition`, `_streamtube_ratios`, `_wall_band_weight`, estado termodinámico.

### 8.2 Estaciones

De 2n+1 a **3n+1 (o 4n+1): añadir una estación de GARGANTA por fila** (Denton 1978; Casey & Robinson 2010). Secuencia por fila: LE → GARGANTA → TE; estaciones internas adicionales con flare fuerte. Las de garganta llevan A_throat(r) y A_out(r). Misma geometría que se fabrica (como Phy-AC).

### 8.3 Cierre del ángulo de salida: **o/s, NUNCA Carter**

```
α₂(r) = arccos( A_throat(r) / A_out(r) ) − δ(M₂(r), β_g(r))    si M₂ ≤ M_crit
α₂(r) = arccos( ṁ*(r) / (ρ₂(r)·C₂(r)·A_out(r)) )               si M₂ > M_crit
```
δ por Aungier (C¹); cociente de ÁREAS, no o/s escalar; guarda dura β_g > 72° → aviso; **fijar la convención de ángulos UNA vez** (desde el eje). Lo que se congela es la **garganta** A_throat(r) (análogo del "álabe fijo en el espacio" de Phy-AC).

### 8.4 Manejo del choking (requisito de v1)

1. **M_crit corregido por pérdidas** (§4.2B) por línea de corriente — una línea de código.
2. **Límite de gasto POR TUBO DE CORRIENTE en la garganta** con redistribución del exceso a las líneas no saturadas (Casey & Robinson 2010) — la capacidad que el meanline no tiene.
3. **Selección de rama**: bisección acotada a M_m ≤ M_crit (rama subsónica); si bloqueada, fijar M_crit en la garganta y resolver la salida en rama supersónica; **registrar qué filas bloquean y en qué fracción del span** (información de diseño de primera clase).
4. **Coherencia con L0**: si el ṁ pedido excede Σṁ*_j de la primera tobera → punto **infactible**, no divergente (§8.8).

### 8.5 Dónde va cada pérdida en el span

| Mecanismo | Modelo | Ubicación radial |
|---|---|---|
| Perfil | KO (1982) por línea, con Re, M, t/c locales | Todo el span **× (1−Z_TE/H)** |
| Borde de fuga | KO / Denton (1993) | Ídem × (1−Z_TE/H) |
| Choque/post-garganta | Y_shock de KO **sin f_hub** + f_Ma=1+60(M₂−1)² | Solo donde M local lo justifica (cubo de tobera HP, punta de rotor) |
| Secundaria | Y_s de Benner Parte II | **Bandas de espesor Z_TE/H de la correlación de penetración — NO constante** |
| Fuga de punta sin banda | Y_cl DC (B≈0.47) + Yaras & Sjolander | Banda de punta (15–30 %) **+ defecto de trabajo local** |
| Fuga con banda (shroud) | Denton (1993) | **Fuente/sumidero de masa que puentea la fila**, no Δs en banda |
| Incidencia | Moustapha (1990) | Por línea de corriente |
| Refrigeración | §8.6 | En el radio de inyección |

Conservar íntegra la mezcla spanwise (`_span_mix`) — referencia turbina: **Lewis (1994), *Spanwise Transport in Axial-Flow Turbines: Part 2*, J. Turbomach. 116:187–193** [V-ABS] — y **subir su importancia** (los gradientes de T₀ por refrigeración son más violentos que los de entropía de pérdida).

### 8.6 Refrigeración: tres lazos anidados (LUAX-T)

```
LAZO EXTERNO   (geometría/líneas)
  LAZO MEDIO   (entropía → equilibrio radial → Cm(r))
    LAZO INTERNO (refrigeración: T_gas(r) → ṁ_c → T_gas(r) …)
```
Cuatro puntos innegociables: (1) **rotalpía, no entalpía, en el rotor** (bombeo centrífugo del refrigerante); (2) Δs_mix separado en térmico + cinético (Young & Wilcock 2002) — el cinético depende del **vector** de velocidad; (3) mezcla a presión estática constante (Hartsel 1972); (4) composición variable: cp(T, FAR), R(FAR) estación a estación. La primitiva `inject(station, ṁ, h₀|I, s, r_distribution)` sirve también para el bypass de shroud y las purgas de cavidad.

### 8.7 Qué reutilizar de Phy-AC y qué cambiar

| Componente | Veredicto | Motivo |
|---|---|---|
| Ecuación de equilibrio radial completa | REUTILIZAR (estátor) | misma ecuación |
| Forma rotálpica para el rotor | AÑADIR | h₀ no se conserva; I sí |
| Término F_r = −F_θ·tanλ | AÑADIR con λ≡0 | hueco para compound lean |
| 9 líneas por fracción de masa | REUTILIZAR, rehacer estudio de independencia | flare (AxSTREAM llega a 49) |
| `_curvature` + curvatura retrasada | REUTILIZAR | Novak 1967 |
| `CURV_MAX = 8.0` | **CAMBIAR/eliminar** | en LPT recortaría física real → suavizar r'', no clampear 1/r_c |
| `_solve_station_cm` bisección | **CAMBIAR** | doble raíz cerca del bloqueo → acotar a M_m≤M_crit + rama supersónica |
| `_beta2_metal` de ley de vórtice | **SUSTITUIR** por arccos(A_t/A_o) | §8.3 |
| `_row_loss_span` (Koch&Smith/Howell/Lieblein) | **SUSTITUIR ENTERO** | modelos de compresor → KO + Benner + Moustapha |
| `_wall_band_weight` con constante | REUTILIZAR función, banda = Z_TE/H | Benner |
| `_tip_clearance_span` | REUTILIZAR sin banda; añadir rama fuente/sumidero con banda | física distinta |
| `_endwall_span` (Koch & Smith) | **SUSTITUIR** | el débito de pared lo lleva Y_s de Benner |
| `_span_mix` | REUTILIZAR y subir importancia | Lewis 1994 |
| `SCMDiverged` + degradación etiquetada | REUTILIZAR — lo mejor del diseño | §8.8 lo extiende |
| `PR_WINDOW` | REUTILIZAR concepto, magnitud = **gasto reducido** | la garganta fija el gasto; el problema será PEOR |
| `physics_core` cp(T) aire | **EXTENDER** a cp(T, FAR), R(FAR) | gas caliente + dilución |
| Test T21 | REUTILIZAR SIN CAMBIOS | agnóstico compresión/expansión |
| — | AÑADIR test actuator disc (V3) y tobera bloqueada (V6) | curvatura y choking |

### 8.8 Modos de fallo honesto específicos de turbina

| Excepción / bandera | Disparador | Semántica | Acción |
|---|---|---|---|
| `TurbineChoked` | fila alcanza ṁ = Σṁ*_j | **NO es fallo** — estado normal de una HPT | resolver en rama bloqueada, registrar fila y fracción de span |
| `MassFlowInfeasible` | ṁ de L0 excede el máximo de la primera tobera | contradicción del contrato L0↔geometría↔L1 | degradar a L0 **etiquetado con el ṁ máximo alcanzable** |
| `LimitLoading` | Δh₀ deja de crecer con la expansión | límite físico (Chen 2018) | trabajo saturado, marcado |
| `SupersonicBranchLost` | oscilación sub/supersónica sin converger | Tiwari 2013 | SCMDiverged con razón |
| `DeviationModelOutOfRange` | β_g > 72° | el cierre o/s pierde sentido | aviso; degradar si >30 % del span |
| `NegativeReactionAtHub` | R(r_hub) < 0 | fallo del DISEÑO, no del solver | resultado con bandera roja — que el optimizador cambie la ley de vórtice |
| `CoolingLoopDiverged` | el lazo interno no cierra | realimentación ṁ_c↔T_gas | SCMDiverged con razón |
| `EnergyImbalance` | balance de energía no cierra | error de contabilidad de refrigeración | **assert duro, nunca degradación silenciosa** |
| `CurvatureClamped` | limitador activo al converger | mismo principio que CM_SPREAD | SCMDiverged |

Principio a preservar: distinguir **fallo numérico** (degradar a L0) de **infactibilidad física** (devolverla como información al optimizador).

### 8.9 Orden de implementación sugerido

1. Esqueleto: estaciones con garganta, equilibrio radial reusado, rotalpía, cierre arccos(A_t/A_o) sin desviación. **V1, V2, V6.**
2. Curvatura con flare + actuator disc **V3**. Eliminar CURV_MAX.
3. Choking: M_crit corregido, límite por streamtube, redistribución, selección de rama. **Validar contra Kofskey E-7776 (estátor abierto).**
4. Pérdidas: KO + Benner (Z_TE/H) + Moustapha. Validar contra Kofskey TN D-6967.
5. Desviación de Aungier.
6. Refrigeración: tres lazos. **V4 antes de cualquier prestación.**
7. Mezcla spanwise recalibrada.
8. Hueco de compound lean (λ≡0).

---

## Fuentes

(Listado completo de URLs por trazabilidad — muchas bloqueadas por el proxy en la sesión; ver hipervínculos in-line. Principales: Denton 1978 J.Eng.Power 100(2); Novak 1967; Smith 1966; Wu NACA TN 2604; Wilkinson 1969; Frost ARC R&M 3687; Casey & Robinson 2010 J.Turbomach 132(3); Tiwari/Stein/Lin 2013 135(4); Petrović & Wiedermann 2013 135(6); Benner 2006 I/II 128(2); Kacker & Okapuu 1982; Moustapha 1990; Young & Wilcock 2002; Hartsel 1972; Genrup 2005; Sammak 2013; Lewis 1994; Hendricks AIAA 2016-0119; OTAC NTRS 20205004138; Kofskey NTRS 19720024422 / 19740018139 / 19740019165; MULTALL; T-AXI; turbigen; TurboFlow JOSS 7588; PCA Vista; AxSTREAM; Concepts NREC; Pacciani IJTPP 2(3):11; Simon & Léonard J.Turbomach 131(3); Cooke 1985; Stodola.)
