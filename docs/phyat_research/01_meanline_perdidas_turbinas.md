# ESTADO DEL ARTE — DISEÑO MEANLINE (L0) DE TURBINAS AXIALES
## Fundamento científico para **Phy-AT (Quasar)**, análogo de Phy-AC

---

## 0. Alcance, método y trazabilidad de fuentes

### 0.1 Cómo se obtuvo este material

La red de la sesión de investigación bloquea casi todos los dominios académicos; los canales útiles fueron `raw.githubusercontent.com` y `github.com`. Consecuencia **favorable**: las ecuaciones no vienen de resúmenes de segunda mano sino de **implementaciones de referencia, abiertas, publicadas y validadas contra medida**:

| Fuente-código | Qué aporta | Estatus |
|---|---|---|
| **TurboFlow** (`turbo-sim/turboflow`), Anderson, Agromayor, Haglind & Nord — JOSS 2025, [DOI 10.21105/joss.07588](https://joss.theoj.org/papers/10.21105/joss.07588) | KO 1982, Moustapha 1990, Benner 1997/2004/2006 completos; desviación (AM, Aungier); choking; definiciones de coeficiente de pérdida; **validación NASA con tablas de error** | Verificado línea a línea |
| **AxialOpt** (`turbo-sim/AxialOpt`), Agromayor & Nord — [IJTPP 4(3):32, 2019](https://doi.org/10.3390/ijtpp4030032) | AM 1951 y DC 1970 completos; **restricciones g con justificación bibliográfica anotada**; correlaciones de stagger y espesor de Kacker | Verificado línea a línea |
| **NASA turbo-design (TD3)** (`nasa/turbo-design`) | Craig–Cox 1970 y Traupel — evidencia directa de que son *chart-based* (~18 figuras digitalizadas en `.pkl`) | Verificado |
| **pyCycle** (`OpenMDAO/pyCycle`) | Algoritmo de refrigeración fila a fila (Gauntner, NASA TM-81453) con constantes explícitas | Verificado |

Marcadores: **[V]** verificado sobre fuente primaria o implementación de referencia; **[S]** de resumen de búsqueda (cita correcta, número de segunda mano); **[?]** no verificado — comprobar antes de usar como ancla.

### 0.2 La tesis del informe en una frase

> Para Phy-AT, la cadena **AMDC-KO + Benner (2006) para secundaria/penetración + Benner (1997) para incidencia** es la única familia que simultáneamente (a) es **cerrada en forma algebraica** — vectorizable en NumPy puro, sin cartas —, (b) tiene **precisión publicada de ±1.5 % en η sobre 33 turbinas sin recalibración por máquina**, y (c) está **validada en abierto contra tres turbinas NASA medidas** con tablas de error reproducibles. Craig–Cox y Traupel, aunque más precisas en nichos, son inaplicables a un L0 de 0.5 ms/punto (sistemas de ~18 cartas interpoladas).

---

## 1. Formulación meanline de turbinas

### 1.1 Convenciones de signo y ángulo — la trampa nº 1

**Convención recomendada (TurboFlow) [V]:** ángulos desde la dirección **axial**, positivos antihorario:

$$v_m = v\cos\alpha,\quad v_\theta = v\sin\alpha,\quad w_\theta = v_\theta - u,\quad \beta = \arctan(w_\theta/w_m)$$

En una turbina de reacción: α₂ > 0 (salida de estátor), β₃ < 0 (salida de rotor). Las correlaciones de AM están tabuladas frente al ángulo **referido a la tangencial**: φ_tan = 90° − |β_ax| (el código de AxialOpt lo marca como *"tricky sign convention"* [V]). **Phy-AT debe encapsular la conversión en una única función y testearla**, como Phy-AC verifica tanβ₁+tanα₁=1/φ a 1e-6.

La pérdida de perfil de KO usa el cociente **con signo** con la construcción `abs(x)*x` para servir a estátor y rotor con el mismo código [V]:

$$Y_p = Y_{p,\text{reaction}} - \left|\frac{\theta_{in}}{\beta_{out}}\right|\left(\frac{\theta_{in}}{\beta_{out}}\right)\bigl(Y_{p,\text{impulse}} - Y_{p,\text{reaction}}\bigr)$$

### 1.2 Triángulos de velocidad a partir de (φ, ψ, R)

Estaciones: 1 = entrada estátor, 2 = salida estátor/entrada rotor, 3 = salida rotor. Con U y c_x constantes:

$$\varphi = \frac{c_x}{U},\qquad \psi = \frac{\Delta h_0}{U^2},\qquad R = \frac{h_2 - h_3}{h_{01} - h_{03}}$$

Euler en expansión: Δh₀ = U(c_θ2 − c_θ3) ⟹ ψ = (c_θ2 − c_θ3)/U; con c_x constante, (c_θ2+c_θ3)/2U = 1−R:

$$\boxed{\;\frac{c_{\theta 2}}{U} = (1-R) + \frac{\psi}{2},\qquad \frac{c_{\theta 3}}{U} = (1-R) - \frac{\psi}{2}\;}$$

**El espejo exacto del par de Phy-AC**: el compresor añade swirl, la turbina lo consume.

$$\tan\alpha_2 = \frac{(1-R)+\psi/2}{\varphi},\quad \tan\alpha_3 = \frac{(1-R)-\psi/2}{\varphi},\quad \tan\beta_2 = \frac{-R+\psi/2}{\varphi},\quad \tan\beta_3 = \frac{-R-\psi/2}{\varphi}$$

**Identidades para tests de verificación**:

$$\tan\alpha_2 - \tan\beta_2 = \tan\alpha_3 - \tan\beta_3 = \frac{1}{\varphi};\qquad \psi = \varphi(\tan\beta_2 - \tan\beta_3);\qquad R = -\frac{\varphi}{2}(\tan\beta_2 + \tan\beta_3)$$

Casos límite: R=0.5 → triángulos simétricos (β₃=−α₂); R=0 → impulso, sin caída de presión estática en rotor.

**Advertencia**: coexisten ψ=Δh₀/U² (Dixon & Hall, Smith 1965 — la de Phy-AC, rangos "1.0–2.5") y ψ=2Δh₀/U². **Fijarla y documentarla.**

### 1.3 Etapa repetitiva vs. etapas reales de HPT/LPT

- **HPT (1–2 etapas)**: R bajo en raíz; HTR alto (0.85–0.92); ψ alto; el estátor 1 casi siempre **estrangulado** en diseño.
- **LPT (3–7 etapas)**: R≈0.5; HTR bajo (0.5–0.65) decreciente; AR 3–6 → Reynolds bajo en crucero; álabes high-lift con Zweifel alto.
- La etapa repetitiva **se rompe en HPT** y se sostiene en LPT. **Recomendación: no imponerla** — usar las pendientes lineales de φ, ψ, R de la fase 8 de Phy-AC (patrón real: φ crece hacia atrás, ψ decrece, R crece).

### 1.4 Diagrama de Smith (1965)

Smith, S. F., *"A simple correlation of turbine efficiency"*, **J. Royal Aeronautical Society 69:467 (1965)** [S]. η_tt de ~70 etapas de banco Rolls-Royce vs (φ, ψ), **±2 %**, etapas de 50 % reacción. KO 1982 validó contra la carta de Smith además de las 33 turbinas [S].

**Rangos de práctica** (orientación, [?]):

| | φ | ψ | R (mean) | AR |
|---|---|---|---|---|
| HPT aeronáutica | 0.45–0.70 | 1.3–2.2 | 0.20–0.45 | 1.0–2.0 |
| LPT aeronáutica | 0.65–1.00 | 1.5–2.5 | 0.40–0.55 | 3.0–6.0 |
| Industrial/vapor | 0.4–0.7 | 0.8–1.6 | 0.45–0.55 | 1.5–4.0 |

**Uso en Phy-AT**: superponer la salida del optimizador sobre el diagrama de Smith como diagnóstico visual, **no** como modelo de pérdidas.

### 1.5 Stage-stacking con expansión y gas real

1. **Marcha de entropía en expansión** con cp(T) — la arquitectura de la fase 9 de Phy-AC se traslada íntegra.
2. **cp(T) hasta 1900 K y dependiente de FAR** — el ajuste JANAF de Phy-AC (250–1000 K) es insuficiente.
3. **La restricción dura no es el bombeo sino el choking** — no existe stall en flujo acelerado; la "línea de trabajo" se sustituye por la ley de Stodola.
4. **El caudal cambia estación a estación** por refrigeración — refactor desde el inicio, no parche.

### 1.6 Definiciones de coeficiente de pérdida — la trampa nº 2

De la documentación de TurboFlow [V]. **Para turbinas la referencia es la presión dinámica de SALIDA**:

$$Y=\dfrac{p_{0,in}-p_{0,out}}{p_{0,out}-p_{out}} \ \text{(turbina)}\qquad \Delta\phi^2 = 1-\left(\frac{v_{out}}{v_{out,s}}\right)^2 \qquad \zeta = \frac{1}{\phi^2}-1 \qquad \varsigma = \frac{T_{out}(s_{out}-s_{in})}{\tfrac12 v_{out}^2}$$

Conversión exacta Δφ² → Y para gas perfecto [V]:

$$Y = \frac{\left[1-\frac{\gamma-1}{2}\mathrm{Ma}_{out}^2\left(\frac{1}{1-\Delta\phi^2}-1\right)\right]^{-\gamma/(\gamma-1)}-1}{1-\left(1+\frac{\gamma-1}{2}\mathrm{Ma}_{out}^2\right)^{-\gamma/(\gamma-1)}}$$

TurboFlow/AxialOpt imponen **Y_definición = Y_modelo** como restricción de igualdad por cascada [V]. Phy-AT puede resolver la marcha explícitamente, pero **el residuo debe reportarse como diagnóstico**.

---

## 2. Modelos de pérdidas — el corazón del L0

### 2.1 Ainley & Mathieson (1951) — ARC R&M 2974 / 2891

$$Y = Y_p + Y_s + Y_{cl} + Y_{te}$$

**Perfil** [V]: interpolación reacción/impulso + corrección de espesor:

$$Y_p = \left[Y_{p,\text{reac}} - \left|\tfrac{\theta_{in}}{\beta_{out}}\right|\tfrac{\theta_{in}}{\beta_{out}}\left(Y_{p,\text{imp}} - Y_{p,\text{reac}}\right)\right]\left(\frac{t_{max}/c}{0.20}\right)^{a},\quad a=\max\!\left(0, -\frac{\theta_{in}}{\beta_{out}}\right)$$

Las curvas de AM son cartas; **Aungier (2006) publicó los ajustes analíticos** que usan TurboFlow y AxialOpt [V]. Con φ ≡ 90°−|β_out,ax| en grados:

*Tobera (reacción):*
$$\left(\tfrac{s}{c}\right)_{\min}=\begin{cases}0.46+\varphi/77 & \varphi<30\\ 0.614+\varphi/130 & \varphi\ge30\end{cases}\qquad X=\tfrac{s}{c}-\left(\tfrac{s}{c}\right)_{\min}$$
$$A=\begin{cases}0.025+\tfrac{27-\varphi}{530}&\varphi<27\\ 0.025+\tfrac{27-\varphi}{3085}&\varphi\ge27\end{cases}\quad B=0.1583-\tfrac{\varphi}{1640}\quad C=0.08\left[\left(\tfrac{\varphi}{30}\right)^2-1\right]\quad n=1+\tfrac{\varphi}{30}$$
$$Y_{p,\text{reac}}=\begin{cases}A+BX^2+CX^3&\varphi<30\\ A+B|X|^n&\varphi\ge30\end{cases}$$

*Impulso:*
$$\left(\tfrac{s}{c}\right)_{\min}=0.224+1.575\tfrac{\varphi}{90}-\left(\tfrac{\varphi}{90}\right)^2,\quad A=0.242-\tfrac{\varphi}{151}+\left(\tfrac{\varphi}{127}\right)^2$$
$$B=\begin{cases}0.30+\tfrac{30-\varphi}{50}&\varphi<30\\ 0.30+\tfrac{30-\varphi}{275}&\varphi\ge30\end{cases}\quad C=0.88-\tfrac{\varphi}{42.4}+\left(\tfrac{\varphi}{72.8}\right)^2\qquad Y_{p,\text{imp}}=A+BX^2-CX^3$$

> **Rango de validez explícito [V]**: *"válidas para 40 < |β_out| < 80; extrapolar puede dar resultados completamente erróneos"*. Y 0.30 ≤ s/c ≤ 1.10. **Restricción dura g en Phy-AT.** Salvaguarda: Y_p ≥ 0.80·Y_p,reac [V].

**Secundaria (AM)** [V]: Y_s = λ(x)·Z con λ tabulada (x=[0,…,0.5] → λ=[0.0055,…,0.0275], saturada).

**Parámetro de carga Zweifel-Ainley Z** (común a Y_s y Y_cl en toda la familia) [V]:

$$\boxed{\;Z = 4\,(\tan\beta_{in}-\tan\beta_{out})^2\,\frac{\cos^2\beta_{out}}{\cos\beta_m}\;},\qquad \beta_m=\arctan\!\left[\tfrac12(\tan\beta_{in}+\tan\beta_{out})\right]$$

**Holgura (AM)**: Y_cl = B·(t_cl/H)·Z, B=0 estátor, 0.25–0.50 rotor [V].
**Borde de fuga (AM)**: factor multiplicativo tabulado; t_te/s = [0, 0.02, 0.04, 0.08, 0.12] → y_te = [0.914, 1.000, 1.105, 1.375, 1.680] [V].
**Flare ≤ 12.5°** (recomendación AM, citada en AxialOpt) [V].

### 2.2 Dunham & Came (1970) — DOI 10.1115/1.3445349

Tres cambios verificados [V]:

$$f_{Re} = \left(\frac{Re}{2\times10^5}\right)^{-0.20}\qquad f_{Ma} = 1 + 60(\mathrm{Ma}_{out}-1)^2\ (\mathrm{Ma}_{out}>1)$$
$$Y_s = 0.0334\,\frac{c}{H}\,\frac{\cos\beta_{out}}{\cos\beta_{in}}\,Z\qquad \boxed{Y_{cl}=B\,\frac{c}{H}\left(\frac{t_{cl}}{H}\right)^{0.78} Z},\ B=0.37\ \text{(rotor con shroud)},\ 0\ \text{(estátor)}$$

(B≈0.47 sin shroud y t_cl→t_cl·N_seals^−0.42 para multi-labio: [?] — verificar en el original.) Y_p·f_Ma·f_Re; Y_s·f_Re; la holgura no se corrige por Re [V].

### 2.3 Kacker & Okapuu (1982) — la referencia de facto — DOI 10.1115/1.3227240

> **Precisión declarada [S]: ±1.5 % en η sobre 33 turbinas modernas de carga convencional** + carta de Smith. El criterio "sin recalibración por máquina" de Phy-AC.

**Factores de compresibilidad** (la aceleración del canal *reduce* la pérdida de perfil) [V]:

$$K_1=\begin{cases}1 & \mathrm{Ma}_{out}<0.2\\ 1-1.25(\mathrm{Ma}_{out}-0.2) & 0.2\le \mathrm{Ma}_{out}<1\end{cases},\qquad K_2=\left(\frac{\mathrm{Ma}_{in}}{\mathrm{Ma}_{out}}\right)^{2},\qquad K_p = \max(0.1,\ 1 - K_2(1-K_1))$$

**Pérdida de choque de raíz** [V]:

$$Y_{shock}=0.75\,[\max(0,\,f_{hub}\mathrm{Ma}_{in}-0.4)]^{1.75}\,r_{ht}\,\frac{p_{0,rel,in}-p_{in}}{p_{0,rel,out}-p_{out}}$$

| r_ht | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 |
|---|---|---|---|---|---|---|
| estátor f_hub | 1.40 | 1.18 | 1.05 | 1.00 | 1.00 | 1.00 |
| **rotor f_hub** | **2.15** | **1.70** | **1.35** | **1.12** | 1.00 | 1.00 |

**r_ht=0.5 es el suelo de la correlación** [V] → justifica HTR ≥ 0.5 como restricción.

**Perfil KO completo** [V]:

$$\boxed{\;Y_p^{KO}=f_{Re}\cdot f_{Ma}\cdot 0.914\left[\tfrac{2}{3}\,Y_p^{AM}\,K_p + Y_{shock}\right]}$$

**Reynolds KO (tres tramos)** [V]:

$$f_{Re}=\begin{cases}(Re/2\times10^5)^{-0.4} & Re<2\times10^5\\ 1 & 2\times10^5\le Re\le 10^6\\ (Re/10^6)^{-0.2} & Re>10^6\end{cases}$$

> **NO trasplantar el f_Re de Phy-AC a Phy-AT** (exponentes y tramos distintos; KO penaliza por encima de 10⁶).

**Secundaria KO** [V]:

$$Y_s = 1.2\,K_s\cdot 0.0334\cdot f_{AR}\cdot Z\cdot\frac{\cos\beta_{out}}{\cos\theta_{in}};\quad K_s = \max(0.1,\ 1-(b/H)^{2}(1-K_p));\quad f_{AR}=\begin{cases}\dfrac{1-0.25\sqrt{|2-H/c|}}{H/c} & H/c<2\\ c/H & H/c\ge2\end{cases}$$

**Borde de fuga KO** (Δφ² vs r_to=t_te/o, interpolación reacción/impulso con la regla del signo; Y_te = 1/(1−Δφ²) − 1; r_to ≤ 0.4) [V].
**Holgura KO** = forma de Dunham–Came [V].

**Correlaciones geométricas auxiliares de Kacker (1982)** — para sintetizar la geometría que los modelos exigen y θ no contiene [V]: carta de **stagger** ξ(β_in, β_out) (tabla digitalizada en AxialOpt) y **espesor máximo**:

$$\frac{t_{max}}{c}=\begin{cases}0.15 & \Delta\beta\le40°\\ 0.15+\tfrac{0.10}{80}(\Delta\beta-40) & 40°<\Delta\beta<120°\\ 0.25 & \Delta\beta\ge120°\end{cases}$$

### 2.4 Moustapha, Kacker & Tremblay (1990) — incidencia off-design — DOI 10.1115/1.2927647

$$\chi = \left(\frac{d_{le}}{s}\right)^{-1.6}\left(\frac{\cos\theta_{in}}{\cos\theta_{out}}\right)^{-2}(\beta_{in}-\beta_{des})\qquad |\chi|\le 800\ \text{(duro)}$$

$$\Delta\phi_p^2=\begin{cases}-5.1734\times10^{-6}\chi + 7.6902\times10^{-9}\chi^2 & -800\le\chi\le0\\ 0.778\times10^{-5}\chi+0.56\times10^{-7}\chi^2+0.4\times10^{-10}\chi^3+2.054\times10^{-19}\chi^6 & 0\le\chi\le800\end{cases}$$

Corrección secundaria por incidencia: χ_s con rango −0.4 < χ_s < 0.3 y Y_corr = e^{0.9χ_s}+13χ_s²+400χ_s⁴ (χ_s≥0) / e^{0.9χ_s} (χ_s<0) [V].

### 2.5 Benner, Sjolander & Moustapha (1997, 2004, 2006) — el estado del arte

- 1997: J. Turbomach. 119(2):193–200, DOI 10.1115/1.2841101 (incidencia con geometría de LE).
- 2004: 126(2):277–287, DOI 10.1115/1.1645533.
- 2006 Parte I: 128(2):273–280, DOI 10.1115/1.2162593 (**desglose + penetración**).
- 2006 Parte II: 128(2):281–291, DOI 10.1115/1.2162594 (**nueva secundaria**). [V todas]

**El nuevo esquema de descomposición** — dentro de la banda del vórtice de pasaje no hay pérdida de perfil 2D:

$$\boxed{\;Y_{tot} = \left(Y_p + Y_{te} + Y_{inc}\right)\left(1-\frac{Z_{TE}}{H}\right) + Y_s + Y_{cl}\;}$$

**Penetración** [V]:

$$\frac{Z_{TE}}{H} = \frac{0.10\,F_t^{0.79}}{\sqrt{CR}\;\mathrm{AR}^{0.55}} + 32.70\left(\frac{\delta^*}{H}\right)^{2};\qquad F_t = 2\,\frac{s}{c_{ax}}\cos^2\beta_m\,(|\tan\beta_{in}|+|\tan\beta_{out}|)$$

**F_t es el coeficiente de Zweifel en notación de Benner** [V] — usar Zweifel como variable de diseño alimenta la correlación sin conversión. CR = cosβ_in/cosβ_out; δ*/H escalado (Re_in/3×10⁵)^(−1/7).

**Secundaria (Parte II)** [V]:

$$Y_s=\begin{cases}\dfrac{0.038+0.41\tanh(1.20\,\delta^*/H)}{\sqrt{\cos\xi}\;CR\;\mathrm{AR}^{0.55}\left(\frac{\cos\beta_{out}}{\cos\xi}\right)^{0.55}} & \mathrm{AR}\le2\\[3ex] \dfrac{0.052+0.56\tanh(1.20\,\delta^*/H)}{\sqrt{\cos\xi}\;CR\;\mathrm{AR}\left(\frac{\cos\beta_{out}}{\cos\xi}\right)^{0.55}} & \mathrm{AR}>2\end{cases}$$

**Incidencia de Benner (1997)** [V]:

$$\chi = \left(\frac{d_{le}}{s}\right)^{-0.05} We^{-0.2}\left(\frac{\cos\theta_{in}}{\cos\theta_{out}}\right)^{-1.4}(\beta_{in}-\beta_{des});\qquad \Delta\phi_p^2=\begin{cases}b_1\chi+b_2\chi^2 & \chi<0\\ \sum_{i=1}^{8}a_i\chi^i & \chi\ge0\end{cases}$$

a = [−6.149e−5, +1.327e−3, −2.506e−4, −1.542e−4, +9.017e−5, +1.106e−5, −5.318e−6, +3.711e−7]; b₁=−8.720e−4, b₂=+1.358e−4 [V].

> **Salvaguardas de TurboFlow a copiar** [V]: extrapolación lineal desde χ̂=5 y limitación suave (`smooth_minimum` logsumexp, α=25) a 0.5. Motivo: **g finito y continuo para θ degenerados** — la dominancia de Deb lo exige.

Y_p, Y_te, Y_cl en Benner son **los de KO sin cambios** [V]: Benner **completa** a KO, no lo sustituye.

### 2.6 Craig & Cox (1970) — DOI 10.1243/PIME_PROC_1970_185_048_02

Estructura verificada en NASA TD3 [V]: η = (Trabajo − Grupo 2)/(Trabajo + Grupo 1) con X_p, X_s, X_a por fila y ~10 factores de figuras (N_pr Fig.3, N_pi Fig.10, N_pt Fig.6, δX_pm Fig.8, δX_p,s/e Fig.9, X_pb Fig.5 — cuya ordenada es X_pb(s/b)sinβ, error clásico de implementación documentado). **Descartado como base del L0**: ~18 figuras digitalizadas en `.pkl`, no diferenciable, lento. Comparación publicada (Lozza, Meccanica, DOI 10.1007/BF02128314) [S]: buen acuerdo con KO en subsónico de φ alto; **diferencias significativas en alta deflexión, sobre todo secundaria**.

### 2.7 Traupel — *Thermische Turbomaschinen* (Springer 2001)

ζ_pr = ζ_p·x_p·x_M·x_δ + ζ_δ + ζ_f con todo leído de Figuras 1–8 [V]. Chart-based → descartado como base. **Nota**: Coull & Hodson (2013) hallaron que su perfil + secundaria de Craig-Cox/Traupel daba los resultados más razonables en LPT high-lift [S] → mantener como **contraste offline**.

### 2.8 Denton (1993) — DOI 10.1115/1.2929299

Marco de **entropía**: fuentes = capas límite, mezcla, choques, y **transferencia de calor** (específica de turbinas refrigeradas, sin análogo en Phy-AC). Rol en Phy-AT: (a) contabilidad en entropía (como la fase 9 de Phy-AC); (b) reporte en ς = T_out·Δs/(½v_out²), la única magnitud aditiva a través de una máquina con T variable.

### 2.9 Coull & Hodson (2013) y Zhu & Sjolander (2005) — perfil moderno

- Coull & Hodson (2013), J. Turbomach. 135(2):021032 (+ 2012, 134(2):021002): carga del álabe y pérdida de perfil en LPT high-lift [S].
- Zhu & Sjolander (2005), GT2005-69077: correlación revisada de perfil + desviación (extiende Islam & Sjolander 1999) [S].

**Si el alcance incluye Zweifel > 1.1, AMDC-KO subestima la pérdida de perfil.** Recomendación pragmática: **empezar con KO/Benner y acotar Zweifel ≤ 1.15 por restricción**; ruta de ampliación documentada. Complementos: Kibsey & Sjolander (2016) GT2016-56410 (Mach en perfil); Yaras & Sjolander (1992) J. Turbomach. 114(1):204–210 (tip leakage físico, alternativa a DC/KO) [V bib].

### 2.10 Aungier (2006) — *Turbine Aerodynamics*, ASME Press

Doble contribución crítica: (1) **ajustes analíticos de las cartas AM** — sin ellos, KO no es implementable sin digitalizar cartas; (2) **modelo de desviación** C² (§4.4).

### 2.11 Comparación de precisión entre sistemas

| Modelo | Precisión reportada | Marca |
|---|---|---|
| Kacker–Okapuu 1982 | ±1.5 % en η, 33 turbinas + Smith | [S] |
| Craig–Cox vs KO | acuerdo en subsónico φ alto; difieren en alta deflexión (secundaria) | [S] |
| Coull–Hodson + CC/Traupel secundaria | "los más razonables" en LPT high-lift | [S] |
| Benner (TurboFlow) | caudal 100 % <2.5 %; par 78–94 % <5 %; ángulo 74–92 % <2.5 % en 3 turbinas NASA | [V] |

Comparativas adicionales: Dahlquist (2008, Lund); Jouybari et al. GT2012-69149; comparativa aire/sCO₂/ORC (S2666202722000210) [S].

### 2.12 Tabla consolidada de rangos de validez (para construir g) — todos [V]

| Cantidad | Límite | Origen |
|---|---|---|
| \|β_out\| | **40°–80°** | cartas AM/Aungier |
| s/c | **0.30–1.10** | sistema AM |
| t_te/o | **0–0.40** | sistema KO |
| r_ht | ≥ **0.50** | tabla f_hub KO |
| β_in estátor | ≤ +30° (AxialOpt: +15°) | KO |
| β_in rotor | ≥ −30° (AxialOpt: −15°) | KO |
| \|χ\| (Moustapha) | ≤ 800 | base experimental |
| χ_s | −0.4 – +0.3 | base experimental |
| Flare | ≤ 12.5° | AM 1951 |
| AR = H/c | ≥ 1.0; 3.0–4.0 "seguros" | Saravanamuttoo |
| R | 0.10–0.90 | AxialOpt |

---

## 3. Criterios de carga y restricciones

### 3.1 Zweifel

$$\psi_Z = 2\,\frac{s}{c_{ax}}\cos^2\alpha_2\,(\tan\alpha_1+\tan\alpha_2)$$

**Idéntico al F_t de Benner** [V]. Regla original: pérdida mínima en 0.8–1.0; rango práctico 0.8–1.1; moderno 0.75–1.15, con un óptimo citado de 0.829 [S]. De ψ_Z se despeja s/c_ax → paso → N_b = 2πr_m/s (entero, preferible coprimo con filas vecinas).

### 3.2 AN² y esfuerzo centrífugo

$$\sigma_{ct,\text{raíz}} = \frac{\rho_b\,\omega^2\,A_{ann}}{2\pi}\,K_{taper}$$

**Derivación auditable** (no cita): con σ_adm=250 MPa (Ni a ~1150 K), ρ_b=8200, K_taper=0.6 → AN²_max ≈ 4.5×10¹⁰ in²·rpm² — en la banda citada de 4–5×10¹⁰ para HPT. **Recomendación firme: usar la fórmula exacta de σ_ct como restricción; AN² solo como reporte** (las cifras de patentes usan unidades incoherentes — no anclar). Conflicto estructural↔η documentado: AN² bajo / ψ alto ⇒ η baja (Smith) — la frontera de Pareto de Phy-AT.

### 3.3 Límites de Mach

Estátor HPT: Ma₂ 0.95–1.2 normal (KO cubre transónico). Rotor: Ma_rel,3 ≲ 1.1. Salida última etapa: Ma ≲ 0.5–0.6 [S]. Difusor: Ma meridional ≤ 0.95 [V]. TurboFlow **no soporta cascadas supersónicas puras** (entrada Y salida supersónicas) [V] — límite a heredar.

### 3.4 Salida, giro, anillo

|β_out| ∈ [40°, 80°] [V]; swirl de salida |α| ≤ 10–20° [S]; flare ≤ 12.5° por defecto (AxialOpt usa ±10°; HPT modernas llegan a 25–35° [S] — exponer el límite como parámetro); deflexión hasta Δβ=120° (saturación de la correlación de espesor) [V].

### 3.5 Conjunto de restricciones de referencia (AxialOpt, anotado) [V]

| Restricción | min | max | Justificación |
|---|---|---|---|
| Flare | −10° | +10° | AM 1951 recomienda ≤12.5° |
| HTR | 0.600 | 0.900 | KO se sostiene hasta 0.5; HTR alto = secundaria grande |
| R | 0.100 | 0.900 | evitar compresión en rotor/estátor |
| Ma meridional difusor | — | 0.95 | no-tobera supersónica |
| β_in estátor | — | +15° | KO |
| β_in rotor | −15° | — | KO |
| h álabe | 0.01 m | — | opcional |
| t_te | 5e-4 m | — | opcional |

Variables con cotas [V]: ω_s, d_s ∈ [0.01,10] (heurística d_s·ω_s = 2/√N_et); velocidades reducidas; ángulo de salida ±[40°,80°]; AR ∈ [1.0,2.0]; s/c ∈ [0.30,1.10]; t_te/o ∈ [0.05,0.40]; holgura fija 0.5 mm ("0.2–0.5 mm razonables").

---

## 4. Off-design en turbinas

### 4.1 Ley de la elipse de Stodola

$$\frac{\dot m\sqrt{T_{01}}}{p_{01}} = K\sqrt{1-\left(\frac{p_{2}}{p_{01}}\right)^{2}}$$

Válida con toberas no estranguladas; satura en bloqueo; **cambia con pocas etapas** [S]. Refinamiento: Cooke (1985), DOI 10.1115/1.3239778 [V].

### 4.2 Choking — tres criterios implementables [V]

**(a) Garganta isentrópica**: Ma_throat = min(Ma_exit, 1) — sobreestima el caudal crítico.
**(b) Máximo caudal**: estado crítico con pérdidas; por encima, β_out = arccos(ṁ*/(ρ·w·A)).
**(c) Mach crítico analítico** (recomendado para L0 — sin iteración):

$$\mathrm{Ma}^*=\left(\frac{2}{\gamma-1}\right)^{1/2}\left[\frac{4\alpha-2}{(2\alpha+\phi^2-3)+\sqrt{(1+\phi^2)^2+4\alpha(1+\alpha-3\phi^2)}}-1\right]^{1/2},\ \alpha=\frac{\gamma}{\gamma-1}$$

**Primera tobera estrangulada ⇒ el mapa de turbina colapsa a casi una curva única** (a diferencia del compresor).

### 4.3 Incidencia off-design

Moustapha (1990) §2.4 y Benner (1997) §2.5. **Paralelismo con Phy-AC**: el diagnóstico de la fase 12/F-02 (sobre-pérdida a incidencia negativa) es el problema que Moustapha/Benner resuelven de forma citada — **Phy-AT nace con la ley correcta**.

### 4.4 Desviación

**Gauging**: β_g = arccos(o/s) = arccos(A_throat/A_out).
**AM** [V]: δ₀ = β_g − [35 + (45/39)(β_g−40)]; δ = δ₀·clip(1−X, 0, 1), X=(Ma−0.5)/(Ma*−0.5); **inexacto si β_g>70°, sin sentido >72°**.
**Aungier (2006)** [V]: δ₀ = arcsin[(o/s)(1+(1−o/s)(β_g/90)²)] − β_g, con p(X)=1−10X³+15X⁴−6X⁵ — **C² en ambos extremos**. **Usar Aungier**: la dominancia de Deb necesita continuidad.

### 4.5 Códigos de referencia off-design

AXOD (Glassman 1994, NTRS 19950004441); AXOD2 (Chen 2014); Chen 2009/2011; Flagg (1967) GE R66FPD258 (NTRS 19670009068); TD2-2 (Glassman 1992); Hendricks (2016) DOI 10.2514/6.2016-0119; Anderson et al. (2024) J. Turbomach.; Esfahanian et al. (2024) IJHFF 107:109370. [V bib]

---

## 5. Refrigeración a nivel meanline

### 5.1 Caudal de refrigerante — algoritmo de Gauntner (NASA TM-81453, 1980)

[NTRS 19800011581]. Implementación exacta (pyCycle) [V]:

$$\phi = \frac{T_{gas}-T_{metal}}{T_{gas}-T_{cool}};\qquad \phi' = \frac{\phi + PF}{1 + PF};\qquad \boxed{\dot m_{cool} = 0.022\,x_{factor}\,\tfrac{4}{3}\,\dot m_{prim}\left(\frac{\phi'}{1-\phi'}\right)^{1.25}}$$

- T_gas = T₀+T_safety (estátor); **0.92·T₀**+T_safety (rotor — ve la temperatura relativa).
- T_safety = 150 °R; T_metal = 2460 °R (≈1367 K) por defecto.
- **PF = 0.30 primera tobera** (perfil del combustor), 0.13 el resto.
- ṁ_cool = 0 si T_gas < T_metal — **discontinuidad a suavizar** para g continua.
- **x_factor = nivel tecnológico** (1 = actual; menor = más avanzada) [V] — el "nivel de refrigeración" pedido, en su forma más simple y calibrable.

Mezcla: h₀_out = (ṁ_prim·h₀_prim + ṁ_cool·h₀_cool)/(ṁ_prim+ṁ_cool) − Δh_rotor/ṁ_prim [V].

### 5.2 Pérdida de mezcla — Hartsel (AIAA 72-11, 1972) y sucesores

Volumen de control a **presión estática constante**; modificaciones de Köllen & Koschel (1985), Urban et al. (1998), Lim et al. (2010); formulación entrópica de Denton (1993) [S]. Mecanismos: pérdida en el orificio, mezcla, y cambio de la secundaria [S]. Revisión moderna: *Review of Efficiency Losses for a Cooled Turbine Stage*, J. Propulsion and Power, DOI 10.2514/1.B39445 [S].

### 5.3 Young & Wilcock (2002) — el marco termodinámico

Partes 1–2, J. Turbomach. 124(2):207–221. **"Para una evaluación precisa hay que dividir la expansión en etapas individuales tratando estátor y rotor por separado"** [S] — la justificación bibliográfica directa de la arquitectura fila a fila de Phy-AT.

### 5.4 La definición de eficiencia de una turbina refrigerada

Young & Horlock (2006), J. Turbomach. 128(4):658. **El E³ HPT dio 90.0 % (definición termodinámica) y 92.5 % (definición de ciclo GE) — 2.5 puntos de diferencia sobre la misma máquina** [S]. **Implicación: validar primero con turbinas de aire frío no refrigeradas (Kofskey); declarar la definición al abordar el E³.**

### 5.5 Niveles tecnológicos de refrigeración

| Nivel | Tecnología | T_metal aprox. | x_factor |
|---|---|---|---|
| 0 | Sin refrigeración | ~1100 K | ṁ_c=0 |
| 1 | Convección simple | ~1200 K | ~1.3 |
| 2 | Multipaso + turbuladores + impingement | ~1300 K | ~1.0 |
| 3 | + film cooling | ~1400 K | ~0.8 |
| 4 | + TBC | ~1450–1500 K ef. | ~0.6 |

(x_factor es la parametrización de pyCycle [V]; los valores por nivel son estimación [?] — calibrar contra los caudales publicados del E³.) Refs adicionales: Sammak et al. GT2013-95469; Genrup et al. GT2005-68716; Shahbazi et al. (2023) ATE 230:120828.

---

## 6. Casos de validación públicos con datos medidos

### 6.1 Tabla de casos

| # | Máquina | Referencia | Qué publica | Idoneidad | Marca |
|---|---|---|---|---|---|
| 1 | **NASA 2 etapas Ø8.00 in** (turbofán pequeño) | Kofskey & Nusbaum (1972), **NASA TN D-6967**, [NTRS 19720024422](https://ntrs.nasa.gov/citations/19720024422) | geometría completa + η de 1ª etapa Y de 2 etapas; 0–110 % velocidad; PR 1.79–5.14; par, caudal, ángulo de salida | ⭐⭐⭐⭐⭐ **El mejor caso.** Dos máquinas en una referencia, aire frío (η sin ambigüedad), geometría ya digitalizada (§6.3) | [S/V] |
| 2 | **NASA 1 etapa** (turborreactor bajo coste) | Kofskey, Roelke & Haas (1974), **TN D-7625**, [NTRS 19740018139](https://ntrs.nasa.gov/citations/19740018139) | diseño + aire frío | ⭐⭐⭐⭐⭐ **Estrangula en el ESTÁTOR** (el caso 1 estrangula en el ROTOR) — juntos cubren ambos modos | [V] |
| 3 | Variante con estátor abierto | Kofskey et al. (1974), E-7776, [NTRS 19740019165](https://ntrs.nasa.gov/citations/19740019165) | garganta modificada | ⭐⭐⭐⭐ Test de sensibilidad de capacidad a la garganta (regla de gauging) | [V] |
| 4 | Estudio aerodinámico de (1) | ASME 73-GT-29 | análisis | ⭐⭐⭐ | [V] |
| 5 | **E³ HPT (NASA/GE)** | Timko (1984), **NASA CR-168289** | 2 etapas, banco caliente escala completa **con refrigeración simulada**; η = 90.0 %/92.5 % según definición | ⭐⭐⭐⭐ Fase 2 — obliga a declarar la definición de η | [S] |
| 6 | Hardware del E³ HPT | [NTRS 19850002687](https://ntrs.nasa.gov/citations/19850002687) | geometría detallada | ⭐⭐⭐⭐⭐ complemento indispensable de (5) | [S] |
| 7 | E³ core | [NTRS 19900019243](https://ntrs.nasa.gov/api/citations/19900019243/downloads/19900019243.pdf) | ciclo, sangrados | ⭐⭐⭐ | [S] |
| 8 | **E³ LPT** | CR no confirmado | LPT multietapa | ⭐⭐⭐⭐ **Localizar el número en NTRS antes de comprometerlo** | [?] |
| 9 | Turbina alta temperatura, etapas 3–4 | NTRS 19680006274, 19720024133 | aire frío | ⭐⭐⭐⭐ | [S] |
| 10 | Turbina de potencia libre | NTRS 19790009688 | aire frío | ⭐⭐⭐ velocidad variable | [S] |
| 12 | **LISA (ETH Zürich)** | Behr, Sell et al.; [LEC/ETH](https://lec.ethz.ch/research/turbomachinery_experimental.html) | 2 etapas, potencia POR etapa, shroud/purge cuantificados | ⭐⭐⭐ para Y_cl; geometría dispersa (curación cara) | [S] |
| 13 | Cascadas VKI / T106 / Pak-B | — | pérdida de perfil en cascada | ⭐⭐⭐ nivel FILA (verificar Y_p aislado) | [?] |
| 14 | Aachen 1.5 | — | — | no verificada aquí — no comprometer | [?] |

### 6.2 Métricas de referencia — el listón a batir (TurboFlow: Benner + Aungier + máximo caudal, sin ajustes por máquina) [V]

| Caso | Magnitud | <2.5 % | <5.0 % | <10 % |
|---|---|---|---|---|
| 1 etapa, **estátor estrangulado** | caudal | **100.0 %** | 100.0 % | 100.0 % |
| | par | 25.6 % | 56.4 % | 82.0 % |
| | ángulo de salida | 87.8 % | 97.6 % | 100.0 % |
| 1 etapa, **rotor estrangulado** | caudal | **100.0 %** | 100.0 % | 100.0 % |
| | par | 77.1 % | 93.8 % | 100.0 % |
| | ángulo | 74.4 % | 94.9 % | 100.0 % |
| **2 etapas** | caudal | **100.0 %** | 100.0 % | 100.0 % |
| | par | 69.9 % | 100.0 % | 100.0 % |
| | ángulo | 92.3 % | 100.0 % | 100.0 % |

(Barrido 70/90/100/110 % de ω; PR total-a-estático variable.)

**Lecturas**: (1) el **caudal se predice esencialmente perfecto** — si Phy-AT no lo logra, hay un bug; (2) el **par es la magnitud difícil**, sobre todo a 70 % — el modelo de incidencia es el eslabón débil, como en Phy-AC; (3) el ángulo valida la desviación de Aungier.

### 6.3 Geometrías ya digitalizadas [V]

Los YAML de TurboFlow contienen la geometría fila a fila (radios hub/tip in/out, pitch, cuerda, stagger, opening, ángulo y diámetro de LE, ángulo de cuña, t_te, t_max, holgura, throat_location_fraction) para **Kofskey1972 1-etapa, Kofskey1972 2-etapas y Kofskey1974**, con `experimental_data/` y `quantify_error.py`. Licencia MIT. **Los tres declaran holgura ABSOLUTA en metros** — la elección que Phy-AC adoptó tras descubrir que ε/h fija borraba el crecimiento trasero. **Phy-AT debe nacer con holgura absoluta.**

---

## 7. Libros y referencias canónicas

| Obra | Uso |
|---|---|
| **Glassman (ed.), *Turbine Design and Application*, NASA SP-290, Vols. 1–3** ([NTRS 19950015924](https://ntrs.nasa.gov/citations/19950015924)) | El manual meanline de turbinas de NASA; base del modelo "TD2" de TD3 [V] |
| **Dixon & Hall** (7ª ed., 2014) | Triángulos, etapa repetitiva, Smith — coherencia de convenciones con Phy-AC garantizada |
| **Saravanamuttoo et al., *Gas Turbine Theory*** (2008) | límites de AR y HTR citados en AxialOpt [V] |
| **Moustapha, Zelesky, Baines & Japikse, *Axial and Radial Turbines*** (Concepts NREC, 2003) | la referencia moderna de pérdidas y off-design — sus autores son los de las correlaciones §2.4–2.5 |
| **Japikse & Baines** | diseño preliminar |
| **Sieverding (VKI)** | física de secundaria y TE |
| **Aungier, *Turbine Aerodynamics*** (2006) | ajustes AM + desviación — imprescindible |
| **Traupel** (2001) | contraste europeo |
| **Denton (1993)** IGTI Scholar Lecture | marco de entropía |
| **Wilson (1987)**, Proc. IMechE A 201(4):279–290 | guías de diseño preliminar |
| **Balje & Binsley (1968)**, J. Eng. Power 90(4) | pérdida-geometría y optimización — precursor de Phy-AT |
| **Pritchard (1985)**, ASME 85-GT-219 | meanline → perfil (capa 5) |
| **Denton (2017)** MULTALL, J. Turbomach. 139(12):121001 | vía de verificación L1/L2 abierta |

---

## 8. Recomendaciones concretas para Phy-AT L0

### 8.1 Cadena de pérdidas — decisión

**`AMDC-KO + Benner`, como el modelo `benner` de TurboFlow:**

$$\boxed{\;Y_{tot} = \bigl(Y_p^{KO} + Y_{te}^{KO} + Y_{inc}^{B97}\bigr)\left(1-\frac{Z_{TE}}{H}\right) + Y_s^{B06} + Y_{cl}^{DC}\;}$$

con desviación de **Aungier** (C²) y choking por **Mach crítico analítico**.

**Seis razones**: (1) cerrada en forma algebraica — 0.5 ms/punto alcanzable en NumPy; (2) ±1.5 % publicado sin recalibración; (3) validación abierta y reproducible contra 3 turbinas NASA; (4) **Benner aporta la incidencia gratis** — el talón de Aquiles de Phy-AC (fase 12/F-02) resuelto de fábrica; (5) la penetración de Benner sustituye el ajuste global de pared que Phy-AC tuvo que recalibrar (K_ENDWALL); (6) Craig–Cox/Traupel descartados por evidencia directa (~18 cartas en `.pkl`), mantenidos como contraste offline (`--loss-model craig_cox`).

**Auditoría en entropía** (Denton 1993): marcha en s, reporte en ς. **Reserva LPT high-lift**: Zweifel ≤ 1.15 por restricción; ruta a Coull–Hodson / Zhu–Sjolander documentada.

### 8.2 Vector de diseño θ propuesto (18 dimensiones)

| # | Variable | Rango | Notas |
|---|---|---|---|
| 0 | `n_stages` | 1–6 | entero dentro de evaluate |
| 1 | `RPM` | 3k–30k | |
| 2 | `HTR_in` | **0.50–0.92** | suelo de la tabla f_hub de KO [V] |
| 3 | `phi_mid` | 0.40–1.00 | |
| 4 | `psi_mid` | **0.80–2.60** | ψ=Δh₀/U² — declarar convención |
| 5 | `psi_slope` | −0.35–0.35 | |
| 6 | `R_mean` | 0.15–0.60 | |
| 7 | `R_slope` | −0.25–0.25 | R crece hacia atrás |
| 8 | `phi_slope` | −0.30–0.30 | |
| 9 | `Zweifel_stator` | **0.70–1.15** | → s/c_ax → N_b; alimenta F_t de Benner |
| 10 | `Zweifel_rotor` | 0.70–1.15 | |
| 11 | `AR_stator` | 1.0–5.0 | |
| 12 | `AR_rotor` | 1.0–5.0 | |
| 13 | `t_te_over_o` | 0.02–0.30 | duro KO ≤0.40 |
| 14 | `alpha_in` | −10°–+10° | |
| 15 | `cool_tech` (x_factor) | 0.5–1.5 | inactivo si T₀<T_metal |
| 16–18 | `T0_in`, `P0_in`, `massflow` | — | pinned por el spec, al final (como Phy-AC) |

**Decisiones**: Zweifel en lugar de solidez libre (variable física, rango óptimo estrecho, = F_t de Benner); la geometría detallada NO va en θ — se **sintetiza** con las correlaciones de Kacker (stagger, t_max) + defaults (d_le/c≈0.03–0.06, We=50° como en los casos NASA [V]); **el radio de punta NO es variable** (lección M1 de Phy-AC); **holgura absoluta en mm**; objetivo min(−η_tt, |PR−PR*|/PR*) s.a. g≤0.

### 8.3 Restricciones g(θ)≤0 propuestas (continuas y finitas para θ degenerados)

**Grupo A — validez de correlaciones**: A1 40°≤|β_out|≤80°; A2 0.30≤s/c≤1.10; A3 t_te/o≤0.40; A4 HTR≥0.50; A5 β_in estátor ≤+30°, rotor ≥−30°; A6 χ de Benner en banda extrapolable (smooth_min); A7 β_g≤70° si desviación AM.

**Grupo B — aerodinámicas**: B1 Ma_rel salida rotor ≤1.15; B2 Ma salida ≤0.60; B3 |swirl salida|≤15°; B4 0.05≤R_i≤0.85; B5 0.70≤ψ_Z≤1.15; B6 |flare|≤12.5°; B7 1.0≤AR≤5.0; **B8 margen de choking** (caudal crítico por fila ≥ diseño × margen; o la fila estranguladora es la deseada) — el análogo estructural del margen de bombeo de Phy-AC.

**Grupo C — mecánicas/refrigeración**: C1 σ_ct raíz exacta ≤ σ_adm(T_metal) — **no AN² tabulado**; C2 T_metal ≤ T_max(tecnología); C3 Σṁ_cool/ṁ ≤ 0.20–0.25; C4 h ≥ 10 mm, N_b entero ≥ mínimo; C5 U_tip, RPM, r_tip del spec.

**Grupo D — spec**: |PR−PR*|, potencia, n_stages.

### 8.4 Campaña de validación — tres fases

**Fase 1 — aire frío NASA** (sin ambigüedad de η, geometría digitalizada, mapas medidos):

| Máquina | Plano | Tolerancia | Base |
|---|---|---|---|
| Kofskey 1974 (D-7625), estátor estrangulado | etapa | ṁ ±2.5 % · par ±5 % · α ±5° | 100/56/98 % en TurboFlow |
| Kofskey 1972 (D-6967), 1ª etapa, rotor estrangulado | etapa | ṁ ±2.5 % · par ±5 % · α ±5° | 100/94/95 % |
| Kofskey 1972, 2 etapas | máquina | ṁ ±2.5 % · par ±5 % · α ±2.5° | 100/100/92 % |
| Kofskey 1974 (E-7776), estátor abierto | etapa | ṁ ±3 % | sensibilidad de garganta |

Mapas: 70/90/100/110 % de ω, PR 1.6–4.5. **Métrica clave adicional: caudal estrangulado** (lo que Phy-AC falla en +6.6 % en compresor — aquí debe salir a 2.5 %, la garganta lo gobierna todo).

**Fase 2 — E³ HPT** (CR-168289 + NTRS 19850002687): añade refrigeración; **declarar la definición de η**; tolerancia ±2 pts sobre la termodinámica.

**Fase 3 — ampliación**: NTRS 19680006274/19720024133, 19790009688, E³ LPT (localizar CR), LISA (Y_cl con/sin shroud, purge).

**Regla invariante**: geometría y punto de operación publicados; **cero recalibración por máquina**; constantes globales congeladas antes de la campaña.

### 8.5 Riesgos y mitigación

| Riesgo | Mitigación |
|---|---|
| Convención de signos | función única ax↔tan + test de identidad 1e-6 + anclas desde el día 1 |
| Y referenciado a la dinámica de SALIDA | test Y_definición vs Y_modelo por fila |
| f_Re de compresor ≠ turbina | usar el de KO |
| cp(T) a 1900 K + FAR | extender JANAF; test de rango |
| ṁ variable por refrigeración | refactor desde el inicio |
| Discontinuidades (refrigeración on/off, choke, saturaciones) | smooth_min/logsumexp + interpolante quíntico de Aungier |
| Geometría que θ no contiene | correlaciones de Kacker + defaults calibrados sobre los 3 casos NASA |
| Definición de η refrigerada (2.5 pts) | validar primero en aire frío |
| Zweifel alto fuera de KO | restricción ≤1.15; ruta Coull–Hodson/Zhu–Sjolander |

### 8.6 Coste computacional

Los únicos bucles por punto: continuidad/Mach por estación (bisección ~15 iter.); desviación **cerrada**; Mach crítico **analítico**; marcha con cp(T) polinómico. **Menos iteración que Phy-AC** (no hay bucle de bombeo). **~0.5 ms/punto en NumPy es realista.**

---

## 9. Bibliografía consolidada

**Modelos de pérdidas**: Ainley & Mathieson (1951) ARC R&M 2974 y 2891 ([Cranfield](https://reports.aerade.cranfield.ac.uk/handle/1826.2/3538)); Dunham & Came (1970) DOI 10.1115/1.3445349; Craig & Cox (1970) DOI 10.1243/PIME_PROC_1970_185_048_02; Kacker & Okapuu (1982) DOI 10.1115/1.3227240; Moustapha, Kacker & Tremblay (1990) DOI 10.1115/1.2927647; Tremblay et al. (1990) 90-GT-314; Benner et al. (1997) DOI 10.1115/1.2841101; (2004) DOI 10.1115/1.1645533; (2006 I/II) DOI 10.1115/1.2162593 / 10.1115/1.2162594; Zhu & Sjolander (2005) GT2005-69077; Coull & Hodson (2012/2013) J. Turbomach. 134(2):021002 / 135(2):021032; Yaras & Sjolander (1992) DOI 10.1115/1.2927987; Denton (1993) DOI 10.1115/1.2929299; Traupel (2001) DOI 10.1007/978-3-642-17469-8; Aungier (2006) ASME Press; Lozza DOI 10.1007/BF02128314; Dahlquist (2008, Lund); Kibsey & Sjolander (2016) GT2016-56410.

**Diseño y eficiencia**: Smith (1965) JRAeS 69:467; Balje & Binsley (1968); Wilson (1987) DOI 10.1243/PIME_PROC_1987_201_035_02; Pritchard (1985) 85-GT-219; Agromayor & Nord (2019) DOI 10.3390/ijtpp4030032; Agromayor, Müller & Nord (2019) DOI 10.3390/ijtpp4030031; Anderson et al. (2025) DOI 10.21105/joss.07588; Esfahanian et al. (2024) DOI 10.1016/j.ijheatfluidflow.2024.109370.

**Off-design y choking**: Cooke (1985) DOI 10.1115/1.3239778; Flagg (1967) NTRS 19670009068; Glassman (1994) NTRS 19950004441; Chen (2014); Hendricks (2016) DOI 10.2514/6.2016-0119; Graham & Kost (1979) 79-GT-37.

**Refrigeración**: Hartsel (1972) DOI 10.2514/6.1972-11; Gauntner (1980) NTRS 19800011581; Young & Wilcock (2002) J. Turbomach. 124(2):207–213 y 214–221; Young & Horlock (2006) 128(4):658; Sammak et al. (2013) GT2013-95469; Shahbazi et al. (2023) DOI 10.1016/j.applthermaleng.2023.120828.

**Validación**: Kofskey & Nusbaum (1972) NTRS 19720024422; Kofskey, Roelke & Haas (1974) NTRS 19740018139; Kofskey, Nusbaum & Haas (1974) NTRS 19740019165; ASME 73-GT-29; Timko (1984) NASA CR-168289; NTRS 19850002687; NTRS 19900019243.

**Libros**: Glassman SP-290 (NTRS 19950015924); Dixon & Hall (2014); Saravanamuttoo et al. (2008); Moustapha et al. (2003); Denton (2017) DOI 10.1115/1.4037819.

**Software de referencia (código verificado)**: [TurboFlow](https://github.com/turbo-sim/turboflow) (MIT); [AxialOpt](https://github.com/turbo-sim/AxialOpt) (MIT); [NASA turbo-design](https://github.com/nasa/turbo-design); [pyCycle](https://github.com/OpenMDAO/pyCycle).

---

**Nota final de honestidad**: las ecuaciones de §2.1–2.5, 3.5, 4.2, 4.4, las geometrías §6.3 y las métricas §6.2 están verificadas línea a línea sobre código abierto y son transcribibles. Los números **[S]** (±1.5 % de KO, rangos de Zweifel, η del E³, contenido de SP-290) tienen cita correcta pero requieren confirmación sobre el PDF. Los **[?]** (E³ LPT CR, B=0.47 sin shroud, x_factor por nivel, rangos φ/ψ por tipo) **no deben usarse como ancla sin verificación**.
