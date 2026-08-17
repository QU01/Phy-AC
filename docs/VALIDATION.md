# Phy-AC — Campaña de validación

Distinción VV&UQ del proyecto: `test_phyac.py` VERIFICA (¿resolvemos bien
las ecuaciones? — 165 checks); `validation/validate.py` VALIDA (¿las
ecuaciones correctas? — contra máquinas NASA medidas). CI corre ambos en
estricto.

## 1. Metodología

Para cada máquina publicada, `validate.py` construye el θ que reproduce
su **annulus** (r_tip desde U_tip/ω publicados; φ1 invertido por
bisección para que la continuidad devuelva ese r_tip) y su **trabajo
medido** (Δh0 implicado por PR y η publicados — derivado con el MISMO
gas caloríficamente imperfecto que usa el meanline desde la fase 9; con
un cp de referencia se le pediría a la máquina una ψ inconsistente con su
propio calor específico). El modelo
recibe el trabajo real y se califica su predicción de PÉRDIDAS → (η, PR).
Sin recalibración por máquina. Desde la fase 8 las multietapa pueden
declarar su distribución por etapa con el campo `slopes`
(phi_slope/Rx_slope del θ 15-D); las monoetapa usan el θ legacy de 13
(paddeado con pendientes 0, bit-exacto).

Casos `rotor` (Rotor 37/67): se califica el PR/η del rotor derivado del
desglose de pérdidas de la etapa (sin estátor).

## 2. Resultados vigentes (2026-08-16, fase 9)

| Máquina | Plano | ΔPR | Δη | Tolerancia |
|---|---|---|---|---|
| NASA Stage 35 | etapa | +1.2% | +1.7 pts | 5% / 2 pts |
| NASA Rotor 37 | rotor | −1.2% | −1.6 pts | 5% / 3 pts |
| NASA Rotor 67 | rotor | −1.2% | −2.5 pts | 5% / 3 pts |
| GE/NASA E³ HPC (10 et.) | máquina | −5.6% | −1.6 pts | 6% / 3 pts |

Las cuatro pasan con el modelo de la fase 9, que sustituye tres ajustes
por física citada (§3). El peor uso de tolerancia es del 93%.

Tabla viva en `validation/RESULTS.md` (regenerar tras tocar
`physics_core.py`).

## 3. Calibraciones ancladas

- **K_SHOCK = 0.70**: el modelo de choque normal a M de entrada
  sobreestima el choque oblicuo real del pasaje; calibrado contra Rotor
  37/67 manteniendo Stage 35 en tolerancia. Documentado en
  `physics_core.py`.
- Conversión arrastre→pérdida del endwall: ω̄ = C_D·σ·cos²β₁/cos³β_m
  (Dixon §3). El bug inicial (factor invertido) dominaba el error de η
  (−17 pts) — corregido y cubierto por las anclas de regresión.
- **Holgura de punta por fila + K_ENDWALL = 1.4 (2026-07-17)**: ε pasa a
  ser ABSOLUTA en mm (ε/h crece hacia las etapas traseras) con los
  regímenes de Sakulkaew 2013, y la validación inyecta la ε publicada de
  cada máquina (R37 0.356 mm AGARD; R67 ≈1.0 mm; Stage 35 0.36 mm y E³
  0.5 mm APROXIMADOS). Al quitar el débito uniforme viejo (ε/h=1.5% ≈ 3
  pts) las máquinas quedaban altas: recalibración GLOBAL documentada de
  K_ENDWALL 1.0→1.4 (el CDa de Howell subestima el endwall — Koch &
  Smith 1976 — y el débito uniforme lo absorbía en silencio). Resultado:
  Δη máx |1.6| pts (antes |2.5|); ΔPR del E³ −4.4%→−5.5% (dentro de su
  tolerancia relajada — es el límite de parametrización por etapa, fase
  8). Ancla REF_AX4 re-congelada: PR 2.7056→2.7160, η 0.8837→0.8871.
- **Corrección de Reynolds (2026-07-16)**: f_Re multiplica perfil y
  endwall (Koch & Smith 1976 nominal a Re_c=10⁶; Re^−0.2 turbulento,
  Re^−0.5 bajo 2×10⁵ — Wassell 1968 / Schäffler 1980). SIN crédito por
  encima de 10⁶: las máquinas NASA de esta tabla corren a Re_c ≳ 10⁶ y
  sus deltas no se movieron (Stage 35 +0.89→+0.87% PR; resto idéntico).
  El ancla REF_AX4 sí se movió (filas traseras con Re < 10⁶:
  η 0.8900→0.8837) y se re-congeló citando esta corrección.

## 4. Anclas de regresión

`REGRESSION_ANCHORS` en `machines.py` congela la salida del meanline para
el θ de referencia (REF_AX4). NO son mediciones: detectan deriva
silenciosa de la física. Actualizarlas (`--freeze-anchors`) es una
decisión consciente que debe citar la corrección que las movió.

- **Pendientes por etapa del E³ (2026-07-17, fase 8)**: la entrada E³
  declara `slopes=dict(phi_slope=-0.10, Rx_slope=0.10)` — dirección de
  la práctica real de HPCs (φ cae hacia atrás, Rx crece), magnitud
  APROXIMADA (fit de 2 parámetros; el CR-165558 no está accesible).
  Recuperan ~0.75 pts del déficit de PR (−5.55% → −4.80%) y la
  tolerancia se endurece 8% → 6%. El déficit restante NO es de
  parametrización (candidatos: cp constante a τ≈2.4, acumulación de
  bloqueo, WDF). Las anclas NO se movieron (padding bit-exacto).

### Fase 9 (2026-08-16) — reconstrucción de la contabilidad de pérdidas

Cuatro cambios de física simultáneos, cada uno con su fuente. Anclas
re-congeladas citando ESTA entrada.

1. **Gas caloríficamente imperfecto.** cp(T) y γ(T) del aire (ajuste
   cuadrático de JANAF, <1% en 250-1000 K) en todo el stacking. Los
   rendimientos pasan a las formas EXACTAS para gas imperfecto vía la
   función de entropía phi(T) = ∫cp/T dT:
   η_poly = R·lnPR/Δphi y η_isen = [h(T2s)−h(T1)]/[h(T2)−h(T1)].
   A τ ≈ 2.4 (E³) cp sube ~7% y γ cae a 1.36: con cp constante el error
   entraba directo en PR.

2. **Pérdida → rendimiento por ENTROPÍA** (Dixon & Hall §5.5, ec.
   5.4-5.9). Antes: Δh_pérdida = ω̄·½W₁², que es ΔP₀/ρ a la entrada de
   la fila. Ahora la etapa se marcha con las presiones totales reales
   (relativa en el rotor, absoluta en el estátor), Δs sale de esa marcha
   y el equivalente en trabajo es T₀₃·Δs. El PR de etapa deja de venir
   de una fórmula y sale de la marcha. Además ω̄ se refiere a la cabeza
   dinámica COMPRESIBLE (P₀−p), que es la definición de Koch & Smith
   1976: usar ½ρW² subestimaba la pérdida un 60% a M_rel ≈ 1.4.

3. **Capacidad de subida de presión al stall: Koch 1981.**
   `CH_STALL_MAX = 0.48` (constante) → correlación real: Ch_ef,stall en
   función del parámetro de difusión L/g₂ de la cascada, corregida por
   Reynolds (Fig. 4), holgura de punta (Fig. 5) y espaciado axial
   (Fig. 6), y comparada contra el Ch dividido por el factor de cabeza
   dinámica efectiva 𝔉_ef (Fig. 13), que es lo que explica el fuerte
   efecto del tipo de triángulo de velocidades. El numerador de Ch pasa
   a ser la subida de presión estática ISENTRÓPICA (definición de Koch),
   no el trabajo real: usar Δh₀ sobreestimaba Ch en 1/η ≈ 1.12.

4. **Endwall y bloqueo: Koch & Smith 1976.** Desaparece el término de
   annulus de Howell (C_Da = 0.020·s/h) y con él el multiplicador
   `K_ENDWALL = 1.4` que lo tapaba. En su lugar, la suma de espesores de
   desplazamiento de pared 2δ*/g de su ec. (3) y Fig. 8 — que crece como
   x³ con la carga relativa x = Ch/Ch_stall y linealmente con ε/g — da
   (a) el BLOQUEO del annulus, que así emerge de la carga en vez de ser
   la recta inventada 0.98 − 0.005·i, y (b) el débito de rendimiento de
   su ec. (2), η = η̃·(1−Σδ*/h)/(1−Σν/h) con Σν ≈ 0.48·Σδ* (Fig. 10).
   Añadido el efecto del Mach de entrada sobre la pérdida de PERFIL
   (su Fig. 6: ω̄ casi se duplica entre M 0.1 y 1.5), que la correlación
   incompresible de Lieblein no conoce.

**Recalibración global documentada.** Con la física anterior había cinco
constantes de ajuste (K_PROFILE, K_ENDWALL, K_SHOCK y las dos del
bloqueo lineal). Ahora quedan TRES, y dos de ellas con respaldo directo:

| Constante | Valor | Origen |
|---|---|---|
| `K_SHOCK` | 1.00 (antes 0.70) | Koch & Smith 1976 Fig. 7: los coeficientes de pérdida de choque medidos en etapas de fan de alta velocidad caen SOBRE la curva de choque normal a M₁, no en el límite oblicuo |
| `K_ENDWALL` | 1.00 (antes 1.40) | vuelve a su valor natural: ya solo multiplica la secundaria de Howell (C_Ds = 0.018·C_L²) |
| `K_MACH_PROFILE` | 0.50 | ajuste de Koch & Smith Fig. 6, calibrado dentro de su incertidumbre de lectura |
| `EW_A` (2δ*/g) | 0.16 | ajuste de la rama ε/g = 0 de su Fig. 8, ídem |
| `K_PROFILE` | 1.24 (antes 1.00) | multiplicador global: el θ/c de Lieblein es de cascada 2-D lisa a incidencia de diseño; la máquina real pierde más |

Efecto medido: η del E³ pasa de −1.4 a −1.6 pts pero su PR de −4.8% a
−5.6%; los transónicos mejoran su PR (Stage 35 +0.8→+1.2%, R37
−1.1→−1.2%, R67 −0.7→−1.2%) a cambio de algo de η. El resultado neto es
un modelo con **menos ajuste libre y más mecanismo**, con las cuatro
máquinas dentro de tolerancia sin recalibración individual.

**Ancla REF_AX4 re-congelada** citando esta entrada:
PR 2.7160→2.7616, η_poly 0.8871→0.9013.

## 5. Pendientes

- **Validar el FUERA DE DISEÑO.** La fase 9 hizo el mapa físico (choke
  que limita gasto, línea de trabajo de Dixon ec. 5.26b, bombeo por el
  criterio de Koch o por pendiente nula, margen en gasto) pero NINGUNA
  de esas curvas está calificada contra medida. Es el hueco nº 1 de la
  campaña. Fuente disponible: AGARD AR-355 (Rotor 37) — digitalizar sus
  speedlines y añadir un criterio `kind="speedline"` a machines.py.
- Digitalizar speedlines completas de Rotor 37 (choke mdot ±2%) cuando
  haya fuente estable (las URL de turbmodels rotaron — 2026-07).
- Verificar la entrada E³ contra el CR original (distribución por etapa
  real → sustituir el fit de `slopes` y endurecer hacia 5%; entrada
  marcada APROXIMADA en machines.py).
- Con validación off-design (speedlines): calibrar el mapa VSV.
- Pares CFX/banco para `HiFiCalibration` (L2) del sesgo L1 (≈0.94 en PR).

### Paridad STL ↔ STEP (fase 10 · G-01)

`validation/parity_stl_step.py` compara las dos rutas de salida: la capa
5c (C#/PicoGK, ruta de fabricación) y la vía CadQuery (ruta de re-CAD).
Medido el 2026-08-16 sobre una máquina de 2 etapas (PR 1.29, r_punta
109 mm) a vóxel 0.6 mm, con el filete de raíz apagado en las dos rutas:

| Conjunto | STL | STEP | Δ | r exterior STL / STEP |
|---|---|---|---|---|
| Eje | 246.5 cm³ | 246.6 cm³ | −0.1 % | 20.00 / 20.00 mm |
| Rotor | 638.2 cm³ | 645.5 cm³ | −1.1 % | 108.66 / 109.05 mm |
| Carcasa | 950.6 cm³ | 957.0 cm³ | −0.7 % | 126.70 / 126.70 mm |

Tolerancias declaradas a partir de esa medida: 5 % en volumen por
conjunto y 1.5 mm en radio exterior. La comparación es por CONJUNTOS y no
pieza a pieza a propósito: las dos rutas parten la máquina en planos
distintos —el STL donde van las juntas apernadas, el STEP donde están las
piezas de re-CAD— y comparar pieza a pieza compararía esa decisión, no la
máquina.

La parte barata de la comparación (perfil de abeto punto a punto,
longitud y pendiente de plataforma por fila, tirantes, rebaje de rim) es
matemática pura y corre en cada pasada de la suite (bloque T20) y en el
CI; la cara (volúmenes, con PicoGK) corre en local.

**Hallazgo abierto**: el paso de filete de raíz de la capa 5c
(`voxWithRootFillets`, `Fillet` = over-offset morfológico con r = 2 mm)
añade **117 cm³ a un anillo de carcasa de 488 cm³ — un 24 %**. Un filete
de 2 mm sobre las raíces de ~70 vanos no puede pasar de unos pocos cm³.
El paso está añadiendo mucho más material del que declara, y ese material
NO está en el STEP ni en el margen estructural. Por eso la comparación se
hace con el filete apagado y por eso queda anotado aquí: pendiente de
investigar.

### Fase 11 (2026-08-16) — L1 pasa a ser un peldaño de verdad

`scm_core.py` sustituye a `turbo-design`: through-flow por curvatura de
líneas de corriente, equilibrio radial COMPLETO sobre 9 líneas,
continuidad por tubo de corriente, cierre por ángulos metálicos del álabe
y pérdidas resueltas en el span. Sin dependencia externa y en proceso
(~3 s por máquina).

**Verificación** (bloque T21 de la suite): con curvatura nula y
Cu ∝ r^n, el ODE integrado numéricamente reproduce la forma cerrada
`vortex_cx` con desvío 2e-16 en vórtice libre y ~2e-3 (error del trapecio
sobre 9 puntos) en n = −0.5, 0, +0.5. Y sobre el θ de referencia:

| Ley de torbellino | ΔPR vs L0 | Δη_poly | Reparto de trabajo en el span |
|---|---|---|---|
| libre (n = −1) | +0.13% | +0.34 pts | 6.9% |
| controlado (n = −0.5) | −6.5% | −0.01 pts | 5.6% |

La fila de vórtice libre es VERIFICACIÓN: es el caso donde las hipótesis
del meanline son exactas, así que coincidir es el resultado esperado. La
de torbellino controlado es el RESIDUAL que la capa 2 necesita.

**Banco de pruebas** (`validation/bench_scm.py` → `BENCH_SCM.md`, 80
diseños factibles por LHS sobre el espacio completo):

| | vórtice libre | torbellino controlado (n = −0.5) |
|---|---|---|
| cobertura | 85% | 75% |
| coste (mediana) | 3.8 s/máquina, ≈0.85 s por etapa | ídem |
| ΔPR mediana vs L0 | −0.08% (p10–p90: −2.2% … +1.3%) | −1.50% (−6.7% … −0.1%) |
| Δη mediana | +0.02 pts | −0.24 pts |
| reparto de trabajo en el span | 6.1% | 6.8% |
| convergencia | 15 iteraciones (p90 18) | ídem |
| independencia de malla (5→13 líneas) | dispersión mediana 0.33%, máx 1.97% | ídem |

**Hallazgo central del banco**: la cobertura CAE con el número de etapas
(100% en 1 etapa, 71-75% en 6-8). La razón es estructural, no numérica —
el annulus lo dimensiona L0 con su Cx uniforme y su densidad media, L1
resuelve un perfil, y el álabe de ángulo fijo convierte esa diferencia en
trabajo, que cambia la densidad, que cambia la siguiente estación. En 1-4
etapas es ruido; en 7-8 se compone hasta ±27% de PR. La cura de fondo es
que el annulus salga del MISMO solver que lo usa (o que L0 lo dimensione
con el perfil de L1); mientras tanto la guarda `PR_WINDOW` (±15%) rechaza
el punto y lo degrada a L0 etiquetado en vez de devolver un número que
nadie debería usar.

**Sigue SIN calificar contra medida.** L1 aporta variación de fidelidad y
distribución radial; no exactitud demostrada. Es el mismo hueco que el
mapa fuera de diseño (F-02) y, con speedlines medidas, la primera
comparación que hay que hacer.

**Hallazgo del solver**: aplicar el exponente de torbellino por igual a
Cu₁ y Cu₂ hace que el trabajo de Euler varíe como r^(n+1). Para n ≠ −1 el
diseño es de trabajo NO uniforme y el meanline solo lo evaluaba en la
línea media. En máquinas de varias etapas ese sesgo se compone y el
limitador de perfil del SCM llega a saltar (el punto degrada a L0
etiquetado). Es una limitación del MODELO de torbellino de la fase 9.1,
no del solver: pendiente de decidir si la ley debe imponerse sobre el
trabajo en vez de sobre Cu.

### Nota de la fase 10 (2026-08-16) — qué significaba «L1» hasta ahora

La campaña califica el meanline L0, que es lo que se compara contra las
máquinas NASA; la vía L1 nunca tuvo comprobación propia. Al añadírsela
(bloque T18 de la suite) aparecieron dos cosas que la degradaban en
silencio a L0:

1. `turbo-design` 1.4.2 no declara `requests` entre sus dependencias.
   Sin ese paquete `_try_load_turbodesign()` devolvía `False` y todo el
   sistema corría en L0 creyendo correr en L1. Fijado en
   `requirements-ci.txt`.
2. `_scm_solve` lanzaba su worker con `multiprocessing`. En Windows
   (método `spawn`) el hijo reimporta el módulo `__main__` del padre:
   cualquier script sin `if __name__ == "__main__"` se reejecutaba
   entero dentro del hijo, agotaba el timeout y degradaba a
   `meanline_L0(L1_unavailable_or_diverged)`. El CLI y `validate.py` sí
   tienen guard, así que las corridas de producto no estaban afectadas;
   la suite de verificación no lo tiene.

Las dos causas dejaron de existir en la fase 11: L1 es propio, no tiene
dependencias que declarar y corre EN PROCESO — no hay subproceso, así que
tampoco hay `__main__` que reimportar.

Con aquel L1 (TD3) el θ de referencia daba PR 2.464 y η_p 0.839 frente a
2.618 / 0.901 de L0 — el sesgo de ≈0.94 en PR anotado arriba. Ese sesgo
era del transformado por etapa de TD3, no de la física: el SCM propio de
la fase 11 da +0.13% sobre L0 en el mismo caso.
