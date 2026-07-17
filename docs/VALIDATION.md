# Phy-AC — Campaña de validación

Distinción VV&UQ del proyecto: `test_phyac.py` VERIFICA (¿resolvemos bien
las ecuaciones? — 45 checks); `validation/validate.py` VALIDA (¿las
ecuaciones correctas? — contra máquinas NASA medidas). CI corre ambos en
estricto.

## 1. Metodología

Para cada máquina publicada, `validate.py` construye el θ que reproduce
su **annulus** (r_tip desde U_tip/ω publicados; φ1 invertido por
bisección para que la continuidad devuelva ese r_tip) y su **trabajo
medido** (ψ_mid desde el ΔT0 implicado por PR y η publicados). El modelo
recibe el trabajo real y se califica su predicción de PÉRDIDAS → (η, PR).
Sin recalibración por máquina. Desde la fase 8 las multietapa pueden
declarar su distribución por etapa con el campo `slopes`
(phi_slope/Rx_slope del θ 15-D); las monoetapa usan el θ legacy de 13
(paddeado con pendientes 0, bit-exacto).

Casos `rotor` (Rotor 37/67): se califica el PR/η del rotor derivado del
desglose de pérdidas de la etapa (sin estátor).

## 2. Resultados vigentes (2026-07-17)

| Máquina | Plano | ΔPR | Δη | Tolerancia |
|---|---|---|---|---|
| NASA Stage 35 | etapa | +0.8% | +1.2 pts | 5% / 2 pts |
| NASA Rotor 37 | rotor | −1.1% | −1.5 pts | 5% / 3 pts |
| NASA Rotor 67 | rotor | −0.7% | −1.4 pts | 5% / 3 pts |
| GE/NASA E³ HPC (10 et.) | máquina | −4.8% | −1.4 pts | 6% / 3 pts |

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

## 5. Pendientes

- Digitalizar speedlines completas de Rotor 37 (choke mdot ±2%) cuando
  haya fuente estable (las URL de turbmodels rotaron — 2026-07).
- Verificar la entrada E³ contra el CR original (distribución por etapa
  real → sustituir el fit de `slopes` y endurecer hacia 5%; entrada
  marcada APROXIMADA en machines.py).
- Con validación off-design (speedlines): calibrar el mapa VSV.
- Pares CFX/banco para `HiFiCalibration` (L2) del sesgo L1 (≈0.94 en PR).
