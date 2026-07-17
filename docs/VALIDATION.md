# Phy-AC — Campaña de validación

Distinción VV&UQ del proyecto: `test_phyac.py` VERIFICA (¿resolvemos bien
las ecuaciones? — 45 checks); `validation/validate.py` VALIDA (¿las
ecuaciones correctas? — contra máquinas NASA medidas). CI corre ambos en
estricto.

## 1. Metodología

Para cada máquina publicada, `validate.py` construye el θ 13-D que
reproduce su **annulus** (r_tip desde U_tip/ω publicados; φ1 invertido por
bisección para que la continuidad devuelva ese r_tip) y su **trabajo
medido** (ψ_mid desde el ΔT0 implicado por PR y η publicados). El modelo
recibe el trabajo real y se califica su predicción de PÉRDIDAS → (η, PR).
Sin recalibración por máquina.

Casos `rotor` (Rotor 37/67): se califica el PR/η del rotor derivado del
desglose de pérdidas de la etapa (sin estátor).

## 2. Resultados vigentes (2026-07-11)

| Máquina | Plano | ΔPR | Δη | Tolerancia |
|---|---|---|---|---|
| NASA Stage 35 | etapa | +0.9% | +1.3 pts | 5% / 2 pts |
| NASA Rotor 37 | rotor | −1.8% | −2.3 pts | 5% / 3 pts |
| NASA Rotor 67 | rotor | −1.2% | −2.5 pts | 5% / 3 pts |
| GE/NASA E³ HPC (10 et.) | máquina | −4.4% | −1.3 pts | 8% / 3 pts |

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

## 5. Pendientes

- Digitalizar speedlines completas de Rotor 37 (choke mdot ±2%) cuando
  haya fuente estable (las URL de turbmodels rotaron — 2026-07).
- Verificar la entrada E³ contra el CR original antes de endurecer su
  tolerancia (entrada marcada APROXIMADA en machines.py).
- Pares CFX/banco para `HiFiCalibration` (L2) del sesgo L1 (≈0.94 en PR).
