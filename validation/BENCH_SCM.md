# Phy-AC — banco de pruebas de la fidelidad L1 (SCM)

Generado por `validation/bench_scm.py` · 80 diseños factibles muestreados por LHS (semilla 71) sobre el espacio de diseño completo.

**No es validación**: no hay ninguna máquina medida aquí. Esto caracteriza el solver contra sí mismo y contra L0. La calificación contra medida es F-02 y sigue abierta.

## 1. Cobertura

| Ley de torbellino | diseños | resueltos | estación bloqueada | limitador activo | deriva vs L0 | sin converger |
|---|---|---|---|---|---|---|
| n = -1.0 | 80 | **85%** (68) | 5 | 7 | 0 | 0 |
| n = -0.5 | 80 | **75%** (60) | 9 | 10 | 1 | 0 |

Cuando no resuelve lo DICE y `evaluate` degrada a L0 con el motivo en `source`; nunca devuelve un número que nadie debería usar.

### Cobertura por número de etapas (vórtice libre)

| etapas | diseños | resueltos |
|---|---|---|
| 1 | 9 | 100% |
| 2 | 12 | 92% |
| 3 | 6 | 83% |
| 4 | 13 | 92% |
| 5 | 13 | 85% |
| 6 | 7 | 71% |
| 7 | 12 | 75% |
| 8 | 8 | 75% |

**Este es el hallazgo central del banco.** La cobertura cae con el número de etapas, y la razón es estructural, no numérica: el annulus lo dimensiona L0 con su Cx uniforme y su densidad media, L1 resuelve un perfil, y el álabe —de ángulo fijo— convierte esa diferencia en trabajo, que cambia la densidad, que cambia la siguiente estación. En 1-4 etapas es ruido; en 7-8 se compone. La cura de fondo es que el annulus salga del MISMO solver que lo usa; mientras tanto la guarda `PR_WINDOW` rechaza el punto en vez de devolver el número.

## 2. Coste

Sobre los 68 diseños resueltos con vórtice libre:

- mediana **3.76 s** por máquina (p10 0.85, p90 6.56, máx 12.48)

| etapas | n | s/máquina (mediana) | s por etapa |
|---|---|---|---|
| 1 | 9 | 0.75 | 0.75 |
| 2 | 11 | 1.69 | 0.84 |
| 3 | 5 | 2.86 | 0.95 |
| 4 | 12 | 3.37 | 0.84 |
| 5 | 11 | 4.03 | 0.81 |
| 6 | 5 | 4.82 | 0.80 |
| 7 | 9 | 5.86 | 0.84 |
| 8 | 6 | 8.24 | 1.03 |

El coste crece de forma aproximadamente LINEAL con el número de etapas: el lazo exterior itera sobre estaciones y cada etapa añade dos.

## 3. Independencia de malla

PR con 5, 7, 9, 11 y 13 líneas de corriente sobre el mismo diseño. La dispersión es respecto al caso más fino que resolvió.

| etapas | PR 5 | PR 7 | PR 9 | PR 11 | PR 13 | dispersión |
|---|---|---|---|---|---|---|
| 2 | 1.5153 | 1.5131 | 1.5120 | 1.5110 | 1.5104 | **0.33%** |
| 5 | 2.0413 | 2.0359 | 2.0333 | 2.0314 | 2.0300 | **0.55%** |
| 1 | 1.3121 | 1.3112 | 1.3107 | 1.3104 | 1.3102 | **0.15%** |
| 7 | 4.1312 | 4.0957 | 4.0779 | 4.0624 | 4.0515 | **1.97%** |
| 2 | 1.1700 | 1.1703 | 1.1704 | 1.1704 | — | **0.03%** |

Dispersión mediana **0.33%**, máxima 1.97%. El default son 9 líneas: por debajo de 7 la derivada radial pierde resolución y por encima de 11 el coste sube sin que el resultado se mueva.

| líneas | s/máquina (mediana) |
|---|---|
| 5 | 1.90 |
| 7 | 1.87 |
| 9 | 1.90 |
| 11 | 1.91 |
| 13 | 1.92 |

## 4. Convergencia

Lazo exterior: mediana **15 iteraciones** (p90 18, máx 23 de 60 permitidas), tolerancia 1e-05 en movimiento radial relativo.

## 5. Residual L1 − L0

Es la razón de existir de la escalera de fidelidad: si L1 devolviera siempre lo mismo que L0, la capa 2 no tendría residual que aprender.

| Ley de torbellino | ΔPR mediana | ΔPR p10–p90 | Δη mediana | Δη p10–p90 | reparto de trabajo en el span |
|---|---|---|---|---|---|
| n = -1.0 | -0.08% | -2.20% … +1.31% | +0.02 pts | -0.55 … +0.66 pts | 6.1% |
| n = -0.5 | -1.50% | -6.72% … -0.12% | -0.24 pts | -1.09 … +0.51 pts | 6.8% |

Con **vórtice libre** el residual pequeño es el resultado ESPERADO y sirve de verificación: es el caso donde las hipótesis del meanline (Cx radialmente constante, trabajo uniforme) son exactas. Con **torbellino controlado** la forma cerrada del meanline es solo aproximada y ahí aparece el residual de verdad.

Correlación de ΔPR con el Mach relativo de punta de la etapa 1: **r = +0.17** (68 puntos).

---

Reproducir: `python validation/bench_scm.py --n 80 --seed 71`  ·  tiempo total de esta corrida: 12.2 min.
