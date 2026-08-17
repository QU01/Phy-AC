# Phy-AC — banco de pruebas de la fidelidad L1 (SCM)

Generado por `validation/bench_scm.py` · 80 diseños factibles muestreados por LHS (semilla 71) sobre el espacio de diseño completo.

**No es validación**: no hay ninguna máquina medida aquí. Esto caracteriza el solver contra sí mismo y contra L0. La calificación de L1 CONTRA MEDIDA la hace `validate.py` desde la fase 12.3 (tabla «Mismas máquinas a L1»), y el mapa fuera de diseño es F-02.

## 1. Cobertura

| Ley de torbellino | diseños | resueltos | estación bloqueada | limitador activo | deriva vs L0 | sin converger |
|---|---|---|---|---|---|---|
| n = -1.0 | 80 | **78%** (62) | 11 | 6 | 1 | 0 |
| n = -0.5 | 80 | **70%** (56) | 13 | 9 | 2 | 0 |

Cuando no resuelve lo DICE y `evaluate` degrada a L0 con el motivo en `source`; nunca devuelve un número que nadie debería usar.

### Cobertura por número de etapas (vórtice libre)

| etapas | diseños | resueltos |
|---|---|---|
| 1 | 9 | 100% |
| 2 | 12 | 75% |
| 3 | 6 | 67% |
| 4 | 13 | 92% |
| 5 | 13 | 85% |
| 6 | 7 | 71% |
| 7 | 12 | 67% |
| 8 | 8 | 50% |

**Este es el hallazgo central del banco.** La cobertura cae con el número de etapas, y la razón es estructural, no numérica: el annulus lo dimensiona L0 con su Cx uniforme y su densidad media, L1 resuelve un perfil, y el álabe —de ángulo fijo— convierte esa diferencia en trabajo, que cambia la densidad, que cambia la siguiente estación. En 1-4 etapas es ruido; en 7-8 se compone. La cura de fondo es que el annulus salga del MISMO solver que lo usa; mientras tanto la guarda `PR_WINDOW` rechaza el punto en vez de devolver el número.

## 2. Coste

Sobre los 62 diseños resueltos con vórtice libre:

- mediana **4.33 s** por máquina (p10 1.01, p90 7.09, máx 12.73)

| etapas | n | s/máquina (mediana) | s por etapa |
|---|---|---|---|
| 1 | 9 | 0.98 | 0.98 |
| 2 | 9 | 2.03 | 1.01 |
| 3 | 4 | 3.15 | 1.05 |
| 4 | 12 | 4.23 | 1.06 |
| 5 | 11 | 4.73 | 0.95 |
| 6 | 5 | 4.98 | 0.83 |
| 7 | 8 | 6.46 | 0.92 |
| 8 | 4 | 9.54 | 1.19 |

El coste crece de forma aproximadamente LINEAL con el número de etapas: el lazo exterior itera sobre estaciones y cada etapa añade dos.

## 3. Independencia de malla

PR con 5, 7, 9, 11 y 13 líneas de corriente sobre el mismo diseño. La dispersión es respecto al caso más fino que resolvió.

| etapas | PR 5 | PR 7 | PR 9 | PR 11 | PR 13 | dispersión |
|---|---|---|---|---|---|---|
| 2 | 1.4953 | 1.4902 | 1.4869 | 1.4837 | 1.4821 | **0.89%** |
| 5 | 2.0119 | 1.9993 | 1.9932 | 1.9875 | 1.9838 | **1.42%** |
| 1 | 1.3114 | 1.3094 | 1.3083 | 1.3073 | 1.3068 | **0.35%** |
| 7 | 3.9634 | 3.8794 | 3.8338 | — | — | **3.38%** |
| 2 | 1.1703 | 1.1699 | 1.1696 | 1.1692 | 1.1689 | **0.12%** |

Dispersión mediana **0.89%**, máxima 3.38%. El default son 9 líneas: por debajo de 7 la derivada radial pierde resolución y por encima de 11 el coste sube sin que el resultado se mueva.

| líneas | s/máquina (mediana) |
|---|---|
| 5 | 2.66 |
| 7 | 2.66 |
| 9 | 2.63 |
| 11 | 2.57 |
| 13 | 2.54 |

## 4. Convergencia

Lazo exterior: mediana **15 iteraciones** (p90 18, máx 21 de 60 permitidas), tolerancia 1e-05 en movimiento radial relativo.

## 5. Residual L1 − L0

Es la razón de existir de la escalera de fidelidad: si L1 devolviera siempre lo mismo que L0, la capa 2 no tendría residual que aprender.

| Ley de torbellino | ΔPR mediana | ΔPR p10–p90 | Δη mediana | Δη p10–p90 | reparto de trabajo en el span |
|---|---|---|---|---|---|
| n = -1.0 | -0.75% | -4.83% … +0.06% | -0.62 pts | -1.75 … +0.15 pts | 11.0% |
| n = -0.5 | -2.30% | -7.38% … -0.42% | -0.97 pts | -2.02 … -0.00 pts | 10.6% |

Con **vórtice libre** el residual pequeño es el resultado ESPERADO y sirve de verificación: es el caso donde las hipótesis del meanline (Cx radialmente constante, trabajo uniforme) son exactas. Con **torbellino controlado** la forma cerrada del meanline es solo aproximada y ahí aparece el residual de verdad.

Correlación de ΔPR con el Mach relativo de punta de la etapa 1: **r = -0.43** (62 puntos).

---

Reproducir: `python validation/bench_scm.py --n 80 --seed 71`  ·  tiempo total de esta corrida: 14.2 min.
