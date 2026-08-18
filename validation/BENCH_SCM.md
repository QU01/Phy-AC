# Phy-AC — banco de pruebas de la fidelidad L1 (SCM)

Generado por `validation/bench_scm.py` · 80 diseños factibles muestreados por LHS (semilla 71) sobre el espacio de diseño completo.

**No es validación**: no hay ninguna máquina medida aquí. Esto caracteriza el solver contra sí mismo y contra L0. La calificación de L1 CONTRA MEDIDA la hace `validate.py` desde la fase 12.3 (tabla «Mismas máquinas a L1»), y el mapa fuera de diseño es F-02.

## 1. Cobertura

| Ley de torbellino | diseños | resueltos | estación bloqueada | limitador activo | deriva vs L0 | sin converger |
|---|---|---|---|---|---|---|
| n = -1.0 | 80 | **84%** (67) | 10 | 0 | 3 | 0 |
| n = -0.5 | 80 | **76%** (61) | 12 | 2 | 5 | 0 |

Cuando no resuelve lo DICE y `evaluate` degrada a L0 con el motivo en `source`; nunca devuelve un número que nadie debería usar.

### Cobertura por número de etapas (vórtice libre)

| etapas | diseños | resueltos |
|---|---|---|
| 1 | 9 | 100% |
| 2 | 12 | 92% |
| 3 | 6 | 83% |
| 4 | 13 | 100% |
| 5 | 13 | 85% |
| 6 | 7 | 71% |
| 7 | 12 | 75% |
| 8 | 8 | 50% |

**Este es el hallazgo central del banco.** La cobertura cae con el número de etapas, y la razón es estructural, no numérica: el annulus lo dimensiona L0 con su Cx uniforme y su densidad media, L1 resuelve un perfil, y el álabe —de ángulo fijo— convierte esa diferencia en trabajo, que cambia la densidad, que cambia la siguiente estación. En 1-4 etapas es ruido; en 7-8 se compone. La cura de fondo es que el annulus salga del MISMO solver que lo usa; mientras tanto la guarda `PR_WINDOW` rechaza el punto en vez de devolver el número.

## 2. Coste

Sobre los 67 diseños resueltos con vórtice libre:

- mediana **7.43 s** por máquina (p10 1.70, p90 13.88, máx 25.42)

| etapas | n | s/máquina (mediana) | s por etapa |
|---|---|---|---|
| 1 | 9 | 1.57 | 1.57 |
| 2 | 11 | 3.30 | 1.65 |
| 3 | 5 | 6.25 | 2.08 |
| 4 | 13 | 7.12 | 1.78 |
| 5 | 11 | 9.11 | 1.82 |
| 6 | 5 | 9.51 | 1.58 |
| 7 | 9 | 12.83 | 1.83 |
| 8 | 4 | 16.27 | 2.03 |

El coste crece de forma aproximadamente LINEAL con el número de etapas: el lazo exterior itera sobre estaciones y cada etapa añade dos.

## 3. Independencia de malla

PR con 5, 7, 9, 11 y 13 líneas de corriente sobre el mismo diseño. La dispersión es respecto al caso más fino que resolvió.

| etapas | PR 5 | PR 7 | PR 9 | PR 11 | PR 13 | dispersión |
|---|---|---|---|---|---|---|
| 2 | 1.4818 | 1.4838 | 1.4831 | 1.4812 | 1.4803 | **0.24%** |
| 5 | 1.9801 | 1.9838 | 1.9837 | 1.9813 | 1.9794 | **0.22%** |
| 1 | 1.3067 | 1.3072 | 1.3070 | 1.3065 | 1.3062 | **0.07%** |
| 7 | 3.7895 | 3.7928 | 3.7813 | 3.7594 | 3.7422 | **1.35%** |
| 2 | 1.1695 | — | 1.1693 | 1.1690 | 1.1688 | **0.06%** |

Dispersión mediana **0.22%**, máxima 1.35%. El default son 9 líneas: por debajo de 7 la derivada radial pierde resolución y por encima de 11 el coste sube sin que el resultado se mueva.

| líneas | s/máquina (mediana) |
|---|---|
| 5 | 3.92 |
| 7 | 9.16 |
| 9 | 3.91 |
| 11 | 3.66 |
| 13 | 3.61 |

## 4. Convergencia

Lazo exterior: mediana **15 iteraciones** (p90 18, máx 24 de 60 permitidas), tolerancia 1e-05 en movimiento radial relativo.

## 5. Residual L1 − L0

Es la razón de existir de la escalera de fidelidad: si L1 devolviera siempre lo mismo que L0, la capa 2 no tendría residual que aprender.

| Ley de torbellino | ΔPR mediana | ΔPR p10–p90 | Δη mediana | Δη p10–p90 | reparto de trabajo en el span |
|---|---|---|---|---|---|
| n = -1.0 | -0.84% | -6.07% … -0.06% | -1.02 pts | -2.18 … -0.29 pts | 7.2% |
| n = -0.5 | -2.70% | -9.11% … -0.45% | -1.28 pts | -2.42 … -0.49 pts | 6.2% |

Con **vórtice libre** el residual pequeño es el resultado ESPERADO y sirve de verificación: es el caso donde las hipótesis del meanline (Cx radialmente constante, trabajo uniforme) son exactas. Con **torbellino controlado** la forma cerrada del meanline es solo aproximada y ahí aparece el residual de verdad.

Correlación de ΔPR con el Mach relativo de punta de la etapa 1: **r = -0.58** (67 puntos).

---

Reproducir: `python validation/bench_scm.py --n 80 --seed 71`  ·  tiempo total de esta corrida: 22.3 min.
