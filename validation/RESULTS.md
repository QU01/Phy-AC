# Resultados de validación — Quasar Phy-AC

Generado por `validation/validate.py` el 2026-08-17 (meanline L0, sin calibración afín). **Regenerar tras cualquier cambio en `physics_core.py`** — este archivo se versiona como evidencia.

Metodología: el θ de cada máquina reproduce su annulus (r_tip vía φ1) y su TRABAJO medido (ψ desde ΔT0 publicado); se califica la predicción de pérdidas → (η, PR). Tolerancias por máquina en `machines.py` (monoetapa 5%/2 pts; transónicos y multietapa relajados — ver notas).

## Máquinas medidas

| Máquina | Plano | PR modelo | PR medido | ΔPR | η modelo | η medida | Δη [pts] | PR | η |
|---|---|---|---|---|---|---|---|---|---|
| NASA Stage 35 | etapa completa (t-t) | 1.841 | 1.820 | +1.17% | 0.845 (eta_isen) | 0.828 | +1.7 | PASS | PASS |
| NASA Rotor 37 | rotor aislado (t-t) | 2.081 | 2.106 | -1.20% | 0.861 (eta_isen) | 0.877 | -1.6 | PASS | PASS |
| NASA Rotor 67 | rotor aislado (t-t) | 1.610 | 1.630 | -1.20% | 0.905 (eta_isen) | 0.930 | -2.5 | PASS | PASS |
| GE/NASA E3 HPC (10 etapas) | máquina multietapa (t-t) | 21.717 | 23.000 | -5.58% | 0.884 (eta_poly) | 0.900 | -1.6 | PASS | PASS |

## Las mismas máquinas a fidelidad L1 (SCM)

Desde la fase 12.3 la campaña también califica L1 contra medida, llamando al SCM DIRECTO (la puerta de factibilidad de `evaluate` es una restricción de diseño, no un requisito para resolver una máquina que existe). Para los rotores aislados se compara el plano `scm["rotor1"]`, que es donde están medidos. El modelo de pérdidas de L1 es el de Koch & Smith 1976 por sección (docs/VALIDATION.md, fase 12.3).

| Máquina | PR modelo | ΔPR | η modelo | Δη [pts] | vs L0 |
|---|---|---|---|---|---|
| NASA Stage 35 | 1.838 | +0.97% | 0.849 | +2.1 | PR mejor, η peor |
| NASA Rotor 37 | 2.105 | -0.04% | 0.881 | +0.4 | PR mejor, η mejor |
| NASA Rotor 67 | — | — | — | — | no resuelve: SCMDiverged: estación BLOQUEADA: ni en el límite sónico pasa el gasto (18.366 < 33.250 kg/s) |
| GE/NASA E3 HPC (10 etapas) | — | — | — | — | no resuelve: SCMDiverged: estación BLOQUEADA: ni en el límite sónico pasa el gasto (30.049 < 54.400 kg/s) |

## Fuera de diseño (F-02)

Primera calificación del MAPA contra medida. Hasta la fase 12 la campaña solo calificaba el punto de diseño, y sin embargo el margen de bombeo es desde la fase 9 la restricción DURA que más recorta el espacio de diseño.

| Caso | Cantidad | Modelo | Medido | Δ | Objetivo | Guarda | |
|---|---|---|---|---|---|---|---|
| NASA Rotor 37 · 100% N | `mdot_choke` | 22.311 | 20.930 | +6.60% | ±5% | ±12% | FAIL |
| NASA Rotor 37 · 100% N | `stall_over_choke` | 0.904 | 0.925 | -2.23% | ±3% | ±6% | PASS |
| NASA Rotor 37 · 100% N | `mdot_stall` | 20.176 | 19.360 | +4.22% | ±5% | ±8% | PASS |

**NASA Rotor 37 · 100% N** — fuente: AGARD AR-355 (Dunham ed., 1998), §2.1.4.1 «Test Conditions»: «This near stall flow rate was experimentally determined to be ṁ/ṁ_choke = 0.925 [...] The experimental ṁ_choke as determined by NASA was 20.93 kg/s». Mismo documento: holgura de punta 0.0356 cm a velocidad de diseño (la que ya inyecta la entrada de MACHINES). Punto de diseño: 20.188 kg/s, PR 2.106, η_ad 0.877, 17 188.7 rpm.

Las dos cantidades acotan el RANGO DE GASTO del mapa al 100% de velocidad: de 19.36 a 20.93 kg/s, un ancho del 7.5% del gasto de choke. El punto de diseño (20.19 kg/s) cae al 96.5% del choke, así que el rotor tiene solo un 3.7% de margen de gasto por encima del diseño y un 4.0% por debajo. Es exactamente lo que el margen de bombeo —restricción dura desde la fase 9— pone en juego, y hasta ahora no estaba contrastado con ningún dato medido.


### Speedline medida — NASA Rotor 37 · 100% N · 13 puntos

Fuente: NASA Turbulence Modeling Resource, paquete experimental del Rotor 37 (`rotor 37 exp data.xlsx`, hoja «map data», filas 22-34; rotor37-exp.zip alojado por NASA, curado por Vogel y Pederson, act. 2026-01-05). Valores VERBATIM. Corroborado por digitalización independiente de AGARD AR-355 Figs. 2.4 y 3.1/3.2 (ΔPR ≤ 0.003, Δη ≤ 0.001).

| Métrica | Modelo | Medido | Δ | Objetivo | |
|---|---|---|---|---|---|
| `pr_max_abs` | -1.11% media | | 2.12% | 0.05 | PASS |
| `eta_max_abs` | 13 puntos | | 4.39p | 0.03 | FAIL |
| `pr_peak` | 2.120/2.144 | | -1.11% | 0.05 | PASS |
| `eta_peak` | 0.862/0.891 | | -2.89p | 0.03 | PASS |
| `slope_rel` | -0.0896/-0.0943 | | -5.0% | 0.20 | PASS |

Punto a punto:

| ṁ [kg/s] | PR med | PR mod | ΔPR | η med | η mod | Δη [pts] | SM modelo | |
|---|---|---|---|---|---|---|---|---|
| 20.880 | 1.995 | 1.990 | -0.24% | 0.890 | 0.846 | -4.39 | +0.064 | choked |
| 20.829 | 1.992 | 1.998 | +0.31% | 0.889 | 0.848 | -4.10 | +0.058 |  |
| 20.811 | 2.018 | 2.001 | -0.85% | 0.891 | 0.849 | -4.24 | +0.056 |  |
| 20.665 | 2.065 | 2.021 | -2.12% | 0.887 | 0.853 | -3.44 | +0.039 |  |
| 20.557 | 2.071 | 2.036 | -1.68% | 0.879 | 0.855 | -2.36 | +0.027 |  |
| 20.520 | 2.084 | 2.041 | -2.06% | 0.879 | 0.856 | -2.28 | +0.024 | peak eff. |
| 20.235 | 2.099 | 2.076 | -1.11% | 0.872 | 0.861 | -1.12 | -0.006 |  |
| 20.135 | 2.110 | 2.085 | -1.16% | 0.875 | 0.862 | -1.33 | -0.014 |  |
| 20.058 | 2.114 | 2.092 | -1.06% | 0.868 | 0.862 | -0.59 | -0.019 |  |
| 19.808 | 2.128 | 2.106 | -1.02% | 0.861 | 0.860 | -0.09 | -0.032 |  |
| 19.805 | 2.135 | 2.106 | -1.34% | 0.862 | 0.860 | -0.19 | -0.032 |  |
| 19.409 | 2.141 | 2.120 | -0.99% | 0.850 | 0.852 | +0.15 | -0.048 | near stall |
| 19.390 | 2.144 | 2.120 | -1.11% | 0.848 | 0.851 | +0.30 | -0.048 |  |

### Detalle por máquina

**NASA Stage 35** — fuente: Reid & Moore, NASA TP-1338 (1978): 'Design and Overall Performance of Four Highly Loaded, High-Speed Inlet Stages for an Advanced High-Pressure-Ratio Core Compressor'. Etapa 35 al 100% de velocidad de diseño.

θ construido: n=1, RPM=17189, HTR=0.700, φ1=0.514, ψ=0.437, Rx=0.65, σr=1.40, σs=1.30, AR=1.19

Notas: Etapa transónica muy cargada (M_rel punta ≈1.4, ψ≈0.44): estresa la pérdida de choque del L0. Rotor 35: 36 álabes, AR 1.19; estátor 35: 46 álabes.

**NASA Rotor 37** — fuente: Reid & Moore, NASA TP-1337 (1978); mediciones láser de Suder (1996) y caso test AGARD AR-355. Punto de calificación al 98% del gasto de choke.

θ construido: n=1, RPM=17189, HTR=0.700, φ1=0.507, ψ=0.525, Rx=0.72, σr=1.50, σs=1.30, AR=1.19

Notas: Rotor aislado transónico (M_rel punta ≈1.48). Se califica contra el PR/η del ROTOR derivado del desglose de pérdidas (sin estátor). Tolerancia η ampliada a 3 pts por el modelo de choque de 1 zona.

**NASA Rotor 67** — fuente: Strazisar, Wood, Hathaway & Suder, NASA TP-2879 (1989): 'Laser Anemometer Measurements in a Transonic Axial-Flow Fan Rotor'. Rotor 67, punto de pico de eficiencia.

θ construido: n=1, RPM=16043, HTR=0.375, φ1=0.625, ψ=0.536, Rx=0.60, σr=1.30, σs=1.20, AR=1.56

Notas: Rotor de fan transónico de baja HTR (0.375, FUERA de los bounds del optimizador — el meanline lo evalúa igual). ψ≈0.53 también fuera de bounds: caso de estrés del espacio de validez declarado.

**GE/NASA E3 HPC (10 etapas)** — fuente: GE Aircraft Engines, NASA CR-168919 / programa Energy Efficient Engine: HPC de 10 etapas, PR 23, ~54.4 kg/s corregidos, velocidad de punta corregida ~456 m/s, radio de cubo/punta de entrada ~0.5. ENTRADA APROXIMADA: verificar contra el CR antes de endurecer tolerancias.

θ construido: n=10, RPM=12300, HTR=0.500, φ1=0.491, ψ=0.418, Rx=0.60, σr=1.30, σs=1.20, AR=1.50

Notas: Multietapa con estátores variables y φ/Rx variables por etapa; desde la fase 8 la distribución se aproxima con pendientes lineales (campo `slopes`, fit de 2 parámetros APROXIMADO — ver comentario). Ancla la acumulación de bloqueo y el work-done factor.

## Anclas de regresión internas (no son mediciones)

| Ancla | Estado | Detalle |
|---|---|---|
| REF_AX4 (4 etapas, θ de referencia del módulo) | PASS | PR=2.7621 (esp. 2.7621); eta_poly=0.9013 (esp. 0.9013); eta_isen=0.8864 (esp. 0.8864); T0_out=396.9713 (esp. 396.9713); U_tip=361.8010 (esp. 361.8010) |

## Cómo añadir una máquina

1. Reunir espec publicada (U_tip, HTR, mdot, RPM, PR, η con tipo declarado) — sin datos reconstruidos de memoria.
2. Añadir la entrada a `validation/machines.py`.
3. `python validation/validate.py` y versionar este archivo.
