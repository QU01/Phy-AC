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
