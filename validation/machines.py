"""
QUASAR Phy-AC · validation/machines.py
======================================
Base de datos de la campaña de VALIDACIÓN física (docs/VALIDATION.md):
compresores axiales reales con geometría y mediciones publicadas.

A diferencia de Phy-CC (θ directo), aquí cada entrada declara la ESPEC
publicada y validate.py construye el θ 13-D:

  * r_tip sale de U_tip publicado y RPM; φ1 se INVIERTE por bisección para
    que el annulus derivado por continuidad reproduzca ese r_tip.
  * ψ_mid sale del TRABAJO medido (ΔT0 desde PR y η medidos): el modelo
    recibe el trabajo real y se califica su predicción de PÉRDIDAS →
    (η, PR). Mismo criterio que la validación del WIF en Phy-CC.

Campos de cada entrada:

  spec      n_stages, RPM, mdot, T0, P0, HTR, U_tip [m/s], sigma_r,
            sigma_s, AR, Rx_est (reacción media estimada de la fuente).
            Puede quedar FUERA de los bounds del optimizador — el meanline
            no depende de ellos.
  eps_tip_mm  holgura de punta de marcha PUBLICADA [mm] (validate.py la
            inyecta en physics_core.TIP_CLEARANCE_MM durante la corrida
            de esa máquina; sin el campo se usa el default del módulo).
  kind      "rotor"   solo rotor (se califica el PR/η del rotor, derivado
                      del desglose de pérdidas de la etapa 1)
            "stage"   etapa completa (rotor + estátor)
            "machine" máquina multietapa completa
  measured  dict con PR y eta_isen O eta_poly (especificar cuál).
  tol       dict PR_rel / eta_pts propios del caso (transónicos y
            multietapa llevan tolerancias más holgadas — ver notes).
  source    de dónde salen geometría y mediciones.
  notes     advertencias del caso (transónico, plano ambiguo, etc.).

Añadir una máquina = añadir una entrada aquí. La tabla de resultados se
regenera con `python validation/validate.py`.
"""

# ---------------------------------------------------------------------------
# Criterios de éxito v1 (plan §validación)
# ---------------------------------------------------------------------------
PR_TOL_REL = 0.05      # monoetapa: |ΔPR|/PR < 5% en el punto de diseño
ETA_TOL_PTS = 0.02     # monoetapa: |Δη| < 2 puntos
# Guarda de regresión interina para CI en η (los transónicos estresan el
# modelo de choque L0; --strict no falla mientras no exceda este techo).
ETA_GUARD_PTS = 0.05

MACHINES = [
    dict(
        name="NASA Stage 35",
        kind="stage",
        spec=dict(n_stages=1, RPM=17_188.7, mdot=20.188, T0=288.15,
                  P0=101_325.0, HTR=0.70, U_tip=454.5,
                  sigma_r=1.40, sigma_s=1.30, AR=1.19, Rx_est=0.65),
        measured=dict(PR=1.82, eta_isen=0.828),
        eps_tip_mm=0.36,   # holgura de marcha familia TP-1337/1338
        #                    (~0.036 cm; ε/h ≈ 0.5%). APROXIMADO: tomado
        #                    del rotor 37 hermano — verificar contra
        #                    TP-1338 antes de endurecer tolerancias.
        tol=dict(PR_rel=0.05, eta_pts=0.02),
        source=("Reid & Moore, NASA TP-1338 (1978): 'Design and Overall "
                "Performance of Four Highly Loaded, High-Speed Inlet "
                "Stages for an Advanced High-Pressure-Ratio Core "
                "Compressor'. Etapa 35 al 100% de velocidad de diseño."),
        notes=("Etapa transónica muy cargada (M_rel punta ≈1.4, ψ≈0.44): "
               "estresa la pérdida de choque del L0. Rotor 35: 36 álabes, "
               "AR 1.19; estátor 35: 46 álabes."),
    ),
    dict(
        name="NASA Rotor 37",
        kind="rotor",
        spec=dict(n_stages=1, RPM=17_188.7, mdot=20.19, T0=288.15,
                  P0=101_325.0, HTR=0.70, U_tip=454.1,
                  sigma_r=1.50, sigma_s=1.30, AR=1.19, Rx_est=0.72),
        measured=dict(PR=2.106, eta_isen=0.877),
        eps_tip_mm=0.356,  # running tip clearance publicada del caso test
        #                    AGARD AR-355 (0.0356 cm; ε/h ≈ 0.5%)
        tol=dict(PR_rel=0.05, eta_pts=0.03),
        source=("Reid & Moore, NASA TP-1337 (1978); mediciones láser de "
                "Suder (1996) y caso test AGARD AR-355. Punto de "
                "calificación al 98% del gasto de choke."),
        notes=("Rotor aislado transónico (M_rel punta ≈1.48). Se califica "
               "contra el PR/η del ROTOR derivado del desglose de "
               "pérdidas (sin estátor). Tolerancia η ampliada a 3 pts por "
               "el modelo de choque de 1 zona."),
    ),
    dict(
        name="NASA Rotor 67",
        kind="rotor",
        spec=dict(n_stages=1, RPM=16_043.0, mdot=33.25, T0=288.15,
                  P0=101_325.0, HTR=0.375, U_tip=429.0,
                  sigma_r=1.30, sigma_s=1.20, AR=1.56, Rx_est=0.60),
        measured=dict(PR=1.63, eta_isen=0.93),
        eps_tip_mm=1.0,    # ≈1.0 mm (0.039 in) — la holgura usada por los
        #                    casos CFD estándar del rotor 67 (ε/h ≈ 0.6%)
        tol=dict(PR_rel=0.05, eta_pts=0.03),
        source=("Strazisar, Wood, Hathaway & Suder, NASA TP-2879 (1989): "
                "'Laser Anemometer Measurements in a Transonic Axial-Flow "
                "Fan Rotor'. Rotor 67, punto de pico de eficiencia."),
        notes=("Rotor de fan transónico de baja HTR (0.375, FUERA de los "
               "bounds del optimizador — el meanline lo evalúa igual). "
               "ψ≈0.53 también fuera de bounds: caso de estrés del "
               "espacio de validez declarado."),
    ),
    dict(
        name="GE/NASA E3 HPC (10 etapas)",
        kind="machine",
        spec=dict(n_stages=10, RPM=12_300.0, mdot=54.4, T0=288.15,
                  P0=101_325.0, HTR=0.50, U_tip=456.0,
                  sigma_r=1.30, sigma_s=1.20, AR=1.50, Rx_est=0.60),
        measured=dict(PR=23.0, eta_poly=0.90),
        eps_tip_mm=0.5,    # APROXIMADO: intención de diseño típica E³
        #                    (~0.020 in); el CR no da un único valor —
        #                    con ε cte las etapas traseras (h pequeña)
        #                    dominan la pérdida de holgura, como en la
        #                    máquina real.
        slopes=dict(phi_slope=-0.10, Rx_slope=0.10),
        #   Distribución por etapa (θ 15-D, fase 8): pendientes en la
        #   dirección de la práctica real de HPCs (φ cae hacia atrás, Rx
        #   crece hacia atrás), magnitud APROXIMADA — el CR-165558
        #   (detailed design) no está accesible para digitalizar la
        #   distribución real. Recuperan ~0.75 pts del déficit de PR
        #   (−5.55% → −4.80%); el resto NO es de parametrización (candidatos:
        #   cp constante a τ≈2.4, acumulación de bloqueo, WDF). Sustituir
        #   por la tabla del CR cuando haya fuente estable.
        tol=dict(PR_rel=0.06, eta_pts=0.03),
        #   PR_rel endurecida 0.08 → 0.06 con las pendientes de fase 8.
        source=("GE Aircraft Engines, NASA CR-168919 / programa Energy "
                "Efficient Engine: HPC de 10 etapas, PR 23, ~54.4 kg/s "
                "corregidos, velocidad de punta corregida ~456 m/s, "
                "radio de cubo/punta de entrada ~0.5. ENTRADA APROXIMADA: "
                "verificar contra el CR antes de endurecer tolerancias."),
        notes=("Multietapa con estátores variables y φ/Rx variables por "
               "etapa; desde la fase 8 la distribución se aproxima con "
               "pendientes lineales (campo `slopes`, fit de 2 parámetros "
               "APROXIMADO — ver comentario). Ancla la acumulación de "
               "bloqueo y el work-done factor."),
    ),
]

# ---------------------------------------------------------------------------
# Anclas de regresión internas: NO son mediciones — congelan la salida
# actual del meanline para detectar deriva silenciosa de la física.
# Actualizarlas es una decisión consciente que debe citar la corrección
# que las movió. Los valores `expect` los llena/regenera validate.py con
# --freeze-anchors tras cada calibración deliberada.
# ---------------------------------------------------------------------------
REGRESSION_ANCHORS = [
    dict(
        name="REF_AX4 (4 etapas, θ de referencia del módulo)",
        theta=[4.0, 12_500.0, 0.62, 0.55, 0.32, -0.10, 0.60,
               1.20, 1.10, 2.20, 288.15, 101_325.0, 25.0],
        expect=dict(PR=2.715996, eta_poly=0.887134, eta_isen=0.870370, T0_out=397.531727, U_tip=361.797379),          # se congela con --freeze-anchors
        rtol=1e-3,
        feasible=False,         # ψ frontal 0.35 viola el margen de Koch
    ),
]
