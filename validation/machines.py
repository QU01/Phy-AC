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
# FUERA DE DISEÑO (F-02, fase 12): calificación del MAPA contra medida.
# ---------------------------------------------------------------------------
# Hasta la fase 12 la campaña calificaba solo el PUNTO DE DISEÑO de cuatro
# máquinas. Ni una speedline, ni un límite de bombeo, ni un gasto de choke
# estaban contrastados con dato medido — y sin embargo el margen de bombeo
# es desde la fase 9 la restricción DURA que más recorta el espacio de
# diseño (en la corrida de PR 4 el optimizador cambió PR 4.04 por 3.35
# para cumplirla). Una restricción con ese peso sin calificar era el
# riesgo abierto más grande del proyecto.
#
# Cada entrada se ENGANCHA a una máquina de MACHINES por nombre (reusa su
# espec, su holgura publicada y su θ invertido) y declara cantidades de
# FUERA DE DISEÑO medidas, con sus propias tolerancias.
#
# Campos:
#   machine     nombre exacto de la entrada de MACHINES
#   speed_frac  fracción de velocidad de diseño del punto medido
#   measured    dict de cantidades medidas. Soportadas:
#                 mdot_choke        gasto de choke [kg/s] a esa velocidad
#                 mdot_stall        gasto de bombeo [kg/s] (DERIVADO de
#                                   los dos anteriores: 0.925 × 20.93)
#                 stall_over_choke  ṁ_stall / ṁ_choke (ancho del mapa).
#                                   Este ratio COMPONE los errores de sus
#                                   dos extremos: conviene mirarlo junto a
#                                   los absolutos, no en lugar de ellos.
#   tol         tolerancia OBJETIVO (lo que un meanline debería lograr)
#   guard       tolerancia INTERINA con la que --strict no falla mientras
#               la brecha esté documentada y abierta (mismo mecanismo que
#               ETA_GUARD_PTS)
#   source      cita literal y verificable
#
# LO QUE FALTA (y por qué): el PR de pico y la PENDIENTE de la speedline.
# El AGARD publica esas dos como la Figura 2.4 (curvas de PR y η frente a
# ṁ/ṁ_choke), no como tabla; los 13 puntos medidos hay que pedírselos a
# NASA («These data ... may be obtained by mailing a request to Dr K. L.
# Suder», §2.1.5). Lo que sí aparece como FRASE explícita —y es lo que se
# califica aquí— son los dos extremos del rango de gasto, que es
# justamente lo que el margen de bombeo pone en juego.
# ---------------------------------------------------------------------------
OFFDESIGN = [
    dict(
        machine="NASA Rotor 37",
        speed_frac=1.0,
        measured=dict(mdot_choke=20.93, stall_over_choke=0.925,
                      mdot_stall=19.36),
        tol=dict(mdot_choke_rel=0.05, stall_over_choke_rel=0.03,
                 mdot_stall_rel=0.05),
        guard=dict(mdot_choke_rel=0.12, stall_over_choke_rel=0.06,
                   mdot_stall_rel=0.08),
        speedline_note=(
            "ṁ_choke = 20.93 kg/s es INFERIDO, no medido: a velocidad de "
            "diseño el difusor de la instalación bloqueaba antes que el "
            "rotor, así que el máximo MEDIDO fue 20.90 kg/s y el 20.93 "
            "sale de CFD del rotor aislado más ensayos de la etapa "
            "completa (Suder, NASA TM 107310 §4.1.1). Se usa el 20.93 "
            "porque es la cifra que canoniza el caso test; la diferencia "
            "(0.14%) es menor que cualquier tolerancia en juego."),
        source=("AGARD AR-355 (Dunham ed., 1998), §2.1.4.1 «Test "
                "Conditions»: «This near stall flow rate was "
                "experimentally determined to be ṁ/ṁ_choke = 0.925 [...] "
                "The experimental ṁ_choke as determined by NASA was "
                "20.93 kg/s». Mismo documento: holgura de punta 0.0356 cm "
                "a velocidad de diseño (la que ya inyecta la entrada de "
                "MACHINES). Punto de diseño: 20.188 kg/s, PR 2.106, "
                "η_ad 0.877, 17 188.7 rpm."),
        notes=("Las dos cantidades acotan el RANGO DE GASTO del mapa al "
               "100% de velocidad: de 19.36 a 20.93 kg/s, un ancho del "
               "7.5% del gasto de choke. El punto de diseño (20.19 kg/s) "
               "cae al 96.5% del choke, así que el rotor tiene solo un "
               "3.7% de margen de gasto por encima del diseño y un 4.0% "
               "por debajo. Es exactamente lo que el margen de bombeo "
               "—restricción dura desde la fase 9— pone en juego, y hasta "
               "ahora no estaba contrastado con ningún dato medido."),
        speedline="R37_100N",
    ),
]


# ---------------------------------------------------------------------------
# SPEEDLINE medida — la característica completa, no solo sus extremos
# ---------------------------------------------------------------------------
# El AGARD publica los 13 puntos como FIGURA, pero NASA los publica como
# TABLA: el paquete de datos experimentales del Turbulence Modeling
# Resource, `rotor 37 exp data.xlsx`, hoja «map data», filas 22-34
# («rotor 37, 100% spd», cabeceras `lbm/sec | pr | tr | eta`). Los valores
# de abajo son VERBATIM de ese fichero; la conversión a kg/s
# (×0.45359237) y la normalización por ṁ_choke son aritmética nuestra.
#
# Las etiquetas «choked», «peak eff.» y «near stall» también son del
# fichero, y caen donde el AGARD dice que caen (0.98 y 0.925).
#
# Verificación independiente: estos 13 puntos se digitalizaron ANTES por
# separado desde las Figuras 2.4 y 3.1/3.2 del AGARD, y el acuerdo con la
# tabla de NASA es ≤0.003 en PR y ≤0.001 en η — dentro de la
# incertidumbre que aquella digitalización se había declarado. Dos
# caminos independientes dan lo mismo.
#
# ATENCIÓN con ṁ_choke = 20.93 kg/s. El AGARD lo presenta como medido
# («as determined by NASA»), pero la tesis de Suder (NASA TM 107310,
# §4.1.1) es explícita en que es INFERIDO: a velocidad de diseño el
# DIFUSOR de la instalación bloqueaba antes que el rotor, así que el
# máximo MEDIDO fue 20.90 kg/s y el 20.93 sale de CFD del rotor aislado
# más ensayos de la etapa completa. «the rotor is not choked at the
# highest mass flow rate in Figure 12». La entrada de OFFDESIGN sigue
# usando 20.93 porque es la cifra que el caso test canoniza, pero la
# diferencia (0.14%) es menor que cualquier tolerancia en juego y queda
# anotada aquí para que nadie la tome por una medida directa.
SPEEDLINES = {
    "R37_100N": dict(
        machine="NASA Rotor 37",
        speed_frac=1.0,
        mdot_choke=20.93,
        source=("NASA Turbulence Modeling Resource, paquete experimental "
                "del Rotor 37 (`rotor 37 exp data.xlsx`, hoja «map data», "
                "filas 22-34; rotor37-exp.zip alojado por NASA, "
                "curado por Vogel y Pederson, act. 2026-01-05). Valores "
                "VERBATIM. Corroborado por digitalización independiente "
                "de AGARD AR-355 Figs. 2.4 y 3.1/3.2 (ΔPR ≤ 0.003, "
                "Δη ≤ 0.001)."),
        #    ṁ [lbm/s]  ṁ [kg/s]   PR      TR       η_ad   etiqueta
        points=[
            (46.032, 20.880, 1.995, 1.2451, 0.890, "choked"),
            (45.920, 20.829, 1.992, 1.2448, 0.889, ""),
            (45.881, 20.811, 2.018, 1.2493, 0.891, ""),
            (45.559, 20.665, 2.065, 1.2595, 0.887, ""),
            (45.320, 20.557, 2.071, 1.2630, 0.879, ""),
            (45.238, 20.520, 2.084, 1.2656, 0.879, "peak eff."),
            (44.610, 20.235, 2.099, 1.2706, 0.872, ""),
            (44.389, 20.135, 2.110, 1.2718, 0.875, ""),
            (44.220, 20.058, 2.114, 1.2747, 0.868, ""),
            (43.670, 19.808, 2.128, 1.2797, 0.861, ""),
            (43.663, 19.805, 2.135, 1.2807, 0.862, ""),
            (42.790, 19.409, 2.141, 1.2858, 0.850, "near stall"),
            (42.747, 19.390, 2.144, 1.2871, 0.848, ""),
        ],
        # El fichero NO trae bandas de incertidumbre. La «±1.0% en PR» que
        # circula citada como tal es la Tabla 5 de Suder, que es el
        # objetivo de PRECISIÓN DESEADA PARA CFD, no el error de medida.
        unc=None,
    ),
}


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
        expect=dict(PR=2.762086, eta_poly=0.901260, eta_isen=0.886434, T0_out=396.971291, U_tip=361.801036),          # se congela con --freeze-anchors
        rtol=1e-3,
        feasible=True,         # ψ frontal 0.35 viola el margen de Koch
    ),
]
