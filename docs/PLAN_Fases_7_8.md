# Plan de trabajo — Fase 7 (dinámica estructural de álabes) y Fase 8 (θ por etapa + STEP de ensamble)

> Documento de tareas con todo el contexto necesario para ejecutar las
> fases 7 y 8 del roadmap de mejora de Phy-AC sin re-descubrir el
> proyecto. Escrito 2026-07-17, tras completar las fases 1–4.

---

## 0. Contexto del proyecto (leer antes de tocar nada)

### 0.1 Arquitectura

```
physics_core.py      L0 meanline (~0.5 ms/punto) + L1 TD3 + calibración L2.
                     evaluate(theta 13-D) → record con stage_table,
                     g (8 restricciones aero), banderas. El L0 NUNCA lanza
                     excepción con θ degenerado y g es CONTINUO (la
                     dominancia de Deb lo exige — ver §0.4).
structures_core.py   L0s estructural: materiales con derating térmico k(T),
                     solver de disco anular 1-D (validado vs Timoshenko
                     <1%), g_struct (4): fluencia a 105% overspeed, burst
                     ≥1.22, AN² ≤ AN²_max(material), σ raíz de álabe.
neural_optimizer.py  Deep ensemble K=5 (residual sobre L0, quality gate
                     R²) + NSGA-II (pop 96, gens 60, dominancia de Deb)
                     + adquisición LCB/k-means. DesignSpec con fixed_vars.
                     Checkpoint save()/load() (warm start).
blade_profiles.py    NACA-65/DCA, incidencia Lieblein/Aungier, desviación
                     de Carter GENERALIZADA con signo del camber
                     (χ₂ = β₂ − sgn(θ)·δ, fase 3).
geometry_generator.py Contrato axial_compressor.json (schema phyac-axial-1):
                     stages[i].rotor/stator + filas top-level igv/ogv
                     (fase 3), annulus, assembly, structural, secciones con
                     `points` (perfil cerrado 60 pts) y `camber_points`
                     (41 pts).
AxialCompressorDesigner/  C#/PicoGK capa 5c: loft sólido del perfil real
                     (fase 4) con fallback a lámina de comba + clamp con
                     WARNING; IGV/OGV cuelgan del primer/último anillo.
test_phyac.py        Suite de verificación (60 checks tras fase 3).
validation/          4 máquinas NASA + anclas de regresión (REF_AX4).
```

### 0.2 Estado de las fases del roadmap

| Fase | Contenido | Estado |
|---|---|---|
| 1 | Controlabilidad CLI (--fix/--seed/--eval-theta/--resume/--pareto-pick) | ✅ commit `684f2eb` |
| 2 | Corrección de Reynolds por fila (f_Re en perfil+endwall) | ✅ commit `ce0085a` |
| 3 | IGV/OGV en contrato y 3D + Carter con signo | ✅ commit `24d7a16` |
| 4 | Loft sólido del perfil real + warnings de clamp | ✅ (commit siguiente a este doc) |
| 5 | Holgura de punta por fila compartida física↔geometría (Denton / Storer-Cumpsty) | ⬜ pendiente |
| 6 | Bucket de incidencia dependiente de Mach + desviación off-design (Aungier), VSV en el mapa | ⬜ pendiente |
| **7** | **Viga rotante + Campbell + K_t de raíz (y fillet en el 3D)** | ⬜ este plan |
| **8** | **φ/ψ/Rx por etapa en θ + STEP de ensamble** | ⬜ este plan |

### 0.3 Flujo VV&UQ obligatorio para CUALQUIER cambio de física

1. `python test_phyac.py` — debe quedar verde (60 checks hoy; añadir
   checks nuevos por cada capacidad nueva).
2. `python validation/validate.py` — las 4 máquinas NASA deben seguir
   PASS. Si el ancla REF_AX4 se mueve por un cambio LEGÍTIMO de física:
   `python validation/validate.py --freeze-anchors` y **citar la
   corrección** en docs/VALIDATION.md §3 (precedente: la corrección de
   Reynolds movió η de REF_AX4 0.8900→0.8837 y se documentó).
3. Actualizar docs/Quasar_PhyAC_Science.md (sección que corresponda +
   referencias) y ambos README (History + conteo de checks).

### 0.4 Gotchas conocidos (costaron tiempo — no re-descubrir)

- **Continuidad de g**: la dominancia restringida de Deb consume la
  MAGNITUD de violación. Toda restricción nueva debe ser continua, finita
  y no lanzar jamás (θ degenerado → violación grande y continua, patrón
  `neural_optimizer.DesignSpec.constraints` líneas del try/except).
- **Índices de θ**: el punto de operación vive en θ[10..12]
  (T0_in, P0_in, massflow) y `fix_operating_point` los escribe por
  índice. Si la fase 8 añade variables, **añadirlas AL FINAL** rompe
  esos índices — ver §8.1.3 para la decisión tomada.
- **Checkpoints**: `phyac_run.json` guarda θ como lista de 13. `load()`
  reconstruye records re-evaluando L0. Cambiar NDIM invalida checkpoints
  viejos — `load()` debe detectar longitud y migrar o rechazar con
  mensaje claro (§8.1.5).
- **N_FEATURES = NDIM + 12**: el MLP del ensemble toma esa dimensión de
  entrada. Cambiar NDIM cambia la red — sin migración posible del
  surrogate (solo del dataset).
- **Anclas CRLF**: data/*.csv y manifest se protegen con .gitattributes
  (eol=lf) y `csv.writer(..., lineterminator="\n")`. No tocar.
- **C# / PicoGK en esta máquina**: 16 GB RAM (a menudo <4 GB libres).
  Vóxel 1.0 mm puede dar "Out of memory" en los anillos de carcasa;
  para smoke tests usar **vóxel 2.0 mm** y el filtro de pieza
  (`Example.dll <json> <outdir> 2.0 StatorRing1`). No recompilar
  mientras un Example.dll esté corriendo (lock del DLL). El proceso
  PicoGK abre viewer OpenGL aunque showViewer=false; si el task muere,
  el proceso queda en "Waiting for task to end" — matarlo.
- **Convención de signos del contrato**: rotor +, estátor − (stagger y
  ángulos). Las secciones de estátor llegan ESPEJADAS y con el orden de
  puntos INVERTIDO (para conservar CCW) — cualquier código que asuma
  "el punto 0 es el LE" falla en estátores (la triangulación de tapas
  de la fase 4 detecta LE/TE geométricamente por eso).
- El IGV es fila "stage=-1" y el OGV "stage=n_st" en el contrato; en C#
  `p.IgvRow`/`p.OgvRow` (RowParams sueltos, NO en `p.Rows`) — el
  emparejamiento rotor/estátor por índice de `AxialCompressor` no los ve.

---

## FASE 7 — Dinámica estructural de álabes

### Motivación (del informe de evaluación)

`structures_core` cubre el DISCO (fluencia, burst, AN²) y la tracción
centrífuga de raíz, pero **nada de dinámica**: sin frecuencias naturales,
sin diagrama de Campbell, sin cribado de flutter, y la σ de raíz se
compara con σ_y **sin factor de concentración K_t** (el 3D además no
tiene fillet de raíz → K_t real 3–5). El margen reportado es optimista
×2–5 en fatiga. Esto es exigencia estándar incluso en diseño preliminar.

### 7.1 Frecuencias naturales del álabe (viga rotante)

**Modelo**: viga en voladizo (empotrada en la raíz) con rigidización
centrífuga. Primer modo de flexión (flap):

```
f₁² = f₁,₀² + S·(N/60)²          (Southwell)
f₁,₀ = (λ₁²/2π)·√(E·I/(ρ·A·h⁴))  λ₁² = 3.516 (Euler-Bernoulli, modo 1)
S ≈ 1.6 (flap 1º; usar 1.45–1.75 como banda)
```

- Sección: usar la sección REAL del hub que ya genera
  `blade_profiles`: A e I_min del polígono de `points` (hay
  `_polygon_centroid` en geometry_generator; añadir segundo momento por
  la misma fórmula de Green). NO inventar una sección rectangular si el
  polígono ya existe. Para el meanline (sin geometría construida),
  fallback rectangular equivalente: A = c·t, I = c·t³/12 con
  t = t/c·c (TC_ROOT_R/TC_ROOT_S en geometry_generator.py:47-48).
- Datos por etapa disponibles en `record["stage_table"]`: `chord_rotor_mm`,
  `h_blade_mm`, `n_blades_rotor/stator`, `T0_in_K` (para E(T) si se
  quiere derating del módulo; con E constante basta en preliminar).
- Material: `structures_core.MATERIALS[mat]` ya tiene `E`, `rho`.
- **Entregable**: función `blade_modes(stage_row, material, RPM) ->
  dict(f1_Hz, f1_static_Hz, ...)` en structures_core, por fila
  (rotor y estátor; el estátor sin término de Southwell).
- Torsión (modo 1): opcional pero barato —
  f_t ≈ (1/4h)·√(G·J/(ρ·I_p)); usar para el cribado de flutter (7.3).

**Verificación**: contra la fórmula analítica de viga empotrada uniforme
(<2%); f₁ crece con RPM (Southwell); f₁ decrece con h (∝1/h²).

### 7.2 Diagrama de Campbell y márgenes de resonancia

**Qué calcular** (por fila de rotor, en el rango 70–105% N_diseño):

1. Cruces de f₁ (y f_t si está) con los **engine orders** k·N/60,
   k = 1..6.
2. Cruce con la **frecuencia de paso de álabes** de las filas vecinas:
   EO = n_blades del estátor aguas arriba y aguas abajo (y el IGV para
   el rotor 1 — ¡n_blades del IGV está en el contrato desde la fase 3!).
   Esta es la excitación dominante en la práctica.

**Criterio** (práctica estándar de preliminar): ningún cruce de los EO
de paso de álabes dentro de ±10% de la velocidad de diseño; para
k = 1..6, margen de frecuencia ≥ 10% a N_diseño.

**Integración**:
- Nueva restricción en g_struct: `g_campbell = 0.10 − margen_min_rel`
  (continua; margen_min_rel = min sobre filas/órdenes de
  |f₁ − EO·N/60| / (EO·N/60) evaluado a N_diseño). Empezar como
  restricción DURA solo para los EO de paso de álabes; los k bajos como
  métrica reportada (los k=1,2 son casi imposibles de esquivar en
  preliminar y bloquearían todo el espacio).
- `N_STRUCT_CONSTRAINTS` en structures_core pasa de 4 a 5 (o 6 con
  Goodman, §7.4). **Actualizar**: test T10 de test_phyac.py (cuenta
  las componentes de constraints), Science.md §5.4, y el reporte HTML
  (report_generator, sección structural).
- El número de álabes por fila sale de `stage_table` — el optimizador ya
  lo "ve" (n_blades deriva de σ y cuerda), así que la restricción de
  Campbell ES optimizable (cambia con σ, AR, RPM, n_stages).

**Entregable extra**: `figures/campbell.png` en visualization.py
(líneas f₁(N) por etapa + abanico de EOs + banda ±10%); sección en el
reporte HTML.

### 7.3 Cribado de flutter (solo bandera, no restricción dura)

Preliminar estándar: **velocidad reducida** en punta
V* = W₁_tip / (b·ω_t) con b = semicuerda de punta, ω_t = 2π·f_t.
Umbral clásico de seguridad V* ≤ ~1.4 para flexión-torsión (citar
Armstrong & Stevenson 1960, o el criterio de reduced frequency
k = ω·b/W ≥ 0.3–0.4 de la práctica de compresores). Reportar por fila
`flutter_margin = 1.4/V*` y bandera si <1. NO meterlo en g todavía —
la incertidumbre del modelo de torsión es alta; documentarlo como
métrica de cribado en Science.md.

### 7.4 K_t de raíz y chequeo de fatiga (Goodman)

1. **K_t**: con fillet de radio r_f en la raíz, K_t ≈ 1 + 2·√(t/r_f)
   es demasiado conservador para este caso; usar la forma estándar para
   filete en T (Peterson): K_t ≈ 1.5–2 con r_f/t ≈ 0.3–0.5. Implementar
   `k_t_root(r_fillet, t_root)` con el fit de Peterson para "shoulder
   fillet in bending/tension" y documentar la fuente exacta en
   Science.md. El fillet por defecto: `blade_fillet_r_mm = 2.0` en el
   bloque `assembly` del contrato (§7.5).
2. La restricción existente de raíz (structures_core, componente 4 de
   g_struct) pasa de `σ_root ≤ σ_y(T)` a `K_t·σ_root ≤ σ_y(T)`.
   **Esto mueve la factibilidad** — puede reducir el % factible del
   espacio; re-correr un LHS(500) (test T4 ya lo hace: exige ≥10%
   factible) y ajustar si hace falta documentándolo.
3. **Goodman** (HCF básico): σ_alt permisible = σ_e·(1 − σ_mean/σ_uts)
   con σ_e ≈ 0.4–0.5·σ_uts (añadir `sigma_e` a MATERIALS con valores
   handbook por aleación). σ_mean = K_t·σ_root; σ_alt de excitación es
   desconocida en preliminar → reportar el **allowable vibratory
   stress restante** (no restricción dura), que es exactamente lo que
   hace la práctica industrial en esta etapa.

### 7.5 Fillet de raíz en el 3D (capa 5c)

El modelo imprime hoy una esquina viva raíz-álabe. Receta vóxel
(sin malla nueva — todo en el dominio de distancia):

1. Contrato: `assembly.blade_fillet_r_mm` (default 2.0; 0 = sin fillet).
   Python lo emite en geometry_generator (bloque assembly, junto a
   root_sink); C# lo lee en PhyACImport → `Parameters.BladeFilletRMm`.
2. En `AxialCompressor.voxRotorStage` / `voxStatorRing`, después de
   unir cuerpo+álabes: el blend clásico de PicoGK es
   `voxSmoothen`/offset dual: `vox.Offset(+r); vox.Offset(−r)`
   (closing morfológico) **restringido a una banda** alrededor de la
   superficie del cuerpo (intersectar el resultado del closing con un
   slab raíz±(r+2·vóxel) y unirlo al original) para no engordar
   LE/TE ni la punta. Verificar API exacta de PicoGK 2.2 (`Offset` /
   `voxOffset` — está en .agent/api-quickref.md).
3. Smoke test: RotorStage1 a vóxel 2.0 con y sin fillet; inspección
   visual del STL (el fillet debe verse en la unión raíz-álabe y NO en
   la punta).
4. Consistencia: el K_t de 7.4 usa `blade_fillet_r_mm` del MISMO bloque
   assembly — física y geometría comparten el parámetro (el patrón que
   la fase 5 pedirá para la holgura).

### 7.6 Checklist de cierre de fase 7

- [ ] `blade_modes` + verificación analítica (3+ checks nuevos).
- [ ] Campbell: cruces + margen continuo en g_struct (test de conteo
      actualizado; check de continuidad del margen con RPM).
- [ ] Flutter V* reportado por fila (check de rango sano).
- [ ] K_t + Goodman en structures_core; MATERIALS con σ_e.
- [ ] Fillet en contrato + C# + smoke STL.
- [ ] report_generator: sección Campbell/fatiga; visualization: figura.
- [ ] Science.md §5 ampliado con las fórmulas y referencias
      (Southwell, Peterson, Goodman, Armstrong & Stevenson).
- [ ] Suite y validación verdes; anclas NO deben moverse (la aero no
      cambia); READMEs actualizados (History + conteo de checks).

---

## FASE 8 — θ por etapa y STEP de ensamble

### Motivación

La parametrización actual (ψ_mid+ψ_slope, un solo φ/Rx/σ/AR para toda
la máquina) es el límite declarado nº1 del proyecto: el E³ HPC valida
con tolerancia relajada (±8%) precisamente porque las máquinas reales
varían φ/Rx por etapa. Además el único export CAD real es un STEP de
demostración de UN álabe (`try_export_step`, loft ruled) — inútil para
re-CAD del ensamble.

### 8.1 Parametrización por etapa (8a: pendientes; 8b: modo experto)

**Decisión de diseño recomendada — dos niveles:**

#### 8.1.1 Nivel 8a: añadir pendientes lineales (13-D → 15-D)

Añadir SOLO 2 variables al espacio de búsqueda:

```python
("phi_slope",   -0.25, 0.25),   # φ_i = phi1·(1 + s_φ·ξ_i)
("Rx_slope",    -0.20, 0.20),   # Rx_i = Rx_mean·(1 + s_Rx·ξ_i)
# con ξ_i = 2i/(N−1) − 1 (mismo esquema que psi_slope)
```

Esto captura el patrón real dominante (φ cae hacia atrás, Rx sube) sin
explotar la dimensión (NSGA-II pop 96 aguanta 15-D; con vectores por
etapa libres serían 8×3 = 24+ dims y el presupuesto de 150–500
evaluaciones no alcanza — NO hacerlo en el espacio de búsqueda).

#### 8.1.2 Nivel 8b: override por etapa para el modo experto

`evaluate()`/`_meanline` aceptan un dict opcional
`per_stage = {"phi": [...], "psi": [...], "Rx": [...]}` (longitud
n_stages) que **sobreescribe** las distribuciones escalares. No entra
al espacio de búsqueda: se usa desde `--eval-theta` (extender el CLI
con `--per-stage archivo.json`), desde la validación (E³ con su
distribución real → endurecer tolerancia, ver 8.1.6) y desde el modo
cascada que pida un usuario experto.

#### 8.1.3 Reglas de compatibilidad de índices (CRÍTICO)

- **Añadir las 2 variables AL FINAL** de DESIGN_VARS
  (índices 13, 14). Así θ[10..12] (T0/P0/ṁ) NO se mueven y
  `fix_operating_point`, el CLI (--eval-theta con 10 o 13 valores) y
  los checkpoints viejos siguen siendo interpretables.
- `denormalize/normalize/BOUNDS` se auto-ajustan (derivan de
  DESIGN_VARS).
- **Padding**: todo consumidor de θ debe aceptar 13-D y completar
  `phi_slope = Rx_slope = 0.0` (equivalencia EXACTA con el
  comportamiento actual — es el criterio de aceptación 8.1.7).
  Implementar `physics_core.pad_theta(theta) -> 15-D` y usarlo en
  evaluate/_meanline/physics_features/evaluate_design/load().
- `--eval-theta` acepta entonces 10, 12 (10+2), 13 (viejo completo) o
  15 valores — documentar el orden en el help.

#### 8.1.4 Cambios en _meanline

En el bucle de etapas (physics_core.py, sección 2):
```python
xi = 2.0*i/(n_st-1) - 1.0  si n_st>1 else 0.0
phi_i = phi1 * (1 + phi_slope*xi)      # clamp a [0.2, 1.2]·phi1? NO —
Rx_i  = Rx    * (1 + Rx_slope*xi)      # dejar que g lo castigue
```
- `tan_a1`, `Cu1/Cu2`, y el annulus de entrada usan φ de la PRIMERA
  etapa (`phi_i(i=0)`) — cuidado: el punto fijo del r_tip de entrada
  (línea ~409) usa `phi1`; con pendiente, la etapa 0 tiene
  φ₀ = phi1·(1−s_φ) — usar ese.
- `Cx_dsg` del área de la siguiente estación (línea ~621) usa
  `phi1·ω·r_mean` — pasar a `phi_i` de la etapa siguiente.
- El dict `frozen` (off-design) ya congela ángulos POR ETAPA — el
  off-design es automáticamente compatible.
- `physics_features`: NDIM sube → N_FEATURES sube 2. Sin cambio de
  código (deriva de VAR_NAMES) pero el surrogate cambia de tamaño.

#### 8.1.5 Migración de checkpoints y CLI

- `AutonomousAxialDesigner.load()`: si `len(theta) == 13`, aplicar
  `pad_theta` y continuar (log informativo). El seed/history se
  conservan.
- `--fix phi_slope=0 --fix Rx_slope=0` reproduce el espacio viejo
  (útil para comparar A/B — mencionarlo en README).

#### 8.1.6 Validación

- Re-correr las 4 máquinas: Stage 35/Rotor 37/67 son 1 etapa/rotor —
  invariantes con slopes=0 (el θ construido por validate.py se paddea).
- **E³ HPC**: buscar en el CR de GE (CR-165558, tabla de diseño por
  etapa) las distribuciones reales de φ/Rx; meterlas vía `per_stage`
  (8b) y **endurecer la tolerancia de 8% hacia 5%**. Si el CR no está
  accesible, ajustar slopes lineales a mano (fit de 2 parámetros) y
  documentar en VALIDATION.md §5 lo que se usó.
- REF_AX4: con padding la salida debe ser BIT-EXACTA → las anclas NO
  se congelan de nuevo. Si se mueven, hay un bug de padding.

#### 8.1.7 Criterios de aceptación 8a/8b

- [ ] θ 13-D paddeado ≡ resultado actual (assert numérico 1e-12 en un
      check nuevo; anclas intactas).
- [ ] LHS(500) en 15-D: 0 excepciones, g finito, ≥10% factible
      (T4 ampliado).
- [ ] `--fix`/checkpoints viejos/`--eval-theta` viejos funcionan.
- [ ] E³ con per_stage real: tolerancia endurecida y PASS.
- [ ] Un run `--quick` completo en 15-D produce diseño factible.

### 8.2 STEP de ensamble (CadQuery)

**Estado actual**: `geometry_generator.try_export_step` exporta UN álabe
del rotor 1 con `loft(ruled=True)` — declaradamente una demo.

**Objetivo**: paquete STEP del ensamble completo, apto para re-CAD
(NX/CATIA/FreeCAD), manteniendo CadQuery como dependencia OPCIONAL con
degradación silenciosa (patrón actual: `try` en el import, False si
falta).

Tareas:

1. **Sólidos de revolución** desde las polilíneas del contrato:
   - Eje: perfil (r_shaft, z) con stubs y barreno → `revolve`.
   - Casco del hub: polilínea hub − pared HubWallMm.
   - Carcasa: tip + holgura y + pared, con bridas (anillos con
     patrón de barrenos `polarArray`).
   Los parámetros están TODOS en el contrato (`assembly`,
   `structural.drum_inner_r_mm`, `annulus`) — no recalcular nada.
2. **Álabes**: loft de las secciones `points` rotadas por su stagger y
   posicionadas en (r, z_center) — mismo pipeline matemático que
   BladeRow.cs pero en CadQuery (`Workplane.polyline().close()` por
   sección sobre planos a r constante… ojo: CadQuery loftea entre
   wires en planos PARALELOS; el wrap cilíndrico real no es
   representable con loft plano → usar la aproximación de planos
   z=r_i (como el try_export_step actual) que es la estándar para
   re-CAD, y documentar que el wrap exacto queda en el STL).
   `ruled=False` (spline) ahora que las secciones son consistentes.
3. **Patrón circular + ensamble**: `cq.Assembly()` con un objeto por
   pieza física (Shaft, RotorStage{i} = disco+álabes×N, StatorRing{i},
   IGV/OGV en sus anillos), nombres = los mismos de los STL, y export
   único `assembly.step` + por-pieza `parts/*.step`.
4. **CLI**: `--step` (independiente de `--stl`); en `emit_deliverables`
   tras la geometría. Aviso de coste: un ensamble de 5 etapas ×
   40-60 álabes/fila puede tardar minutos y producir cientos de MB —
   añadir `--step-mode {assembly,parts,blade0}` con default `parts`.
5. **Tests**: check condicionado a `import cadquery` (skip silencioso
   si no está, patrón del check de matplotlib si existe); si está:
   exportar blade0 + 1 anillo de un contrato mini y verificar que el
   archivo existe y pesa >0; no validar contenido STEP (sin lector).
6. **Docs**: README (tabla de outputs + dependencia opcional `step`),
   Science.md §8.1.

### 8.3 Checklist de cierre de fase 8

- [ ] 8.1.7 completo (padding bit-exacto, anclas intactas, E³
      endurecido).
- [ ] Science.md §1.1 (tabla de variables 15-D) y §7 (DesignSpec)
      actualizados; README "Current Status" — quitar el bullet de
      "repeating-stage parameterization" o reescribirlo como resuelto
      a nivel de pendientes.
- [ ] STEP de ensamble con CadQuery opcional + `--step` + tests
      condicionales.
- [ ] Suite/validación verdes; History en ambos README.

---

## Orden de ejecución sugerido y estimación

1. **7.1 → 7.2 → 7.4** (núcleo estructural puro Python, sin tocar aero:
   riesgo bajo, anclas intactas). Después 7.3 (métrica) y 7.5 (C#).
2. **8.1a** (pendientes 15-D — toca physics_core en profundidad: hacerlo
   con la suite abierta y el criterio de padding bit-exacto desde el
   primer momento). Después 8.1b (per_stage), 8.1.6 (E³) y al final 8.2
   (STEP, independiente de todo lo demás — puede paralelizarse).
3. Un commit por sub-bloque coherente (patrón de las fases 1–4:
   commits `684f2eb`, `ce0085a`, `24d7a16`), con el dance de archivos
   compartidos si un archivo mezcla fases.

Referencias a añadir en Science.md §11 cuando se implementen:
Southwell (1921); Peterson, *Stress Concentration Factors*; Armstrong &
Stevenson (1960), *Some Practical Aspects of Compressor Blade
Vibration*; GE/NASA CR-165558 (E³ HPC detailed design report).
