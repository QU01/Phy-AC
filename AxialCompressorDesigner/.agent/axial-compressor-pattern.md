# Patrón: compresor axial multietapa (Phy-AC, capa 5c)

Receta PicoGK del `AxialCompressorDesigner` — el análogo axial del
`compressor-pattern.md` de Phy-CC. Tres piezas:

## 1. Filas de álabes: loft SÓLIDO del perfil real (default 2026-07-17)

El contrato lleva por sección el PERFIL CERRADO (`points`: 60 pts CCW,
distribución de espesor NACA-65/DCA real, sin autointersección por
invariante analítico de blade_profiles.py). Cada álabe se construye como
malla estanca: paredes regladas entre costillas densificadas + tapas
hub/punta trianguladas LE→TE (zipper entre las cadenas presión/succión,
independiente del índice de inicio y del winding — las secciones de
estátor llegan espejadas y reordenadas), y se voxeliza directo con
`new Voxels(msh)`. Esto conserva el radio de LE, la posición del espesor
máximo y la asimetría presión/succión que la receta de lámina destruía.

Sobre la advertencia v3 de Phy-CC ("nunca coser presión+succión"): ese
modo de fallo es de mallas construidas OFFSETEANDO la lámina de comba a
nivel de malla, que se pliegan donde el radio de comba < ½·espesor. Aquí
el contorno cerrado viene analítico de Python con el invariante de no
autointersección, así que el loft sólido es seguro.

Extremo libre retraído `holgura + ½·vóxel` (el redondeo de superficie);
raíz hundida `RootSinkMm` en su cuerpo.

## 1b. Fallback: lámina de comba + voxMeshShell (patrón Phy-CC v3)

Se usa para contratos legacy sin `points` y para filas con espesor medio
< 2 vóxeles (el loft sólido no sobreviviría la voxelización). Receta:

1. Rotar la línea de comba por su stagger FIRMADO, envolver en el
   cilindro (φ = φ0 + y'/r), coser UNA lámina de quads.
2. `Voxels.voxMeshShell(msh, ½·espesor_medio)`: LE/TE/punta redondeados,
   imposible autointersecar. Espesor uniforme (pierde la distribución).
3. Raíz hundida; extremo libre retraído `holgura + ½·espesor`.
4. Espesor clampeado a 2 vóxeles con WARNING en el log — la pieza sale
   MÁS GRUESA que el diseño verificado; usar vóxel más fino.

## 2. Tambor del rotor

`LatticeUtils.voxRevolveZ` sobre la polilínea del hub (interpolación
clampeada; stubs de montaje en ambos extremos) menos el cilindro del
barreno (`DrumInnerRadiusMm`, del bloque structural del contrato).

## 3. Carcasa con estátores colgados (+ IGV/OGV)

Revolución exterior = línea de tip + holgura + espesor de pared; cavidad
interior = tip + holgura EXTENDIDA axialmente más allá de ambos extremos
(deja el casco abierto, annulus pasante). Los estátores se unionan a la
pared interior. Los guide vanes del contrato (`igv`/`ogv` — las filas que
la física asume: pre-swirl α₁ de entrada y salida axial) cuelgan del
PRIMER y ÚLTIMO anillo respectivamente, como cualquier estátor.

## Contrato

Todo llega por `axial_compressor.json` (schema phyac-axial-1, ver
geometry_generator.py). El C# es deliberadamente "tonto": no hay
matemática de perfiles ni correlaciones aquí — solo transformar, coser,
voxelizar y exportar `<Name>_Rotor.stl`, `<Name>_Casing.stl` y el
ensamble.

## Gotchas

- Conteo de puntos de comba IGUAL en todas las secciones de una fila
  (Validate() lo exige — el loft necesita correspondencia 1:1).
- Máquinas grandes (r_tip 300+ mm) a vóxel 0.3 mm son pesadas: iterar a
  0.8 mm y bajar solo para el export final; rotor y carcasa se construyen
  en el MISMO Library.Go (los campos son independientes).
- El viewer debe quedar apagado en modo subprocess (bEndAppWithTask).
