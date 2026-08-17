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

### Recorte contra la vena (desde la fase 10)

El álabe NO se posa en los radios que trae su sección: se loftea con
4 mm de margen radial a cada extremo y se INTERSECA con el sólido de
revolución de la vena. El annulus se contrae también A LO LARGO DE LA
CUERDA, así que una punta que encaja en el centro de la fila se sale de
la carcasa en el borde de fuga. Es la misma receta que
`geometry_generator._cq_blade_trimmed`, para que las dos rutas
construyan el mismo álabe.

Holguras, con la convención del contrato: el ROTOR nace en el cubo y
termina en la línea de punta — el hueco de marcha vive entre él y el
barreno de la carcasa, que está en `tip + ε`; el ESTÁTOR cuelga del
barreno y deja ε sobre el tambor. El extremo libre se retrae además
½ vóxel (el redondeo de superficie); si ε < 1 vóxel el log avisa de que
el hueco impreso lo fija el vóxel y no el diseño.

### Raíz de ABETO (desde la fase 10)

`FirTree.cs` es un puerto línea a línea de `firtree_profile` y
`_cq_firtree_solid`. Perfil en (u = tangencial, v = profundidad desde la
plataforma), envuelto sobre el cilindro de cada radio y lofteado entre
`nSections = 5` caras desplazadas ±(L/2)·tan γ (brochado inclinado). La
plataforma sigue la línea de cubo con pendiente constante: raíz, ranura y
rim del disco comparten esa ley, así que el juego de brochado es el único
hueco entre ellos. `RootSinkMm` queda SOLO para contratos sin bloque
`firtree`.

Quien edita un lado edita el otro: `test_phyac.py` (bloque T20) compara
los dos perfiles punto a punto y `validation/parity_stl_step.py` compara
volúmenes por conjunto.

## 1b. Fallback: lámina de comba + voxMeshShell (patrón Phy-CC v3)

Se usa para contratos legacy sin `points` y para filas con espesor medio
< 2 vóxeles (el loft sólido no sobreviviría la voxelización). Receta:

1. Rotar la línea de comba por su stagger FIRMADO, envolver en el
   cilindro (φ = φ0 + y'/r), coser UNA lámina de quads.
2. `Voxels.voxMeshShell(msh, ½·espesor_medio)`: LE/TE/punta redondeados,
   imposible autointersecar. Espesor uniforme (pierde la distribución).
3. También se recorta contra la vena (el margen radial incluye la
   inflación de la lámina).
4. Espesor clampeado a 2 vóxeles con WARNING en el log — la pieza sale
   MÁS GRUESA que el diseño verificado; usar vóxel más fino.

## 2. Tambor del rotor y DISCOS

`LatticeUtils.voxRevolveZ` sobre la polilínea del hub (interpolación
clampeada; stubs de montaje en ambos extremos) menos el cilindro del
barreno (`DrumInnerRadiusMm`, del bloque structural del contrato).

`RotorDrum.voxDisc` (fase 10) es el puerto de `_cq_disc`. Tres cosas que
el chequeo de interferencias del ensamble encontró en la vía CadQuery y
que aquí se evitan por construcción:

- el **rim sigue la LÍNEA DE CUBO** con `rim_relief_mm` por debajo de la
  plataforma del álabe, no un cilindro al radio del centro: la línea de
  cubo sube y un cilindro asoma dentro del álabe aguas arriba;
- el **círculo de tirantes es COMÚN** a toda la pila — calculado disco a
  disco se corre hacia fuera etapa a etapa y los tirantes atraviesan el
  alma de los discos delanteros;
- el **casco del tambor se detiene en la banda del disco** en vez de
  atravesar el rim de lado a lado.

## 3. Carcasa con estátores colgados (+ IGV/OGV)

Revolución exterior = línea de tip + holgura + espesor de pared; cavidad
interior = tip + holgura EXTENDIDA axialmente más allá de ambos extremos
(deja el casco abierto, annulus pasante). Los estátores se unionan a la
pared interior. Los guide vanes del contrato (`igv`/`ogv` — las filas que
la física asume: pre-swirl α₁ de entrada y salida axial) cuelgan del
PRIMER y ÚLTIMO anillo respectivamente, como cualquier estátor.

Desde la fase 10 CADA anillo lleva brida en SUS DOS extremos, no solo los
dos extremos de la máquina: los anillos se apernan entre sí, y un anillo
impreso sin brida interna no se puede montar — los planos de corte
existían pero las juntas no.

## Contrato

Todo llega por `axial_compressor.json` (schema phyac-axial-2, ver
contract_schema.py y schemas/phyac-axial-2.schema.json). El entero final
del identificador es la versión MAYOR: `PhyACImport.nSCHEMA_MAJOR` la
comprueba y RECHAZA lo que no entienda en vez de leerlo a medias
rellenando defaults. Al cambiar el contrato hay que subir las dos
constantes a la vez — la suite comprueba que coinciden.
El C# es deliberadamente "tonto": no hay
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
