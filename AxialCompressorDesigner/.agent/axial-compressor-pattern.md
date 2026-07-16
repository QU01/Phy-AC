# Patrón: compresor axial multietapa (Phy-AC, capa 5c)

Receta PicoGK del `AxialCompressorDesigner` — el análogo axial del
`compressor-pattern.md` de Phy-CC. Tres piezas:

## 1. Filas de álabes: lámina de comba + voxMeshShell (patrón Phy-CC v3)

NUNCA construir el perfil cerrado presión+succión como malla cosida: donde
el radio local de curvatura baja de ½·espesor las dos láminas offset se
autointersecan y OpenVDB lo renderiza como esquirlas. La receta robusta
(heredada de `Blades.cs` v3 de Phy-CC):

1. La capa Python (geometry_generator.py) emite por sección la LÍNEA DE
   COMBA (25 pts, marco de cuerda, centroide en el origen) + stagger
   FIRMADO (rotor +, estátor −) + espesor.
2. Por álabe: rotar cada polilínea por su stagger, envolver en el cilindro
   de su radio de stacking (φ = φ0 + y'/r — conserva longitud de arco),
   coser UNA lámina de quads entre secciones (winding irrelevante).
3. `Voxels.voxMeshShell(msh, ½·espesor_medio)` engorda la lámina en el
   dominio de distancia con signo: LE/TE/punta redondeados, imposible
   autointersecar.
4. Raíz hundida `RootSinkMm` en su cuerpo (tambor para rotores, carcasa
   para estátores); extremo libre retraído `holgura + ½·espesor` para que
   la inflación del shell no se coma el gap de marcha.
5. Espesor mínimo 2 vóxeles o los álabes traseros desaparecen del campo.

## 2. Tambor del rotor

`LatticeUtils.voxRevolveZ` sobre la polilínea del hub (interpolación
clampeada; stubs de montaje en ambos extremos) menos el cilindro del
barreno (`DrumInnerRadiusMm`, del bloque structural del contrato).

## 3. Carcasa con estátores colgados

Revolución exterior = línea de tip + holgura + espesor de pared; cavidad
interior = tip + holgura EXTENDIDA axialmente más allá de ambos extremos
(deja el casco abierto, annulus pasante). Los estátores se unionan a la
pared interior.

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
