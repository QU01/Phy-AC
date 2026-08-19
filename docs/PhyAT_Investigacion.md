<!--
Quasar Phy-AT · Investigación de recursos — turbinas axiales
Documento de investigación previo al arranque del proyecto Phy-AT.
Generado sobre la rama claude/phy-at-axial-turbines-q0gfh4.
-->

# 🔬 Quasar Phy-AT — Investigación de recursos para el análogo de Phy-AC en turbinas axiales

**Estado**: investigación previa (no hay código Phy-AT todavía).
**Objetivo**: reunir la literatura, los métodos, los sistemas y los
principios necesarios para construir un sistema autónomo de diseño
inverso de **turbinas axiales** con la lógica de Phy-AC: *prior físico
calibrado → ensemble profundo residual con puerta de incertidumbre →
NSGA-II restringido → informe autocontenido → contrato de geometría →
STEP/STL imprimibles*, con **escalera de fidelidades** (meanline → SCM →
CFD/datos) e **incertidumbre calibrada**.

> **Pendiente declarado — Phy-CB**: el contrato de geometría de Phy-AC
> (`phyac-axial-2`) se reutilizará en un sistema nuevo de cámaras de
> combustión (Phy-CB), al que hoy no hay acceso. El contrato de Phy-AT
> (`phyat-axial-1`) debe diseñarse con la misma disciplina (schema JSON
> versionado, validador sin dependencias, consumidor que rechaza lo que
> no entiende) y con la **interfaz combustor→turbina explícita** (perfil
> radial de temperatura OTDF/RTDF, swirl residual, sangrados de
> refrigeración), de modo que cuando Phy-CB exista, ambos contratos
> casen en la frontera. Ver §H.

<!-- ÍNDICE — se completa en la síntesis final -->

## A. El patrón Phy-AC: qué se replica y qué cambia con una turbina

*(sección redactada desde el propio repo — README.md,
docs/Quasar_PhyAC_Science.md, docs/VALIDATION.md, contract_schema.py,
AxialCompressorDesigner/.agent/axial-compressor-pattern.md)*

<!-- PENDIENTE: se rellena en la síntesis -->

## B. Fidelidad L0 — meanline de turbina axial y sistemas de pérdidas

<!-- PENDIENTE: informe del agente 1 -->

## C. Fidelidad L1 — through-flow / streamline curvature para turbinas

<!-- PENDIENTE: informe del agente 2 -->

## D. Fidelidades altas (L2/L3) y estado del arte en diseño y CFD de turbinas

<!-- PENDIENTE: informe del agente 3 -->

## E. PINNs, operadores neuronales, surrogates multifidelidad e incertidumbre calibrada

<!-- PENDIENTE: informe del agente 4 -->

## F. Geometría (perfiles, 3D, STEP/STL) y núcleo estructural de turbina

<!-- PENDIENTE: informe del agente 5 -->

## G. Recursos de software existentes

<!-- PENDIENTE: se rellena en la síntesis (TD3, MULTALL, etc.) -->

## H. El contrato `phyat-axial-1` y la frontera con Phy-CB (pendiente)

<!-- PENDIENTE: se rellena en la síntesis -->

## I. Plan de validación y máquinas de referencia

<!-- PENDIENTE: se rellena en la síntesis -->

## J. Hoja de ruta propuesta y riesgos

<!-- PENDIENTE: se rellena en la síntesis -->

## K. Referencias

<!-- PENDIENTE: consolidado de todas las secciones -->
