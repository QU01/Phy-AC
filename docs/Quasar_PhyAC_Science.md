<!--
Quasar Phy-AC · Scientific Foundations
Same format as Quasar_PhyCC_Science.md (the pattern project).
-->

# 🔬 Quasar Phy-AC — Scientific Foundations
## The physics, mathematics, and algorithms behind the Computational Engineering Model

Phy-AC designs **multistage axial compressors** autonomously. This
document collects the equations, correlations, algorithms and pseudocode
implemented in the code, with references. It is the axial counterpart of
`Quasar_PhyCC_Science.md`.

## Table of Contents

1. [Inverse Design Problem Formulation](#1-inverse-design-problem-formulation)
2. [Layer 1 · Physics L0: Stage-Stacking Meanline](#2-layer-1--physics-l0-stage-stacking-meanline)
3. [Layer 1 · Physics L1: TD3 Axial Spool and Patches](#3-layer-1--physics-l1-td3-axial-spool-and-patches)
4. [Layer 1 · Multi-Fidelity Calibration (L2)](#4-layer-1--multi-fidelity-calibration-l2)
5. [Layer 1s · Structural Core](#5-layer-1s--structural-core)
6. [Layer 0 · Public Data Policy](#6-layer-0--public-data-policy)
7. [Layers 2–4 · Surrogate, Search and Active Learning](#7-layers-24--surrogate-search-and-active-learning)
8. [Layer 5 · Geometry: from Meanline to Printable Parts](#8-layer-5--geometry-from-meanline-to-printable-parts)
9. [Complexity Analysis and Evaluation Budgets](#9-complexity-analysis-and-evaluation-budgets)
10. [Verification, Validation, and Uncertainty Quantification](#10-verification-validation-and-uncertainty-quantification)
11. [References](#11-references)

---

## 1. Inverse Design Problem Formulation

### 1.1 Design Space and Variables

The design vector is **13-dimensional**, using scalar distributions
instead of per-stage free parameters (`physics_core.DESIGN_VARS`):

| # | Variable | Range | Meaning |
|---|---|---|---|
| 0 | `n_stages` | 1–8 | stage count (rounded to integer inside `evaluate`) |
| 1 | `RPM` | 5k–25k | shaft speed |
| 2 | `HTR_in` | 0.40–0.80 | inlet hub-to-tip radius ratio |
| 3 | `phi1` | 0.35–0.80 | flow coefficient φ = Cx/U_m |
| 4 | `psi_mid` | 0.22–0.45 | mean stage loading ψ = Δh₀/U_m² |
| 5 | `psi_slope` | −0.30–0.30 | linear front→rear tilt of ψ |
| 6 | `Rx_mean` | 0.50–0.85 | mean-line degree of reaction |
| 7–8 | `sigma_r`, `sigma_s` | 0.9–1.6 / 0.8–1.5 | rotor / stator solidity |
| 9 | `AR` | 1.2–3.5 | rotor aspect ratio h/c (stator = 1.1·AR) |
| 10–12 | `T0_in`, `P0_in`, `massflow` | — | pinned by the spec (`fix_operating_point`) |

**The inlet tip radius is NOT a design variable.** Given (φ₁, HTR, RPM,
ṁ), continuity determines the inlet annulus uniquely; `_meanline` solves
r_tip by a contractive fixed point. With r_tip free the space was
over-determined: continuity imposed a φ different from the requested one
and the whole Smith-chart parameterization broke (found and fixed during
M1).

The per-stage loading is a scalar distribution:

$$\psi_i = \psi_{mid}\left(1 + s_\psi\,\frac{2i - (N-1)}{N-1}\right)$$

### 1.2 The Constrained Multi-Objective Inverse Problem

$$\min_\theta \; \bigl(-\eta_{poly}(\theta),\; |PR(\theta)-PR^*|/PR^*\bigr)
\quad \text{s.t.} \quad g(\theta) \le 0$$

with $g$ the concatenation of 8 aerodynamic constraints (§2.4), the spec
constraints (U_tip, RPM, r_tip, n_stages, optional power) and 4 exact
structural constraints (§5). Deb's constrained dominance (§7) consumes the
continuous violation magnitudes — degenerate θ must return finite,
continuous g, never raise.

### 1.3 Why the Problem is Hard (and What Makes it Tractable)

The response surface is discontinuous at choke, plateaued by the integer
stage count, and the feasible region is ~20% of the box. Tractability
comes from the same three Phy-CC pillars: a **fast physics prior** (~0.5
ms/point), **residual learning** (the surrogate learns truth − L0, so it
degrades toward physics, not noise) and **constraint-aware search** (Deb
dominance over the exact g of the prior).

---

## 2. Layer 1 · Physics L0: Stage-Stacking Meanline

### 2.1 Kinematics: Repeating-Stage Velocity Triangles

Mean-line triangles from (φ, ψ, Rx), angles measured from the axis, design
Cx constant through the machine (Dixon & Hall §5):

$$\frac{C_{u1}}{U} = 1 - R_x - \frac{\psi}{2}, \qquad
  \frac{C_{u2}}{U} = 1 - R_x + \frac{\psi}{2}$$
$$\tan\beta = \frac{U - C_u}{C_x}, \qquad \tan\alpha = \frac{C_u}{C_x},
  \qquad \Delta h_0 = \psi U^2 \;(\text{Euler})$$

The identity $\tan\beta_1 + \tan\alpha_1 = 1/\varphi$ is verified to 1e-6
in the test suite. Howell's work-done factor enters when **carving the
blade** (layer 5a metal angles), not the achieved work.

Per-station continuity solves the subsonic axial Mach by bisection of
$\dot m = \rho(M)\,A\,K_B\,C_x$ with swirl; no subsonic root ⇒ choke with
a **continuous** penalty folded into g (Deb ranking needs continuity —
verified through the choke boundary in the suite). Endwall blockage
$K_B = 0.98 - 0.005\,i$ (floor 0.90), Koch–Smith style. Annulus modes:
`const_hub | const_mean | const_tip` (default mean); each station's area
maintains the design Cx with the local state, **including the
rotor-exit station** (needed by the L1 passage, §3).

### 2.2 Loss Correlations (ω̄ referenced to row-inlet dynamic head)

* **Profile — Lieblein (1959) / NASA SP-36 ch. VII.** Equivalent
  diffusion and wake momentum thickness:

  $$D_{eq} = \frac{\cos\beta_2}{\cos\beta_1}\left[1.12 +
    0.61\,\frac{\cos^2\beta_1}{\sigma}\,|\tan\beta_1-\tan\beta_2|\right],
    \qquad \frac{\theta}{c} = \frac{0.0045}{1 - 0.95\ln D_{eq}}$$
  $$\bar\omega_p = 2\,\frac{\theta}{c}\,\frac{\sigma}{\cos\beta_2}
    \left(\frac{\cos\beta_1}{\cos\beta_2}\right)^2$$

* **Secondary + annulus — Howell (1945).** $C_{Ds} = 0.018\,C_L^2$,
  $C_{Da} = 0.020\,(s/h)$, with the drag→loss conversion

  $$\bar\omega_{ew} = C_D\,\sigma\,\frac{\cos^2\beta_1}{\cos^3\beta_m}$$

  (An inverted cosine factor here was the dominant −17 pt η error found
  during the NASA validation; the regression anchors now guard it.)

* **Shock (M_rel > 1).** Two-point normal-shock average between tip and
  mean-line inlet Mach (Miller-style — the shock only spans the outer
  blade), scaled by **K_SHOCK = 0.70**, calibrated against NASA Rotor
  37/67 (the real passage shock is oblique; the normal-shock model
  overestimates).

* **Tip clearance.** $\Delta\eta = 2.0\,(\varepsilon/h)$ with
  $\varepsilon/h = 0.015$ — the classical 2–3% η per 1% clearance
  sensitivity.

Stage efficiency $\eta_{tt} = 1 - \sum \Delta h_{loss}/\Delta h_0$; stage
PR from the local polytropic relation; machine
$\eta_{poly} = \frac{\gamma-1}{\gamma}\,\ln PR / \ln\tau$. The OGV
(last row) removes the residual swirl and pays its ω̄ as a P0 drop.

### 2.3 L0 Pseudocode

```
solve inlet r_tip by fixed point (continuity with phi1, HTR)
for stage i in 0..N-1:
    psi_i from (psi_mid, psi_slope)
    solve axial Mach at stage inlet (bisection, swirl alpha1)
    triangles from (phi, psi_i, Rx)  →  W1, W2, C1, C2, angles
    losses: Lieblein profile + Howell endwall + shock avg + clearance
    eta_tt_i, PR_i, T0/P0 update; Koch Ch and stall margin SM_i
    next-station area maintaining design Cx (local density, blockage)
OGV loss; machine totals (PR, eta_poly, power, AN², length)
g vector (8 aero constraints, continuous)
```

### 2.4 Constraints and Stability

$g \le 0$ (all continuous): max rotor DF ≤ 0.55 and stator DF ≤ 0.60
(Lieblein); de Haller ≥ 0.72 both rows; first-rotor tip M_rel ≤ 1.35
(choke penalty folded in); **Koch (1981) stall margin**: per-stage static
pressure-rise coefficient

$$C_h = \frac{\Delta h_{static}}{\tfrac12\,(W_1^2 + C_2^2)}
\quad (\Delta h_{static} = \Delta h_0 \text{ for a repeating stage})$$

with limit 0.48 (Koch's 0.45–0.55 band) and minimum margin 0.10; exit
Mach ≤ 0.55; last-stage blade height ≥ 8 mm (manufacturability).

**Off-design** (`compressor_map`): frozen metal angles and areas;
$\psi_{od} = 1 - \varphi\,(\tan\alpha_1 + \tan\beta_2)$; parabolic
incidence bucket $\omega(i) = \omega^*(1 + (i/10°)^2)$ capped ×4; surge =
Koch SM ≤ 0 or the positive-slope branch left of the PR peak; choke from
station continuity.

---

## 3. Layer 1 · Physics L1: TD3 Axial Spool and Patches

### 3.1 The Method

`turbo-design` (NASA TD3) solves the compressible marching of a
multi-row spool over the real annulus. Phy-AC builds the passage from the
L0 stations (including rotor-exit areas) and one rotor + one stator row
per stage, extracts PR and τ, and derives η_poly. L1 overwrites only
PR/η in the record; everything else (g, stage table) stays L0.

### 3.2 The Patches and Findings (De Facto Proprietary Knowledge)

Verified on TD3 **1.4.2** (spike M7, 2026-07-11), applied lazily at
import with graceful degradation to L0:

1. **Frustum area patch** (inherited from Phy-CC):
   `compute_streamline_areas` squares dx and allows negative areas; the
   conical-frustum form degenerates correctly to a cylinder for axial
   passages (dr ≈ 0), so one patch serves both machine types.
2. **`stator_calc` Yp coercion**: with a non-Pressure loss on stator
   rows, `row.Yp[:] = 0` crashes when Yp arrives as a 0-d scalar. The
   patch coerces to a 1-d array — and must be applied to BOTH
   `compressor_math` and the **by-name import** in `compressor_spool`.
3. **Meanline mode is mandatory** (`num_streamlines=1`): with more
   streamlines the radial-equilibrium ODE (`radeq`, RK45) collapses its
   step size and hangs indefinitely.
4. **Work-preserving stage transform**: TD3's pressure-balance marching
   does **not** propagate upstream swirl into the rotor's Euler work
   (each rotor effectively sees axial inflow). Each stage is therefore
   transformed conserving work: axial inlet, rotor exit swirl
   $C_{u2}^{eff} = \psi U$, stator back to axial. A damped outer loop
   (relaxation 0.5, ≤6 iterations) corrects the equivalent exit angle
   with TD3's actual Cx — the work↔density↔Cx coupling diverges without
   relaxation. τ then matches L0 to 3 decimals by construction.
5. **Loss imposition**: `FixedPressureLoss` with the L0 ω̄ rescaled to
   TD3's reference (the *upstream row's* dynamic head; relative frame for
   rotors). The Polytropic loss type is unreachable by TD3's internal Yp
   optimizer on axial spools (the upstream head bounds the achievable
   ΔP0 loss). Sign convention: rotor metal exit angles **negative**,
   stator positive (cf. TD3's own turbine example).
6. **Hard per-solve timeout** (`PHYAC_L1_TIMEOUT`, default 180 s): the
   massflow-convergence loop can hang; each solve runs in a spawned
   subprocess that is terminated on timeout and degrades to tagged L0.

### 3.3 Protocol for L1 Usage in Phy-AC

L1 runs only on **feasible** points (the L0 g gates it), results outside
the sanity window $[0.70, 1.40]\cdot PR_{L0}$ are discarded as
divergence. Measured on feasible LHS samples: **95% solve success, PR
correlation r = 0.955, systematic ratio PR_L1/PR_L0 ≈ 0.94** — a level
bias absorbed by the affine calibration (§4), exactly the role L1 plays
in the multi-fidelity ladder.

---

## 4. Layer 1 · Multi-Fidelity Calibration (L2)

Identical to Phy-CC. `HiFiCalibration` fits per-output affine corrections
$y_{cal} = a\,y + b$ on (model, hi-fi) pairs for PR and η_poly; identity
below 2 pairs; CLI `--hifi-pairs`. The nonlinear correction lives in the
surrogate's residual head once ≥15 pairs exist. `register_hifi_pair`
re-injects CFD/rig results Noyron-style.

---

## 5. Layer 1s · Structural Core

### 5.1 Contract

`evaluate_structural(theta, record, material)` — the record's
`stage_table` carries radii, chords, blade counts and local temperature
per stage; every stage is evaluated at its own metal temperature and the
**worst stage** governs.

### 5.2 Materials and Thermal Derating

Handbook-typical library (Al-2618, Al-7075, Ti-6Al-4V, Inconel-718,
17-4PH) with piecewise-linear derating k(T) on σ_y/σ_uts and a
material-dependent **AN² limit** (Ti: 4.5×10¹⁰ in²·rpm²).

### 5.3 Per-Stage Disc Model

Each stage is an annular disc (web thickness = 0.5·chord) from the bore
(0.30·r_hub) to the hub, loaded at the rim by the centrifugal blade pull
(tapered-blade volume, factor 0.55). The displacement-form axisymmetric
solver (plane stress, variable thickness, O(n) Thomas solve) is inherited
**verbatim from Phy-CC** and validated against the exact Timoshenko
annular-disc solution to <1%.

### 5.4 Product Margins (g_struct ≤ 0, hard constraints)

1. Yield at 105% overspeed (σ scales with ω²).
2. Burst margin ≥ 1.22 (area-weighted mean tangential stress, API-617
   style practice).
3. **AN² ≤ AN²_max(material)** — the industry disc-sizing metric.
4. Blade-root centrifugal stress
   $\sigma = k_{taper}\,\rho\,\omega^2\,(r_t^2 - r_h^2)/2 \le \sigma_y(T)$.

---

## 6. Layer 0 · Public Data Policy

No large public **axial**-compressor dataset exists (DATED on Zenodo is
centrifugal — 22M samples *generated by a meanline*, i.e. the same
flywheel strategy Phy-AC uses with `dataset.csv`). The public aggregates
that do exist are small NASA report tables:

* **Validation machines** (hard-coded with citations in
  `validation/machines.py`): NASA Stage 35 (TP-1338), Rotor 37 (TP-1337 /
  AGARD AR-355), Rotor 67 (TP-2879), GE/NASA E³ HPC.
* **Correlation anchors** versioned in-repo (`data/cascade/*.csv`, SHA-256
  manifest, regenerated by `data_pipeline.build()`): Lieblein θ/c(D_eq),
  Carter m(γ), Howell work-done factor. They are regression anchors of
  the implemented correlations, not digitized measurements.

The remote whitelist is deliberately **empty**: turbmodels.larc.nasa.gov
URLs rotted (verified 2026-07), so validation is offline by design. The
download+manifest machinery is kept for a future stable mirror.

---

## 7. Layers 2–4 · Surrogate, Search and Active Learning

Ported from Phy-CC **nearly verbatim** — the optimizer core is
domain-agnostic. What Phy-AC provides is the domain adapter:

* **Physics embedding** (`physics_features`): normalized θ (13) + 12
  closed-form features (U_tip, tip Mach, max rotor/stator DF, first-rotor
  M_rel, min Koch SM, exit HTR, last blade height, exit Mach, min hub
  reaction, log PR_L0, η_L0). The ensemble learns the **residual over
  L0**.
* **Deep ensemble**: K=5 NumPy MLPs (96, 96, 64) SiLU, bootstrap +
  early stopping; μ±σ across members; **quality gate** R²(PR) ≥ 0.90,
  R²(η) ≥ 0.85, ±2σ coverage ∈ [0.80, 1] before the surrogate may guide
  search (space-filling recovery rounds otherwise).
* **NSGA-II** (pop 96, gens 60) with Deb's constrained dominance over
  the exact g (aero + spec + structural), conservative LCB objectives.
* **Acquisition**: predicted-Pareto exploitation + high-σ exploration,
  k-means de-duplication, batch of 14 physical evaluations per round.
* `DesignSpec(PR_target, massflow, RPM_max, U_tip_max, r_tip_max_mm,
  n_stages_max, power_max_W, material, fixed_vars)`; `n_stages` gets the
  same continuous-variable/integer-physics treatment as `N_blades` in
  Phy-CC. `fixed_vars` pins any of the 10 design variables (e.g.
  `{"n_stages": 5, "phi1": 0.60}`): every evaluation path goes through
  `fix_operating_point`, so the pinned dimension collapses for LHS,
  NSGA-II and acquisition alike.
* **Controllability** (CLI): `--seed` (default 71), `--fix VAR=VALUE`,
  `--eval-theta` (direct evaluation of a given design, no optimization),
  `--resume` (warm start from a checkpoint — the saved dataset replaces
  the LHS and the surrogate retrains on it), `--list-pareto` /
  `--pareto-pick N` (re-verify the front of a checkpoint with the
  physics and regenerate the full deliverables for any point).

---

## 8. Layer 5 · Geometry: from Meanline to Printable Parts

### 8.1 Layer 5a (Python — turbodesigner-inspired)

Free-vortex spanwise construction ($r\,C_u = const$, $C_x = const$) at
**13 span stations per row** (exact twist at every rib — with fewer
stations the layer-5c linear interpolation left visible creases at the
original ribs). Local solidity $\sigma(r) = Z c(r)/2\pi r$ with chord
taper 0.85. Metal angles by fixed point:

$$\chi_1 = \beta_1 - i^*, \qquad \chi_2 = \beta_2 - \delta, \qquad
  \delta = m\,\theta/\sqrt{\sigma},\;\; m = 0.23(2a/c)^2 + \gamma/500$$

with Lieblein design incidence (Aungier 2003 fits of the SP-36 charts).
Profiles: circular-arc camber with NACA 65-010 tabulated thickness, or
biconvex DCA when the inlet Mach exceeds 0.8. Sections are emitted as 60
closed CCW points **plus the 41-point camber line** (both centroid
centered) with **signed** stagger (rotor +, stator −).

**Exact axial placement**: each row's slot is the envelope
$2\max_j(\max|x_{rot}|_j + t_j/2)$ over its sections' rotated camber —
the hub runs at far lower stagger than the tip, and the sheet is centered
at the profile centroid, so the naive mean-stagger placement overlapped
rotor and stator by up to ~4.5 mm (found in STL review; now
regression-tested).

### 8.2 Layer 5c (C#/PicoGK — real machine construction)

Blade rows use the **Phy-CC v3 vane recipe**: mesh only the camber
mid-surface (one quad sheet per blade, densified to ~25 ribs) and thicken
with `Voxels.voxMeshShell` — the signed-distance offset cannot
self-intersect, LE/TE/tip come out rounded at ½·thickness. Roots sink
into their body; free ends retract clearance + shell radius.

Parts (split where real bolted joints sit — mid-gap between a stator TE
and the next rotor LE):

* **Shaft**: stubs past both ends, center bore.
* **RotorStage i** (bladed disc): disc web (0.5·axial chord, clamped
  4–14 mm — the §5 model) + hub flow-path shell segment (5 mm wall,
  incl. the spacer under its stator) + its rotor blades; bore slides
  over the shaft.
* **StatorRing i**: casing shell segment + stator vanes; end rings carry
  the bolted flanges (12 × M6 on the mid bolt circle); the bleed-stage
  ring carries the **bleed port** (radial boss + through hole).
* Union views (Rotor / Casing / Assembly) for inspection.

A bill of materials (`bom.csv`) lists every part with quantities.

---

## 9. Complexity Analysis and Evaluation Budgets

Measured on a laptop CPU:

| Component | Cost |
|---|---|
| L0 meanline | ~0.5 ms/point |
| g_struct (per-stage discs, Thomas O(n)) | ~1 ms/point |
| L1 TD3 spool (subprocess incl. spawn) | ~4–6 s/point |
| Ensemble training (150–500 samples) | 5–20 s |
| NSGA-II round over surrogate (96×60, with L0 per individual) | ~15 s |
| Full quick run (L0) | ~1–2 min |
| STL build, 5 stages Ø500 mm @ 0.8 mm voxel (all parts) | ~20–25 min |

Budget comparison: the active-learning loop reaches a verified feasible
design in ~150–500 physical evaluations vs several thousand for blind
LHS at equal quality — the same ~15× ratio measured in Phy-CC, because
the machinery is the same.

---

## 10. Verification, Validation, and Uncertainty Quantification

* **Verification** (`test_phyac.py`, 47 checks): triangle identities to
  1e-6, Euler/enthalpy conservation to 0.1%, isentropic limit (η → 1
  with losses off), g continuity through choke (max jump < 0.6 over an
  81-point ṁ sweep), LHS(500) with zero exceptions and finite g,
  profile/contract invariants (counts, winding, signed staggers, **no
  axial row overlap including shell inflation**), disc solver vs
  Timoshenko < 1%, Deb dominance rules, ensemble gate on a smooth
  function, checkpoint round-trip, seed-71 reproducibility.
* **Validation** (`validation/`, no per-machine recalibration; θ rebuilt
  from published annulus + measured work, grading the loss prediction):
  Stage 35 +0.9% PR / +1.3 pts η; Rotor 37 −1.8% / −2.3; Rotor 67
  −1.2% / −2.5; E³ HPC (10 stages) −4.4% / −1.3. Regression anchors
  freeze the meanline output (`--freeze-anchors` is a deliberate act).
* **UQ**: ensemble σ gates the search (LCB) and is reported per round;
  the L1 systematic bias is quantified (≈0.94 in PR) and handled by the
  affine calibration, not hidden.

---

## 11. References

- Lieblein, S. (1959). *Loss and Stall Analysis of Compressor Cascades*.
  ASME J. Basic Eng. — and NASA SP-36 (1965), *Aerodynamic Design of
  Axial-Flow Compressors*, ch. VI–VII.
- Howell, A. R. (1945). *Fluid Dynamics of Axial Compressors*. Proc.
  IMechE 153.
- Koch, C. C. (1981). *Stalling Pressure Rise Capability of Axial Flow
  Compressor Stages*. ASME J. Eng. Power 103.
- Koch, C. C., & Smith, L. H. (1976). *Loss Sources and Magnitudes in
  Axial-Flow Compressors*. ASME J. Eng. Power 98.
- Carter, A. D. S. (1950). *The Low Speed Performance of Related
  Aerofoils in Cascade*. ARC CP-29.
- Aungier, R. H. (2003). *Axial-Flow Compressors: A Strategy for
  Aerodynamic Design and Analysis*. ASME Press.
- Dixon, S. L., & Hall, C. A. (2014). *Fluid Mechanics and Thermodynamics
  of Turbomachinery*, 7th ed., ch. 3 & 5.
- Reid, L., & Moore, R. D. (1978). NASA TP-1337 / TP-1338 (Rotors 37/35,
  Stages 35–38).
- Strazisar, A. J., Wood, J. R., Hathaway, M. D., & Suder, K. L. (1989).
  NASA TP-2879 (Rotor 67).
- Suder, K. L. (1996). NASA TM-107310; AGARD AR-355 (Rotor 37 test case).
- GE Aircraft Engines. *Energy Efficient Engine (E³) High-Pressure
  Compressor* reports, NASA CR series.
- Deb, K., et al. (2002). *A Fast and Elitist Multiobjective Genetic
  Algorithm: NSGA-II*. IEEE Trans. Evol. Comput. 6(2).
- Lakshminarayanan, B., et al. (2017). *Simple and Scalable Predictive
  Uncertainty Estimation using Deep Ensembles*. NeurIPS.
- NASA `turbo-design` (TD3) v1.4.2 — L1 spool solver.
- OpenOrion `turbodesigner` (MIT) — layer-5a/5c design reference.
- LEAP 71 `PicoGK` v2.2.0 — voxel geometry kernel (layer 5c).
