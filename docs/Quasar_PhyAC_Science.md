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
3. [Layer 1 · Physics L1: Streamline-Curvature Through-Flow](#3-layer-1--physics-l1-streamline-curvature-through-flow)
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

The design vector is **15-dimensional** (13 legacy + the phase-8 per-stage
slopes; `physics_core.DESIGN_VARS`):

| # | Variable | Range | Meaning |
|---|---|---|---|
| 0 | `n_stages` | 1–8 | stage count (rounded to integer inside `evaluate`) |
| 1 | `RPM` | 5k–25k | shaft speed |
| 2 | `HTR_in` | 0.40–0.80 | inlet hub-to-tip radius ratio |
| 3 | `phi1` | 0.35–0.80 | mean flow coefficient φ = Cx/U_m |
| 4 | `psi_mid` | 0.22–0.45 | mean stage loading ψ = Δh₀/U_m² |
| 5 | `psi_slope` | −0.30–0.30 | linear front→rear tilt of ψ |
| 6 | `Rx_mean` | 0.50–0.85 | mean-line degree of reaction |
| 7–8 | `sigma_r`, `sigma_s` | 0.9–1.6 / 0.8–1.5 | rotor / stator solidity |
| 9 | `AR` | 1.2–3.5 | rotor aspect ratio h/c (stator = 1.1·AR) |
| 10–12 | `T0_in`, `P0_in`, `massflow` | — | pinned by the spec (`fix_operating_point`) |
| 13 | `phi_slope` | −0.25–0.25 | linear front→rear tilt of φ (phase 8) |
| 14 | `Rx_slope` | −0.20–0.20 | linear front→rear tilt of Rx (phase 8) |

**Per-stage distributions (phase 8).** φ_i = φ1·(1+s_φ·ξ_i),
Rx_i = Rx·(1+s_Rx·ξ_i), ξ_i = 2i/(N−1)−1 — the same scheme ψ always used.
The slopes were appended **at the end** of θ so the operating-point
indices 10–12 never move: legacy 13-D vectors (old checkpoints, anchors,
`--eval-theta` inputs) are padded by `pad_theta` with zero slopes, which
reproduces the pre-phase-8 physics **bit-exactly** (verified by a 1e-12
check and by the untouched regression anchors). For arbitrary per-stage
distributions there is the expert-mode `per_stage` override
({"phi"/"psi"/"Rx": [...]}, `--eval-theta --per-stage file.json`,
meanline-only) — it is deliberately NOT part of the search space (a
free 24-D+ per-stage vector is intractable for the 150–500-evaluation
budget; the linear slopes capture the first-order pattern of real
machines: φ falls rearward, Rx rises).

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

* **Tip clearance (per-row, 2026-07-17).** The clearance is an ABSOLUTE
  $\varepsilon$ in mm (`TIP_CLEARANCE_MM`, default 0.4, env-overridable
  via `PHYAC_TIP_CLEARANCE_MM`) — mechanically set by tolerances and
  thermal growth, so $\varepsilon/h$ **grows toward the rear stages** as
  h falls, which a fixed ratio erased. Per rotor row,
  $\Delta h_{cl}/\Delta h_0 = f(\varepsilon/h)$, a continuous piecewise
  ramp with the regimes of Sakulkaew, Tan et al. (2013): half slope below
  $\varepsilon/h \approx 0.8\%$ (clearance optimum — casing shear
  competes with leakage), linear slope `K_TIP_CLEARANCE` = 1.8 in
  0.8–3.4% (the measured ~1.6 pts η per 1%), half slope above (tip
  unloads); $\varepsilon/h$ clamped at 8% for degenerate θ (g stays
  continuous). Mechanism per Denton (1993) / Storer & Cumpsty (1991-94).
  The geometry contract emits **the same** $\varepsilon$
  (`annulus.tip_clearance_mm`), and the validation injects each NASA
  machine's published running clearance. Removing the old uniform 1.5%
  debit required a documented **global** recalibration of
  `K_ENDWALL` 1.0 → 1.4 (Howell's simple $C_{Da}$ underpredicts endwall
  loss — Koch & Smith 1976 — and the uniform debit was silently
  absorbing it); no per-machine tuning.

* **Reynolds correction (2026-07-16).** Profile and endwall losses (the
  friction-driven terms — not shock, not clearance) are multiplied by

  $$f_{Re} = \begin{cases}
    1 & Re_c \ge 10^6\\
    (10^6/Re_c)^{0.2} & 2\times10^5 \le Re_c < 10^6\\
    5^{0.2}\,(2\times10^5/Re_c)^{0.5} & Re_c < 2\times10^5
  \end{cases}$$

  with $Re_c = \rho_1 W_1 c/\mu$ per row (stage-inlet density for both
  rows — slightly conservative for the stator). Nominal point and the
  turbulent exponent follow Koch & Smith (1976); the machine-level
  sensitivity $(1-\eta) \propto Re^{-n}$, $n \approx 0.1$–$0.2$, is the
  Wassell (1968) / Schäffler (1980) band; the laminar branch below
  $Re \approx 2\times10^5$ reflects laminar-separation loss growth. **No
  credit above $10^6$** — deliberately conservative so the NASA
  calibration (machines at $Re_c \gtrsim 10^6$) is preserved; the
  correction activates for the small machines in the design space
  (measured: −4.8 pts η at ṁ = 2.5 kg/s vs 25 kg/s at equal θ). The
  REF_AX4 regression anchor was re-frozen citing this change.

Stage efficiency, stage PR and machine efficiency: see 2.2b - since
phase 9 they come from an entropy march with a calorically imperfect gas,
not from the constant-cp polytropic formulas. The IGV (first row,
axial -> alpha_1) and the OGV (last row, alpha_1 -> axial) both pay their
omega_bar as a P0 drop and their chord as machine length.

### 2.2b Real gas, entropy bookkeeping, stall and end walls (phase 9)

**Calorically imperfect gas.** $c_p(T)$ is a quadratic fit of the JANAF
tables (<1% over 250-1000 K) and $\gamma = c_p/(c_p-R)$. Three state
functions carry the stacking:

$$h(T)=\int c_p\,dT,\qquad \varphi(T)=\int \frac{c_p}{T}\,dT,\qquad
s_2-s_1 = \varphi(T_2)-\varphi(T_1) - R\ln\frac{p_2}{p_1}$$

from which the efficiencies are **exact** for an imperfect gas - the
generalization of Dixon & Hall 1.11:

$$\eta_{poly}=\frac{R\ln PR}{\varphi(T_2)-\varphi(T_1)},\qquad
\eta_{isen}=\frac{h(T_{2s})-h(T_1)}{h(T_2)-h(T_1)},\quad
\varphi(T_{2s})=\varphi(T_1)+R\ln PR$$

At $\tau \approx 2.4$ (the E3 HPC) $c_p$ rises ~7% and $\gamma$ falls to
1.36; with a constant $c_p$ that error goes straight into $PR$.

**Loss to work by entropy** (Dixon & Hall 5.5, Eq. 5.4-5.9). Each row's
$\bar\omega$ is referred to the **compressible** dynamic head $P_0-p$
(the Koch & Smith 1976 definition - using $\tfrac12\rho W^2$ understates
the loss by 60% at $M_{rel}\approx1.4$). The stage is then marched with
real total pressures: relative total across the rotor (conserved at
constant radius), absolute across the stator, so $PR_i$ *falls out of the
march* instead of a formula. The work-equivalent of a loss is

$$\Delta h_{loss} = T_{03}\,\Delta s,\qquad
\Delta s = R\ln\frac{P_{0,in}}{P_{0,out}}$$

evaluated at the **stage exit** temperature. The previous form
$\Delta h = \bar\omega\,\tfrac12 W_1^2$ is $\Delta P_0/\rho$ at the row
inlet, which understates hot rear stages by the factor $T_{03}/T_{0,row}$
(15-25% per stage, compounding through a multistage machine).

**Stalling pressure rise - Koch (1981).** The constant
`CH_STALL_MAX = 0.48` is replaced by the real correlation. The stage's
enthalpy-equivalent static-pressure-rise coefficient uses the
**isentropic** static rise from $p_1$ to $p_3$ (Koch's definition; using
$\Delta h_0$ overstates $C_h$ by $1/\eta \approx 1.12$):

$$C_h=\frac{h(T_s)-h(T_1)}{\tfrac12 (W_1^2+C_2^2)},\qquad
\varphi(T_s)=\varphi(T_1)+R\ln\frac{p_3}{p_1}$$

and stall is reached when $C_h/\mathfrak{F}_{ef}$ meets

$$C_{h,stall}^{ef} = \bigl[0.35+0.145\ln(L/g_2)\bigr]\cdot f_{Re}\cdot
f_\epsilon \cdot f_{\Delta z}$$

with $L/g_2 = \frac{\theta/2}{\sin(\theta/2)}\cdot\frac{\sigma}{\cos\beta_2}$
the cascade diffusion parameter (circular-arc meanline length over exit
staggered spacing), the three correction factors from Koch's Figs. 4-6
(Reynolds, tip clearance $\epsilon/g$, axial spacing $\Delta z/s$), and

$$\mathfrak{F}_{ef} = \frac{V^2 + 2.5\,V_{min}^2 + 0.5\,U^2}{4V^2}$$

the **effective dynamic head factor** (his Fig. 13), where $V_{min}$ is
the lowest inlet velocity the row can see in the presence of upstream
wakes and wall boundary layers ($V\sin(\alpha+\beta)$ when
$\alpha+\beta\le90$ deg). This is the term that explains why low-stagger
stages stall early: they cannot re-energize the low-momentum fluid that
reaches them.

**End walls - Koch & Smith (1976).** Howell's annulus drag
($C_{Da}=0.020\,s/h$) and the `K_ENDWALL = 1.4` multiplier that
compensated for it are gone. In their place, the sum of wall displacement
thicknesses of their Eq. (3) and Fig. 8,

$$\frac{2\delta^*}{g} = 0.16\,x^3 + 2\frac{\epsilon}{g}x,
\qquad x=\frac{C_h^{ef}}{C_{h,stall}^{ef}}$$

gives **both** the annulus blockage (which now *emerges* from loading
instead of being the invented line $0.98-0.005i$) and the efficiency
debit of their Eq. (2):

$$\eta = \tilde\eta\,\frac{1-\sum\delta^*/h}{1-\sum\nu/h},
\qquad \sum\nu \approx 0.48\sum\delta^*$$

A Mach factor on the **profile** loss closes the set: their Fig. 6 shows
the rotor loss coefficient nearly doubling between $M_1=0.1$ and 1.5, an
effect Lieblein's incompressible correlation cannot know:
$f_M = 1+0.50\,(M^2-0.0625)$.

### 2.2c Off-design that answers the engineer's question (phase 9)

A compressor map without a **working line** cannot state a surge margin -
there is no denominator. With a fixed-area choked exit nozzle, Dixon &
Hall 5.9 Eq. (5.26b) ties inlet non-dimensional flow to pressure ratio:

$$\frac{\dot m\sqrt{c_pT_{01}}}{D^2p_{01}} = C\,
\left(\frac{p_{0e}}{p_{01}}\right)^{1-\frac{\gamma-1}{2\gamma\eta_p}}$$

so $\dot m \propto PR^{\,k}$, with $C$ fixed by requiring the line to
pass through the design point. Each speedline is then bounded by its own
physics rather than swept over an arbitrary flow range:

* **choke** is found by bisection on the choke flag and *limits the mass
  flow* - the speedline goes vertical. (Previously the meanline kept
  computing past choke and returned $PR<1$: a compressor that expands.)
* **surge** is whichever comes first walking down from choke: Koch's
  capability exhausted ($\min SM \le 0$) **or** zero/positive speedline
  slope (Dixon 5.11, Fig. 5.14). The second criterion is essential -
  incidence-driven loss growth can flatten the line before $C_h$ reaches
  its limit, and the first criterion alone would then never fire.

The reported margin is the one an engineer writes in a specification:

$$SM = \frac{PR_s/\dot m_s}{PR_{wl}/\dot m_{wl}} - 1,\qquad
SM_{\dot m} = \frac{\dot m_{wl}-\dot m_s}{\dot m_{wl}}$$

and $SM_{\dot m} \ge 15\%$ is a **hard constraint** of the optimizer.
The measured correlation between the design-point Koch margin and the
working-line flow margin is only $r = 0.57$ over feasible LHS designs -
which is precisely why constraining the former did not deliver the
latter, and why the map of a "verified" design used to have zero width.

**Bleed** is now physics, not just a hole: a per-stage extracted fraction
reduces the mass flow the downstream stages see, and the geometry
contract places the casing port **behind the stage that stalls first**,
so the printed part and the aerodynamic model finally point at the same
place.

### 2.3 L0 Pseudocode

```
solve inlet r_tip by fixed point (continuity with phi1, HTR)
for stage i in 0..N-1:
    psi_i from (psi_mid, psi_slope)
    solve axial Mach at stage inlet (bisection, swirl alpha1)
    triangles from (phi, psi_i, Rx)  →  W1, W2, C1, C2, angles
    cp_i, gamma_i = gas properties at local T0         # phase 9
    losses: Lieblein profile (x Mach factor) + Howell secondary
            + shock avg + tip clearance                # all as omega_bar
    entropy march: P01rel -> P02rel -> P02 -> P03      # real totals
    eta_fs = 1 - T03*ds_fs/dh0                         # Dixon Eq. 5.4
    Koch: Ch (isentropic static rise)/F_ef vs Ch_stall(L/g2, Re, eps, dz)
    Koch & Smith: 2*delta*/g -> blockage for next station AND eta debit
    T0/P0 update; mdot *= (1 - bleed_i)
    next-station area maintaining design Cx (local density, blockage)
IGV and OGV losses; machine totals (PR, eta_poly via phi(T), power, AN²)
g vector (9 aero constraints, continuous)
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

**Off-design** (`compressor_map`, upgraded 2026-07-17): frozen metal
angles and areas; $\psi_{od} = 1 - \varphi\,(\tan\alpha_1 +
\tan\beta_2)$; surge = Koch SM ≤ 0 or the positive-slope branch left of
the PR peak; choke from station continuity. Three physical effects on
top of the frozen-geometry model:

* **Mach-dependent incidence bucket** (`_incidence_bucket`): parabolic
  loss multiplier $\omega(i) = \omega^*(1 + (i/W)^2)$ capped ×4, with
  half-width $W(M_1)$ falling linearly from 10° at $M_1 \le 0.2$ (the
  classical low-speed bucket) to 3.5° at $M_1 \ge 0.8$ — the useful
  incidence range narrows drastically with inlet Mach (Aungier 2003,
  off-design performance; his full ranges also depend on θ and σ — the
  Mach dependence is the first-order term). The negative-incidence
  (choke-side) branch is 1.5× wider than the stall side, standard
  practice.
* **Off-design deviation** (Creveling 1968 / SP-36 practice): with
  positive incidence the flow underturns progressively,
  $\Delta\delta = 0.30\,i^+$ (capped 10°), applied to the frozen rotor
  exit angle **before** the emergent ψ — the achieved work drops toward
  stall and the speedline flattens realistically. At design $i = 0$, so
  the design point is exactly invariant (the regression anchors did not
  move — no refreeze).
* **Variable stators** (`vsv="auto"`, CLI `--map-vsv`): without VSVs the
  part-speed lines of a PR ≳ 4 fixed-geometry machine are **unphysical
  extrapolation** — the front stages sit in deep stall (the historical
  reason for the J79's variable stators). The auto schedule closes the
  front VSVs by $\Delta = 50°\,(1 - N/N_d)$ (capped 35°, only below
  design speed), full angle on stage 1 fading linearly to zero at
  mid-machine: closing a front stator adds pre-swirl to the following
  rotor and unloads it. Off-design shock loss (Boyer & O'Brien 2003)
  remains unmodeled — declared limit.

---

## 3. Layer 1 · Physics L1: Streamline-Curvature Through-Flow

### 3.1 Why this replaced TD3 (phase 11)

Until phase 11, "L1" was `turbo-design` 1.4.2 with three monkey-patches,
run with `num_streamlines=1` **because with more the library's
radial-equilibrium ODE collapsed its step and hung**. That is: L1 was
another meanline. It ran in a timeout-guarded subprocess because of those
hangs, and in practice it was dead in silence — the package does not
declare `requests`, so the import failed and the whole system ran at L0
believing it ran at L1 (found by the first test ever written for that
path; see docs/VALIDATION.md).

A fidelity ladder with one real rung is not a ladder. The residual deep
ensemble of layer 2 had no residual to learn (it short-circuits under
L0), and the affine L2 calibration was an API waiting for a user.

`scm_core.py` replaces it: a through-flow solver of our own, no external
dependency, running in-process.

### 3.2 The equation

Along a radial quasi-orthogonal, with $h_0$ the stagnation enthalpy, $s$
the entropy and $r_c$ the meridional radius of curvature:

$$C_m\frac{\partial C_m}{\partial r}
  = \frac{\partial h_0}{\partial r}
  - T\frac{\partial s}{\partial r}
  - \frac{C_u}{r}\frac{\partial (rC_u)}{\partial r}
  - \frac{C_m^2\cos\gamma}{r_c}$$

which follows from combining $dh = T\,ds + dp/\rho$ with the normal
equilibrium $\frac{1}{\rho}\frac{\partial p}{\partial r}
= \frac{C_u^2}{r} + \frac{C_m^2\cos\gamma}{r_c}$. Sign convention:
$1/r_c > 0$ when the centre of curvature sits at *smaller* radius, i.e.
$1/r_c = -r''/(1+r'^2)^{3/2}$ and $\cos\gamma = 1/\sqrt{1+r'^2}$.

With zero curvature and $C_u \propto r^n$ this integrates to the closed
form of §2.2 (`physics_core.vortex_cx`),

$$C_x^2(r) = C_{x,m}^2 - C_{u,m}^2\frac{n+1}{n}
  \left[\left(\frac{r}{r_m}\right)^{2n} - 1\right],$$

and the verification test T21 checks exactly that: the numerically
integrated ODE reproduces the analytic profile to $2\times10^{-16}$ for
free vortex and to $\sim2\times10^{-3}$ (trapezoid discretisation over 9
points) for $n = -0.5,\,0,\,+0.5$. The SCM is therefore an **extension**
of the phase-9.1 model, not a parallel one.

### 3.3 Closure and what actually gets solved

The blade is a **fixed object in space**: its metal angles $\beta_2(r)$
and $\alpha_{out}(r)$ come from the design vortex law and are frozen. The
solver then finds the field $(C_m, C_u, \rho)$ that simultaneously
satisfies

1. radial equilibrium at every station,
2. global continuity **and** continuity per streamtube (the streamlines
   are placed by mass fraction, not by radius),
3. the angle the blade imposes at each trailing edge,
   $C_{u2} = U - C_m\tan\beta_2$,
4. **Euler per streamline**, $\Delta h_0 = \omega\,(r_2C_{u2} -
   r_1C_{u1})$ — the work stops being uniform across the span, which is
   precisely what the meanline cannot see.

Stations sit at the leading and trailing edge of every row ($2n+1$ of
them for $n$ stages); the annulus is the same one layer 5a builds, so the
SCM solves the flow in the machine that gets manufactured.

Losses are resolved across the span: each streamline pays its own
diffusion, its own Reynolds and its own Mach through the same
`_row_losses` the meanline uses; the shock is only paid where the
relative Mach justifies it; tip clearance is concentrated in the outer
25% of span and the Koch & Smith end-wall debit in the wall bands,
instead of being smeared uniformly.

### 3.4 Numerics and honest failure

Outer loop: reposition streamlines by mass fraction, update meridional
curvature (lagged one iteration — implicit curvature is what
destabilises), repeat until the radii move less than $10^{-5}$
relative. Inner loop per station: the radial-equilibrium integral fixes
the *shape* of $C_m$, and the *level* comes from bisecting global
continuity — bracketed by the **sonic limit**, because mass flow only
grows with $C_m$ on the subsonic branch and a blind bisection lands on
the supersonic one (measured once: PR = 0.19 with an exit $T_0$ below
inlet).

Two failure modes are reported, never swallowed:

* **choked station** — not even at the sonic limit does the annulus pass
  the required mass flow;
* **profile limiter still active at convergence** — the numerical clamp
  on the $C_m$ profile (0.30–2.40 × its mean) is holding the answer
  rather than the physics; the field would be kneaded by the limiter, so
  it is refused.

Either raises `SCMDiverged`; `evaluate` then degrades to L0 **with the
reason in `source`**.

### 3.5 What it costs and what it buys

About 3 s per machine (4 stages, 9 streamlines, 15 outer iterations) —
comparable to the old subprocess but now solving something real, and
cacheable to `phys_cache.jsonl` like any L1 record. Measured against L0
on the reference θ:

| Vortex law | ΔPR vs L0 | Δη_poly | Span work spread |
|---|---|---|---|
| free ($n=-1$) | +0.13% | +0.34 pts | 6.9% |
| controlled ($n=-0.5$) | −6.5% | −0.01 pts | 5.6% |

The free-vortex row is **verification**: it is the case where the
meanline's assumptions are exact, so agreement is the expected result.
The controlled-vortex row is the **residual** the surrogate needs — and
it also exposes something the meanline never accounted for: applying the
vortex exponent to $C_{u1}$ and $C_{u2}$ alike makes the Euler work vary
as $r^{n+1}$, so for $n \neq -1$ the design is a non-uniform-work design
and the meanline only ever evaluated it at mid-span.

Not yet qualified against measurement — the same open gap as the
off-design map (docs/VALIDATION.md §5).

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
4. Blade-root centrifugal stress **with the fillet K_t** (§5.5):
   $K_t\,\sigma_{root} \le \sigma_y(T)$,
   $\sigma_{root} = k_{taper}\,\rho\,\omega^2\,(r_t^2 - r_h^2)/2$.
5. **Campbell margin** (§5.5): first-flap frequency at least ±10% away
   from the blade-passing engine orders of both neighbouring vane rows
   at design speed.

### 5.5 Blade Dynamics (phase 7, 2026-07-17)

* **Natural frequencies — rotating cantilever.** First flap mode of an
  Euler-Bernoulli clamped-free beam with Southwell's centrifugal
  stiffening (rotor rows only):

  $$f_1^2 = f_{1,0}^2 + S\,(N/60)^2, \qquad
    f_{1,0} = \frac{\lambda_1^2}{2\pi}\sqrt{\frac{E\,I_{min}}{\rho A h^4}},
    \quad \lambda_1^2 = 3.516,\; S = 1.6$$

  Section: the real hub polygon (`geometry_generator.
  polygon_section_props`, Green's theorem — principal moments) when
  available, else the rectangular equivalent with the layer-5a root t/c
  (same source as the printed part). First torsion
  $f_t = \tfrac{1}{4h}\sqrt{G J/(\rho I_p)}$ (thin-strip J) feeds the
  flutter screen. Verified against the analytical uniform cantilever to
  <2% in the suite.

* **Campbell diagram.** Hard constraint (5th g_struct component):
  relative margin $|f_1 - EO|/EO \ge 0.10$ against the blade-passing
  orders of the upstream and downstream vane rows (for rotor 1 the
  upstream row is the IGV) at design speed. Low engine orders k = 1..6
  are reported as metrics only — enforcing them all would block the
  space at preliminary level. Measured on LHS(500): structural
  feasibility unchanged (65% → 65%) — f₁ (10²–10³ Hz) sits far below
  blade passing (10³–10⁴ Hz) except for extreme geometries, which is
  exactly what the constraint should catch. `figures/campbell.png` plots
  f₁(N) per stage against the EO fans.

* **Root K_t and Goodman.** The root constraint pays the Peterson
  shoulder-fillet concentration factor (polynomial fit, Peterson via
  Shigley tab. A-15) evaluated with the printed fillet radius
  (`BLADE_FILLET_R_MM` = 2.0 — the SAME parameter layer 5c uses to
  build the root fillet by restricted morphological closing, phase 7.5):
  K_t ≈ 1.4–1.6 for the design space (sharp-corner 1+2√(t/r) ≈ 3–4 is
  over-conservative once the fillet is printed). LHS(500): blade-ratio
  p90 rose 0.55 → 0.93 with feasibility intact. **Goodman** is a
  reported metric (not g): remaining allowable vibratory stress
  $\sigma_{alt} = \sigma_e\,(1 - K_t\sigma_{root}/\sigma_{uts}(T))$
  with handbook σ_e per alloy — the standard preliminary-design
  deliverable when the excitation amplitude is unknown.

* **Flutter screen.** Tip reduced velocity $V^* = W/(b\,\omega_t)$
  against the classical bending-torsion threshold 1.4 (Armstrong &
  Stevenson 1960). Reported with a warning flag only — the L0s torsion
  model is too coarse to gate on.

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
  Phy-CC. `fixed_vars` pins any of the 12 design variables (the 10
  legacy ones plus `phi_slope`/`Rx_slope`, e.g. `{"n_stages": 5,
  "phi_slope": 0.0}`): every evaluation path goes through
  `fix_operating_point`, so the pinned dimension collapses for LHS,
  NSGA-II and acquisition alike (pinning both slopes to 0 reproduces
  the pre-phase-8 13-D space exactly).
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

**IGV and OGV (2026-07-16).** The meanline starts stage 1 with pre-swirl
α₁ ("implicit IGV", its loss charged to stage 1) and removes the exit
swirl with an OGV row whose loss and length it does count — but neither
row existed in the geometry, so the printed machine could not meet the
verified triangles (real axial inlet ≠ assumed α₁) nor the axial-exit
condition behind M_exit. The contract now carries both as top-level
`igv`/`ogv` rows (same section schema as stators, built by the same
free-vortex machinery): the IGV turns 0 → α₁(r) (accelerating row) ahead
of rotor 1, the OGV turns α₁(r) → 0 after the last stator. Carter's
deviation is generalized to follow the camber sign (metal always turns
PAST the target flow angle: χ₂ = β₂ − sgn(θ)·δ) — for the accelerating
IGV, Carter overpredicts δ by a few degrees (declared first-order
approximation). Layer 5c hangs the IGV from the first casing ring and
the OGV from the last one.

### 8.2 Layer 5c (C#/PicoGK — real machine construction)

Blade rows are built (since 2026-07-17) as **solid profile lofts**: the
contract's closed sections (real NACA-65/DCA thickness distribution, 60
CCW points, analytically free of self-intersection) become one watertight
mesh per blade — ruled side walls between densified ribs plus hub/tip
caps triangulated LE→TE (a zipper between the pressure and suction
chains, independent of start index and winding) — voxelized directly.
This preserves the LE radius, the max-thickness position and the
pressure/suction asymmetry that the previous camber-sheet recipe
destroyed (uniform thickness). The **Phy-CC v3 vane recipe** (camber
mid-surface thickened with `Voxels.voxMeshShell`) remains as the fallback
for legacy contracts without `points` and for rows thinner than 2 voxels,
where the shell clamps the thickness **and logs a warning** (the printed
blade is thicker than the verified design — previously a silent
distortion). Roots sink into their body; free ends retract clearance
plus the surface rounding.

Parts (split where real bolted joints sit — mid-gap between a stator TE
and the next rotor LE):

* **Shaft**: stubs past both ends, center bore.
* **RotorStage i** (bladed disc): disc web (0.5·axial chord, clamped
  4–14 mm — the §5 model) + hub flow-path shell segment (5 mm wall,
  incl. the spacer under its stator) + its rotor blades; bore slides
  over the shaft.
* **StatorRing i**: casing shell segment + stator vanes; end rings carry
  the bolted flanges (12 × M6 on the mid bolt circle); the bleed-stage
  ring carries the **bleed port** (radial boss + through hole); the
  first/last rings additionally carry the **IGV/OGV** vanes.
* Union views (Rotor / Casing / Assembly) for inspection.

A bill of materials (`bom.csv`) lists every part with quantities.

---

## 9. Complexity Analysis and Evaluation Budgets

Measured on a laptop CPU:

| Component | Cost |
|---|---|
| L0 meanline | ~0.5 ms/point |
| g_struct (per-stage discs, Thomas O(n)) | ~1 ms/point |
| L1 streamline-curvature through-flow (in-process) | ~3 s/point |
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
- Creveling, H. F. (1968). *Axial Flow Compressor Computer Program for
  Calculating Off-Design Performance*. NASA CR-72427.
- Boyer, K. M., & O'Brien, W. F. (2003). *An Improved Streamline
  Curvature Approach for Off-Design Analysis of Transonic Axial
  Compression Systems*. ASME J. Turbomach. 125(3) — off-design shock
  loss reference (not yet modeled).
- Koch, C. C., & Smith, L. H. (1976). *Loss Sources and Magnitudes in
  Axial-Flow Compressors*. ASME J. Eng. Power 98.
- Wassell, A. B. (1968). *Reynolds Number Effects in Axial Compressors*.
  ASME J. Eng. Power 90(2).
- Schäffler, A. (1980). *Experimental and Analytical Investigation of the
  Effects of Reynolds Number and Blade Surface Roughness on Multistage
  Axial Flow Compressors*. ASME J. Eng. Power 102 (79-GT-2).
- Denton, J. D. (1993). *Loss Mechanisms in Turbomachines*. ASME
  J. Turbomach. 115(4) (IGTI Scholar Lecture).
- Storer, J. A., & Cumpsty, N. A. (1991). *Tip Leakage Flow in Axial
  Compressors*. ASME J. Turbomach. 113(2); and (1994) *An Approximate
  Analysis and Prediction Method for Tip Clearance Loss in Axial
  Compressors*, 116(4).
- Sakulkaew, S., Tan, C. S., Donahoo, E., Cornelius, C., & Montgomery,
  M. (2013). *Compressor Efficiency Variation with Rotor Tip Gap from
  Vanishing to Large Clearance*. ASME J. Turbomach. 135(3).
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
- Southwell, R. V. (1921). *On the Free Transverse Vibrations of a
  Uniform Circular Disc Clamped at its Centre; and on the Effects of
  Rotation*. Proc. Roy. Soc. A 101 (Southwell coefficient).
- Peterson, R. E. *Stress Concentration Factors*, 2nd ed. (shoulder
  fillet fit via Shigley, *Mechanical Engineering Design*, tab. A-15).
- Armstrong, E. K., & Stevenson, R. M. (1960). *Some Practical Aspects
  of Compressor Blade Vibration*. J. Roy. Aero. Soc. 64.
- Deb, K., et al. (2002). *A Fast and Elitist Multiobjective Genetic
  Algorithm: NSGA-II*. IEEE Trans. Evol. Comput. 6(2).
- Lakshminarayanan, B., et al. (2017). *Simple and Scalable Predictive
  Uncertainty Estimation using Deep Ensembles*. NeurIPS.
- OpenOrion `turbodesigner` (MIT) — layer-5a/5c design reference.
- LEAP 71 `PicoGK` v2.2.0 — voxel geometry kernel (layer 5c).
