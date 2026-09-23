# Structural Parameters — Reference
**Last updated:** September 21, 2026

Complete parameter set for the Gabrovski-Silva endogenous exit model. Every parameter
belongs to one of three groups: externally calibrated (Θ_ext), recovered per draw by the
calibration stages (Θ_d), or estimated via SMM (Θ_e). The table maps each to its economic
meaning and to the data moment(s) or model mechanism that disciplines it.

Cross-references: `steady_state.jl` for the calibration code, `estimation_design.md` for
the estimation design, `data_and_files.md` for numerical targets, `findings.md` for results.

---

## External parameters (Θ_ext) — fixed before estimation

| Symbol | Code name | Value | Economic meaning | Source / target |
|--------|-----------|-------|-----------------|----------------|
| r | `r_ann` | 4%/yr | Discount rate | Standard macro calibration |
| η_L | `eta_L` | 0.6 | Matching elasticity w.r.t. unemployment | Petrongolo & Pissarides (2001); Hosios condition |
| ε | `epsilon` | 4.3 | Elasticity of substitution across varieties (CES aggregator) | BGM (2012) / Compustat markup μ ≈ 1.30. Controls strength of **variety channel**: N→ρ→w^int→θ. Lower ε → stronger variety externality → more amplification from firm turnover |
| τ̄ | `sep` | 3.1%/mo | Total separation rate | Shimer (2005/2012), JOLTS. Sets the level of match turnover |
| δ̄_e | `dest_ann` | 3.2%/yr | Permanent establishment exit rate | BLS BED Deaths, emp-weighted, 1993–2019. D1 settled at 0.0320. Sets δ_e/τ ≈ 0.087. ⚠️ **A directly observed *lower bound*, not a point estimate** — the literature brackets δ at [3.2%, 10%]/yr. Candidate to move into Θ_e with a prior on that interval once [D11](decisions.md) is resolved. See [`../Notes/delta_calibration_and_the_reposting_margin.md`](../Notes/delta_calibration_and_the_reposting_margin.md) |
| f̄ | `f` | 0.41/mo | Gross job-finding rate | JOLTS. Pins matching efficiency A via q(θ) |
| q̄ | `q` | 0.80/mo | Gross vacancy-filling rate | JOLTS. Pins θ̄ via θ = f/q |
| N̄ | `N` | 1.0 | Steady-state firm mass | Normalization |
| w̄ | `w` | 1.0 | Steady-state wage | Normalization |
| X/Y | `X_Y` | 1.5% | Recruiting cost share of GDP | CET (2016). Pins the vacancy value Q̄ in Stage 4 |
| Xc/Y | `Xc_Y` | 10% | Fixed cost share of GDP (PATH A target) | Abraham et al. (2019). Pins π_s. **Superseded by `dest_elast_target`** under PATH B |

---

## Recovered parameters (Θ_d) — pinned per draw by the four-stage calibration

These are not free parameters. Given (Θ_ext, Θ_e), `calibrate_shares()` solves for them
in four sequential stages. They are functions of the estimated and external parameters.

| Symbol | Code name | Economic meaning | How pinned | Intuitive link to moments |
|--------|-----------|-----------------|-----------|--------------------------|
| z̄ | `z` | Technology level | z̄ = w̄/ρ̄ from wage normalization (Stage 1) | Sets the scale of output; normalized out |
| ϕ | `phi` | Worker bargaining power (Nash) | Residual from wage equation (Stage 1) given b, w̄, θ̄ | Higher ϕ → higher wages → lower firm surplus → less entry. Affects the level of θ̄ and amplification |
| A | `A` | Matching function scale | From q̄ = A·θ̄^(-η_L) (Stage 1) | Pins the job-finding and vacancy-filling rates at steady state |
| f_e | `f_e` | Sunk entry cost | Free entry: Q̄ = f_e · G'(ē) (Stage 4) | Higher f_e → fewer entrants → slower replacement of destroyed firms → more persistence |
| s̄ | `s` | Match separation rate (SS level) | From identity s̄ = (τ̄ − δ̄_e)/(1 − δ̄_e) (Stage 1) | Sets the base level of non-destructive separations |
| δ̄ | `delta` | Exogenous component of destruction | From ω_δ: δ̄ = δ̄_e·(1 − ω_δ) (Stage 3) | The primitive shock; endogenous exit amplifies it by 1/(1−ω_δ) |
| ψ | `psi` | Shape of continuation cost Pareto tail | PATH B: from dest_elast_target (Stage 3). PATH A: from Xc/Y | Controls curvature of the endogenous exit margin F(χ^c). Higher ψ → sharper exit threshold → more responsive δ_e to shocks |
| χ_m (x_m) | `x_m` | Upper bound of sunk entry cost distribution | From distribution normalization (Stage 3) | |
| f_m | `f_m` | Scale of continuation cost distribution | From ζ and ψ (Stage 3) | |
| κ | `kappa` | Per-period flow matching cost | From X/Y target (Stage 4) | Affects the flow cost of maintaining a vacancy |

---

## Estimated parameters (Θ_e) — identified by SMM (Blocks M + B)

These are the free parameters the posterior adjudicates. Each has a prior anchored at the
default value listed. The "Key moments" column says what data variation disciplines it —
i.e., what would change in the model if this parameter moved.

### Preference and technology

| Symbol | Code name | Default | Economic meaning | Key moments | Intuition |
|--------|-----------|---------|-----------------|-------------|-----------|
| b/w | `b_ratio` | 0.71 (PATH A) / 0.90 (PATH B) | Replacement ratio (unemployment benefit / wage) | σ(θ), σ(u), amplification σ(θ)/σ(LP) | **The Hagedorn-Manovskii lever.** Higher b → smaller surplus → larger % surplus fluctuations → more amplification. The single most powerful parameter for the Shimer puzzle. Also affects ϕ (recovered) |
| σ | `sigma` | 1.0 | Inverse IES (risk aversion / intertemporal substitution) | cor(u, LP), ρ₁(u), consumption dynamics | σ = 1 is log utility. Higher σ → more consumption smoothing → dampens hours response to shocks |
| ξ⁻¹ | `xi_inv` | 1.0 | Inverse elasticity of entry cost G(e) = e^(1+ξ)/(1+ξ) | σ(N_e), δ→u persistence tail, σ(v) | **Entry friction lever.** Higher ξ⁻¹ → steeper entry cost → slower replacement of destroyed firms → N stays depressed longer → more persistent IRFs. Near free entry (ξ⁻¹→0) attenuates δ→u (Comparison D). **Sept 21 sweep: ξ_inv ∈ [0.5, 8] shifts the δ→u peak only from h=1 to h=2 — persistence improves (unit-root limit) but the peak never becomes hump-shaped. Structural limitation, not parametric.** |
| x_v | `x_v` | 1.0 (PATH A) / 0.5 (PATH B) | Sunk share of total recruiting cost (x_v·f_e vs κ) | σ(v), cor(v, u) | Split between sunk and flow vacancy costs. Higher x_v → more of the cost is sunk → stronger option-value dynamics in vacancy creation |

### Endogenous exit margin

| Symbol | Code name | Default | Economic meaning | Key moments | Intuition |
|--------|-----------|---------|-----------------|-------------|-----------|
| ω_δ | `dest_end_frac` | 0.5 | Endogenous share of total destruction: ω_δ = (δ_e − δ)/δ_e | **cor(u, v)** (Beveridge curve), cor(δ_e, u), cor(δ_e, z) | **The Beveridge curve lever.** Higher ω_δ → more of what the data calls "separation" actually destroys product lines in the model → more negative u-v comovement from the δ mechanism. This is the key parameter for reconciling the s-shock problem (Prop 5 Part 1) |
| p_0 | `p_0` | 0.5 | Mass on continuous (Pareto) part of continuation cost distribution F | Interacts with ω_δ: shapes the nonlinearity of δ_e response | p_0 = 0: all firms face zero continuation cost → no endogenous exit margin. p_0 = 1: all firms draw from Pareto → maximum exit-margin responsiveness. Together with ψ, controls the curvature of the exit response |

### Shock processes (6 parameters)

| Symbol | Code name | Calibrated value | Economic meaning | Key moments | Intuition |
|--------|-----------|-----------------|-----------------|-------------|-----------|
| ρ_z | `rho_z` | 0.902 | Monthly persistence of technology shock | ρ₁(u), ρ₁(v), shape of z→u IRF | From VAR(1) on HP-1600 log(z), part6b. Governs how long a productivity shock propagates |
| σ_z | `sigma_z` | 0.0092 | Monthly innovation SD of technology shock | σ(LP), σ(u) from z channel | Reduced-form residual SD (z ordered first in Cholesky) |
| ρ_δ | `rho_delta` | 0.592 | Monthly persistence of destruction shock | Shape of δ→u IRF, ρ₁(δ_e) | From VAR(1), Cholesky z-first. Low persistence (quarterly ρ = 0.207) reflects noisy BED Deaths |
| σ_δ | `sigma_delta` | 0.0669 | Monthly innovation SD of destruction shock (orthogonal to z) | σ(δ_e), cor(δ_e, u), δ→u IRF scale | The only Cholesky-dependent number: stripped of contemporaneous z covariation |
| ρ_s | `rho_s` | 0.874 | Monthly persistence of match separation shock | ρ₁(u), ρ₁(v), ρ₁(τ) | From AR(1) on s = (τ−δ_e)/(1−δ_e), part6b. The dominant separation shock |
| σ_s | `sigma_s` | 0.0854 | Monthly innovation SD of match separation shock | **σ(u), σ(v), cor(u, v)**, σ(θ)/σ(LP) | 8.5× larger than the 0.010 placeholder. The s shock drives most of the business cycle volatility and the (wrong-sign) Beveridge correlation. See Prop 5 Part 1 |

---

---

## Prospective parameters — not in the model yet

Conditional on [D11](decisions.md). Full argument in
[`../Notes/delta_calibration_and_the_reposting_margin.md`](../Notes/delta_calibration_and_the_reposting_margin.md).

| Symbol | Economic meaning | How disciplined | Intuition |
|---|---|---|---|
| λ̄ = G(Q̄) | Steady-state **reposting rate**: fraction of match separations after which the firm reposts the vacancy | LD→vacancy IRF (sign and magnitude); cor(u,v) in Block M | λ̄ = 1 is the current model. λ̄ < 1 makes s destroy positions as well as matches, which is the lever on the Beveridge curve. **The binding one** — raising δ̄_e does not fix cor(u,v), this does |
| elasticity of λ_t to Q_t | Responsiveness of reposting to the value of a position | Shape of the LD→vacancy IRF across horizons | Makes reposting procyclical and position destruction countercyclical. **Cannot be calibrated at all** — the object exists only inside the model |

Two modeling constraints worth recording, because both are easy to get wrong:

- A **constant λ is not structural**. Reposting falls when position values fall, i.e. in
  recessions, which is exactly when the Beveridge curve shifts. A fixed λ holds it constant
  there and would be a reduced form dressed as a mechanism.
- A **homogeneous reposting fee does not produce λ ∈ (0,1)**. All positions being identical,
  the decision is all-or-nothing. An interior rate requires *dispersion* in the reposting
  cost, exactly as the firm-level exit margin requires the distribution F to deliver an
  interior χ^c.

## Parameter-to-moment intuition map

The two estimation blocks identify different subsets of Θ_e:

**Block M (unconditional moments)** — σ(x), cor(x,y), ρ₁(x) for {u, v, s, f, δ, N^e, z}:
- b/w: main driver of amplification (σ(θ)/σ(LP))
- ω_δ: main driver of cor(u, v) sign
- ξ⁻¹: affects σ(N_e) and the speed of recovery after shocks
- σ: affects cor(u, LP) and consumption dynamics
- Shock (ρ, σ): each shock's own volatility and persistence in the moment table

**Block B (conditional IRF moments)** — δ→u and δ→v paths, h=0..20:
- ξ⁻¹: peak horizon of δ→u (how quickly entry replaces destroyed firms)
- ω_δ: IRF scale (how much destruction per unit of the δ shock)
- p_0 + ψ: nonlinearity of the exit response; shape of the IRF
- ε (external, but matters): strength of the N-variety channel that generates persistence

**Parameters that pull in opposite directions on different moments** (the identification tension):
- ω_δ: higher → better Beveridge curve, but also higher δ_e volatility which may overshoot σ(δ_e)
- b/w: higher → more amplification, but also makes wages less responsive to θ, affecting IRF shapes
- ξ⁻¹: higher → more persistent δ→u, but also less entry → lower σ(N_e). **Cannot shift peak horizon** (Sept 21 sweep)
- ε: lower → stronger variety channel, but ε < 3 hits BK violation (determinacy ceiling for the persistence channel)

---

## Notation concordance

| Draft / equations | Code | This file |
|---|---|---|
| δ̄_e | `cal.δ` (exog) + endogenous part | Total SS destruction rate |
| ω_δ | `dest_end_frac` | Endogenous share |
| ξ | 1/`xi_inv` | Entry cost elasticity |
| χ^c | `x_c` (computed in SS) | Continuation cost cutoff |
| F(χ^c) | `zeta_ss` (ζ) | Survival probability |
| ψ | `psi` | Pareto shape of continuation costs |
| Xc/Yc | `Xc_Yc` | Fixed cost / value-added ratio |
