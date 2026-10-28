# Structural Parameters — Reference
**Last updated:** October 19, 2026 (Λ̄_r replaces α in Θ_e, Q1)

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
| x_m | `x_m` | Upper bound of the **sunk entry cost** distribution G(x) = (x/x_m)^ξ (draft `eq` at L887) | From distribution normalization (Stage 3) | ⚠️ **Not** χ_m. This row previously read "χ_m (x_m)", conflating two different bounds — corrected Sept 30, 2026 |
| χ_m | **`f_m`** | Upper bound of the **continuation cost** distribution F (draft `eq:F_chi`) | From ς and ψ (Stage 3) | **The draft's χ_m is the code's `f_m`, not `x_m`.** This is the bound Λ_r = F(Q/α) depends on, so getting it wrong silently breaks the reposting margin (R8) |
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
| ω_δ | `dest_end_frac` | 0.5 | Endogenous share of total destruction: ω_δ = (δ_e − δ)/δ_e | **cor(u, v)** (Beveridge curve), cor(δ_e, u), cor(δ_e, z) | **The Beveridge curve lever.** Higher ω_δ → more of what the data calls "separation" actually destroys product lines in the model → more negative u-v comovement from the δ mechanism. This is the key parameter for reconciling the s-shock problem (`lem:accounting` Part 1) |
| p_0 | `p_0` | 0.5 | Mass on continuous (Pareto) part of continuation cost distribution F | Interacts with ω_δ: shapes the nonlinearity of δ_e response | p_0 = 0: all firms face zero continuation cost → no endogenous exit margin. p_0 = 1: all firms draw from Pareto → maximum exit-margin responsiveness. Together with ψ, controls the curvature of the exit response. **Also caps the reposting margin: max position destruction from non-reposting = p_0** (see α below) |

### Reposting margin (in code since R8, Oct 14–19; go/no-go in [D11](decisions.md))

Adds **one** estimated parameter and **no new distribution** (reuses F). Full equations:
[`model_equations.md`](model_equations.md) Block 2.

**The estimated parameter is the steady-state reposting rate Λ̄_r, not α** (Q1, Oct 14).
`calibrate_shares` takes `Λr_target` and backs out α from Λ̄_r = F(Q̄/α); α is then a Θ_d
object like ϕ or f_m. Reason: α is not a natural unit (it ranges from 4e10 at Λ̄_r = 0.85 to
2e28 at 0.70, `R8_test_steady_state.txt`), while Λ̄_r is a share with a direct reading. The
interior requires 1 − p_0 < Λ̄_r < 1 (`steady_state.jl` errors otherwise).
`Λr_target = nothing` (the `TARGETS` default) is costless reposting, implemented by the
switch `rp = 0` in `run_solution_core.jl` (α = 0 would divide by zero, and Λ̄_r = 1 is a kink).

| Symbol | Code name | Default | Economic meaning | Key moments | Intuition |
|--------|-----------|---------|-----------------|-------------|-----------|
| Λ̄_r | `Λr_target` | `nothing` (costless); grid 0.95 / 0.85 / 0.70 in `R8_step7_baseline.jl` | Steady-state share of vacated positions that are reactivated | **cor(u, v)**, s-shock responses of u and v (Block B′), σ(v) | Lower Λ̄_r → more position destruction per separation → s shocks lower v → more negative cor(u,v). Step 7 (part6b shocks): cor(u,v) = +0.995 / +0.968 / +0.812 / +0.461 at costless / 0.95 / 0.85 / 0.70, against −0.804. Bounded below by 1 − p_0 |
| α | `α` (in PAR via `repost_par`) | backed out | Scale of the per-position reactivation cost α·χ, χ ~ F | — (Θ_d) | Reactivate iff Q ≥ α·χ, so Λ_r = F(Q/α). Threshold is the stock value Q, not the flow K (note §7.3b). Per-position draws |

**Shared F makes the margin act through its level.** With reactivation drawn from F,
d ln Λ_r/d ln Q ≈ 0.006 at the defaults: Λ_r barely moves over the cycle. The margin works
through Λ̄_r, not through cyclical reactivation (MS, Oct 14: no separate reactivation shape).

**Λ_r** (reposting rate, `Lambda_r`) is a *derived* object like δ_e, not a free parameter:
Λ_r = F(Q/α), a tracked control with its own equation. Symbol is capital Λ_r, **not** λ (which
is marginal utility of consumption).

### Shock processes (6 parameters)

> **Estimated S14 values (Oct 28, `shock_calibration.txt`), monthly:** ρ_z 0.9121, σ_z 0.0128; ρ_δ 0.5824, σ_δ 0.0663; ρ_s 0.8054, σ_s 0.0711. Method: [`session_handout_20261028.md`](session_handout_20261028.md) §2.
> The "Calibrated value" column shows the S14 estimate; the code (`SHOCKS_BASE`) still holds the
> pre-S13 value in parentheses until MS approves the swap. See [decisions.md](decisions.md) S13, S14, D14.

| Symbol | Code name | Calibrated value | Economic meaning | Key moments | Intuition |
|--------|-----------|-----------------|-----------------|-------------|-----------|
| ρ_z | `rho_z` | 0.9121 (code 0.902) | Monthly persistence of technology shock | ρ₁(u), ρ₁(v), shape of z→u IRF | Univariate AR(1) on HP-1600 log z, 1951–2019 (part6b step 1). Governs how long a productivity shock propagates |
| σ_z | `sigma_z` | 0.0128 (code 0.0092) | Monthly innovation SD of technology shock | σ(LP), σ(u) from z channel | AR(1) residual SD, 1951–2019. 39 % above the 1992–2019 value: the window adds pre-1985 volatility |
| ρ_δ | `rho_delta` | 0.5824 (code 0.592) | Monthly persistence of destruction shock | Shape of δ→u IRF, ρ₁(δ_e) | Own lag of the (z, δ) VAR(1), 1992Q3–2019Q4 (part6b step 2). Low persistence (quarterly ρ = 0.198) reflects noisy BED Deaths |
| σ_δ | `sigma_delta` | 0.0663 (code 0.0669) | Monthly innovation SD of destruction shock (orthogonal to z) | σ(δ_e), cor(δ_e, u), δ→u IRF scale | The only Cholesky-dependent number: stripped of contemporaneous z covariation, because the model generates z → χ^c → δ_e |
| ρ_s | `rho_s` | 0.8054 (code 0.8741) | Monthly persistence of match separation shock | ρ₁(u), ρ₁(v), ρ₁(τ) | Univariate AR(1) on s = (τ−δ_e)/(1−δ_e), 1951–2019, spliced from τ before 1992Q3 (part6b step 3) |
| σ_s | `sigma_s` | 0.0711 (code 0.0854) | Monthly innovation SD of match separation shock (not purged, D14) | **σ(u), σ(v), cor(u, v)**, σ(θ)/σ(LP) | 8.5× larger than the 0.010 placeholder. The s shock drives most of the business cycle volatility and the (wrong-sign) Beveridge correlation. See `lem:accounting` Part 1 |

---

---

## Reposting margin — background

**Status October 19, 2026.** In the draft (recruiter block, `eq:v_lom`, `eq:agg_repost_costs`,
`eq:rc`/`eq:gdp`, `def:equilibrium`) and, since R8, in `steady_state.jl` and
`run_solution_core.jl` (see the Θ_e table above: Λ̄_r is estimated, α backed out). Still
**not in §5.2**: no prior and no `tab:calib_targets` row for Λ̄_r (Q3/Q6). The September text
below is about α; read "α" as "the Λ̄_r it implies".

Full argument in
[`../Notes/delta_calibration_and_the_reposting_margin.md`](../Notes/delta_calibration_and_the_reposting_margin.md);
equations in [`model_equations.md`](model_equations.md) Block 2.

**One new estimated parameter, α** (reposting-cost scale). The responsiveness of Λ_r to Q is
*not* a separate parameter: it is the slope of F at Q/α, and F's shape (p_0, ψ, χ_m) is already
identified by the exit margin. So the margin asks the data for one number plus a consistency
check that F's shape fits both margins.

| Symbol | Code | Economic meaning | How disciplined |
|---|---|---|---|
| Λ̄_r = F(Q̄/α) | `Λr_target` | Steady-state reactivation rate. **Estimated since Q1 (Oct 14)**; α is backed out from it | Block B′ s-shock responses of u and v; cor(u,v) in Block M. See R2/Q4 for the identification check |

**Why Λ_r < 1 is the lever.** Λ_r = 1 is costless reposting (the pre-R8 model; `Λr_target = nothing`). Λ_r < 1 makes an s shock destroy
positions as well as matches, which is what moves cor(u,v). Raising δ̄_e does *not* fix
cor(u,v) — it improves it monotonically but by 0.108 across the whole defensible range against
a gap near 1.8, and the relation is convex with the leverage concentrated above the range
measurement allows. This does.

**Λ_r ∈ [1−p_0, 1], so p_0 caps the margin.** The atom at zero always reactivates, so at most
fraction p_0 of separations can be destroyed. p_0 therefore does double duty: it shapes the exit
margin *and* bounds the reposting margin. If cor(u,v) needs more than p_0 = 0.5 allows, the
estimator pushes p_0 → 1, which removes the zero-cost atom and changes the exit margin
materially.

Three modeling constraints, all easy to get wrong:

- **A constant reposting rate is not structural.** Reactivation falls when position values fall,
  i.e. in recessions, which is exactly when the Beveridge curve shifts. Fixing the rate holds it
  constant there and would be a reduced form dressed as a mechanism. Hence Λ_r = F(Q_t/α),
  endogenous and procyclical.
- **A homogeneous reposting fee does not produce Λ_r ∈ (0,1).** All positions being identical,
  the decision is all-or-nothing. An interior rate requires *dispersion* in the reactivation
  cost, exactly as the firm-level exit margin requires F to deliver an interior χ^c.
- **Reactivation draws from the same F as continuation costs**, so a single α is the only new
  parameter — and `p_0 = 0` shuts both margins at once, which is what leaves
  `prop:independence` intact with three conditions rather than four.

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
| Λ = F(χ^c) | **`surv_prob`** | Survival probability |
| ς | `zeta_ss` (ζ) | **Survival quantile**, the normalized Pareto-part factor (χ^c/χ_m)^ψ — *not* the survival probability. `steady_state.jl:708` sets `ζ_ss = (surv_prob − (1 − p_0))/p_0`, so Λ = (1 − p_0) + p_0·ς. ⚠️ This row previously mapped F(χ^c) → ζ, conflating the two; corrected Sept 30, 2026. They are numerically close only because Λ ≈ 0.9986 here, and diverge as p_0 falls |
| ψ | `psi` | Pareto shape of continuation costs |
| Xc/Yc | `Xc_Yc` | Fixed cost / value-added ratio |
