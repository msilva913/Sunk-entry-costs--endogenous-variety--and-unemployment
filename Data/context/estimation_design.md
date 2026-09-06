# Estimation Design — Conditional + Unconditional Moments
**Last updated:** September 5, 2026 · **Status:** design settled in the draft; **no code written**

The design lives in the draft's §5.2 prose (`sec:calib`). This file states it compactly,
records what is missing to run it, and points at the decisions that must be settled first.
Open decisions live in [`decisions.md`](decisions.md), not here.

---

## The design (settled — `Draft.tex` §5.2, `eq:posterior`)

Bayesian simulated method of moments over **two blocks**.

**Block M — unconditional moments.** m̂ = the σ(x), contemporaneous correlations, and ρ₁(x)
of `tab:smm_moments`, over seven series {u, v, s, f, δ, N^e, z}, in log-levels, HP λ=1,600.
Weight Ω_m: HAC, block bootstrap, block length 8 quarters (HP filtering induces serial
correlation in cyclical residuals at all lags, so this is needed regardless of λ).

**Block B — conditional (IRF) moments.** β̂ = the stacked δ→u and δ→v Bartik LP paths,
h = 0..20, so **K_β = 42**. Gaussian quasi-likelihood in the Christiano-Eichenbaum-Trabandt
(2016) sense, `eq:ql_irf`. Weight Ω_β: sampling covariance of the LP estimator, wild cluster
bootstrap, 50 state clusters.

**Combination rule — degrees-of-freedom normalization.** Each quadratic form is divided by
its own dimension:

```
log p(θ|data) ∝ log π(θ) − (1/2K_β)(β̂−β(θ))' Ω_β⁻¹ (β̂−β(θ))
                         − (1/2K_m)(m̂−m(θ))' Ω_m⁻¹ (m̂−m(θ))
```

Under correct specification each form is ≈ χ²(K) with mean K, so dividing by K makes the two
blocks contribute equally in expectation and stops the 42-dimensional IRF block from
mechanically swamping Block M (Dridi-Guay-Renault 2007). Robustness to rescaling Block M by
½ and 2 is promised in `app:weighting_robustness` — **that appendix does not exist yet**
(undefined reference in `Draft.log`).

**Filtering asymmetry — resolved on the model side, not the data side.** The blocks target
data objects built differently (HP-filtered second moments vs. level-difference LP
coefficients). Symmetry is imposed within each block:

- m(θ): simulate → aggregate to quarterly → HP(1600) → compute moments. Mirrors the data.
- β(θ): read straight off the linearized state-space impulse response, **unfiltered**.
  Mirrors the LP, which is estimated in level differences.

Justification for HP moments as a valid indirect-inference target: Canova-Ferroni (2011) —
applying the identical filter to data and simulation neutralizes filter-induced distortion.

**Why Block B is needed at all.** With three shocks, unconditional second moments cannot
separately discipline δ vs. s propagation: both raise u, and their vacancy responses — the
qualitative distinction — are collinear in second moments. The Bartik design is the only
clean exogenous variation in δ.

---

## Parameter blocks

**Estimated Θ_e** (12 as documented): b, x_v, ξ⁻¹, ω_δ, p_0, σ, and (ρ_x, σ_x) for
x ∈ {z, δ, s}. ⚠️ The contents depend on [D2](decisions.md) — under PATH B, ψ is pinned by
`dest_elast_target` rather than estimated, which changes the dimension of Θ_e.

**External:** r = 4%/yr, η_L = 0.6, ε = 4.3 (μ ≈ 1.30), τ = 3.1%/month, and δ_e — whose
value is [D1](decisions.md), unresolved.

**Dependent Θ_d:** recovered per draw by the four-stage calibration in
`steady_state.jl::calibrate_shares` → z, ϕ, f_e, A, s, ψ, χ_m, x_m, κ, δ.

---

## Determinacy — not a binding constraint

`Notes/Baseline_Blanchard_Kahn.md` (June 6, 2026) works out Blanchard-Kahn for the full
baseline (DS-CES, p_0 = 0): two jump variables (Ψ̃_t, Ĉ_t), leading-order upper-triangular
backward map, and BK satisfied under **ε > 1**. `Notes/AGS_Blanchard_Kahn.md` does the σ = 0
case. Two consequences for the sampler:

1. At the externally fixed ε = 4.3, determinacy holds comfortably — BK needs ε > 1 and
   Prop. 5 Part 2's condition `eq:gN_cond` needs ε > 2. **ε > 2 is *not* required for BK**;
   that is a Proposition 5 requirement only.
2. §12 of the note shows **endogenous exit only makes BK easier to satisfy** (Λ̄ < 1 shrinks
   both persistence coefficients, raising 1/a and 1/b). So draws with larger ω_δ or p_0 do
   not threaten determinacy.

The sampler should still check the solver's return code per draw and reject failures, but no
prior truncation on ε is needed for determinacy.

---

## What is missing to actually run it

1. **Ω_β does not exist.** `part5_lp.py` estimates each horizon in a separate regression and
   saves only a scalar clustered `se` per h. `eq:ql_irf` needs the full **42×42** covariance
   across horizons *and* across the two outcomes — LP coefficients are strongly correlated
   across h by construction. `part5_wcrb.py` must persist the matrix, not bootstrap SEs.
2. **Ω_m does not exist.** No block-bootstrap covariance of the empirical moments is computed
   anywhere; `observables_moments.py` produces point estimates only.
3. **The m(θ) pipeline is stale.** `Programs baseline/second_moments.jl` still uses
   λ = 100,000 and covers only {u, v, θ, labor_prod, ls, z, δ, w_R} — missing s, f, and N^e,
   and using the filter the paper explicitly rejects (S1 in `decisions.md`).
4. **No β(θ) extractor.** Nothing maps the state-space solution to a model IRF conforming to
   the LP object (1-SD δ shock, u and v in pp, quarterly, h = 0..20, levels difference vs.
   t−1). Blocked on [D3](decisions.md).
5. **No sampler.** `run_solution_core.jl:112` still reads `estimate = []  # filled in when
   SMM is wired up`; `priors = (;)`. `run_solution.jl` carries placeholder ρ_s = 0.90,
   σ_s = 0.010. The repo's `posterior_mode.mat` belongs to the **old Matlab/Dynare
   generation** of the model — do not reuse it.

---

## Suggested build order

Each step is testable on its own, and the first three are unblocked once D1–D3 are settled.

1. **D1, D2, D3 settled** → freeze one canonical target set in `data_and_files.md`.
2. **`part5_wcrb.py`** → `omega_beta.npy` (42×42) + a diagonal-vs-clustered-SE sanity check.
3. **`moments_bootstrap.py`** → `omega_m.npy` via block bootstrap, block length 8.
4. **Refresh `second_moments.jl`** → λ=1,600, all seven series, returning m(θ) in exactly
   the order of `tab:smm_moments`. Verify against the empirical table on simulated data at
   the calibrated point.
5. **`model_irf.jl`** → β(θ) in the LP's units and normalization (per D3).
6. **Objective + priors** → evaluate the log-posterior kernel once at the calibrated point;
   confirm both quadratic forms are O(K) before sampling.
7. **Sampler** → then `sec:posterior`, `app:weighting_robustness`, and the AGS counterfactual.
