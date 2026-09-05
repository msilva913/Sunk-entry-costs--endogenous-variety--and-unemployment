# Estimation Design — Conditional + Unconditional Moments
**Last updated:** September 5, 2026 · **Status:** design settled in draft; **code not written**

This file records the resolution of the question "how do conditional (IRF) moments and
unconditional (second) moments get combined?", which is currently documented only in the
draft prose of Section~\ref{sec:calib}. It also records exactly what is missing to run it.

## The design (settled — Draft.tex §5.2, `eq:posterior`)

Bayesian simulated method of moments with **two blocks**:

**Block M — unconditional moments.** $\hat m$ = the $\sigma(x)$, contemporaneous
correlations, and $\rho_1(x)$ of Table `tab:smm_moments` (7 series: $u,v,s,f,\delta,N^e,z$),
log-levels, HP $\lambda=1{,}600$. Weight $\Omega_m$: HAC, block bootstrap, block length 8
quarters (HP filtering induces serial correlation in cyclical residuals at all lags).

**Block B — conditional (IRF) moments.** $\hat\beta$ = the stacked $\delta\to u$ and
$\delta\to v$ Bartik LP paths, $h=0,\dots,20$ (so $K_\beta = 42$). Gaussian
quasi-likelihood in the CET (Christiano-Eichenbaum-Trabandt 2016) sense, `eq:ql_irf`.
Weight $\Omega_\beta$: sampling covariance of the LP estimator, **wild cluster bootstrap,
50 state clusters**.

**Combination rule — degrees-of-freedom normalization.** Each quadratic form is divided
by its own dimension:

    log p(θ|data) ∝ log π(θ) − (1/2K_β)(β̂−β(θ))'Ω_β^{-1}(β̂−β(θ))
                              − (1/2K_m)(m̂−m(θ))'Ω_m^{-1}(m̂−m(θ))

Rationale: under correct specification each form is ≈ χ²(K) with mean K, so dividing by K
makes the two blocks contribute equally in expectation and stops the 42-dimensional IRF
block from mechanically swamping the moment block (Dridi-Guay-Renault 2007).
Robustness to rescaling Block M by ½ and 2 is promised in `app:weighting_robustness` —
**that appendix does not exist yet** (undefined reference in Draft.log).

**Filtering asymmetry — resolved on the model side, not the data side.** The two blocks
target data objects built differently (HP-filtered second moments vs. level-difference LP
coefficients). The fix is symmetry of treatment within each block:
- $m(\theta)$: simulate → aggregate to quarterly → HP(1600) → compute moments (mirrors data).
- $\beta(\theta)$: read straight off the linearized state-space IRF, **no filtering**
  (mirrors the LP, which is estimated in levels differences).
Justification for HP moments as a valid indirect-inference target: Canova-Ferroni (2011) —
applying the identical filter to data and simulation neutralizes filter-induced distortion.

**Why the IRF block is needed at all (identification argument, already written):** with
three shocks, unconditional second moments cannot separately discipline $\delta$ vs. $s$
propagation — both raise $u$, and their vacancy responses are collinear in second moments.
The Bartik design is the only clean exogenous variation in $\delta$.

## Parameter blocks

- **Estimated $\Theta_e$:** $b$, $x_v$, $\xi^{-1}$, $\omega_\delta$, $p_0$, $\sigma$, and
  $(\rho_x,\sigma_x)$ for $x\in\{z,\delta,s\}$ → **12 parameters**.
- **External:** $r=4\%$/yr, $\eta_L=0.6$, $\varepsilon=4.3$ ($\mu\approx1.30$),
  $\tau=3.1\%$/month, $\delta_e$ (see open issue below).
- **Dependent $\Theta_d$:** recovered per draw by the 4-stage calibration in
  `steady_state.jl::calibrate_shares` → $z,\phi,f_e,A,s,\psi,\chi_m,x_m,\kappa,\delta$.

## What is missing to actually run it

1. **$\Omega_\beta$ does not exist.** `part5_lp.py` estimates each horizon in a separate
   regression and saves only a scalar clustered `se` per $h$. `eq:ql_irf` needs the full
   $42\times42$ covariance across horizons *and* across the two outcomes — LP coefficients
   are strongly correlated across $h$ by construction. No script computes or stores it.
   → `part5_wcrb.py` must save the full covariance, not just bootstrap SEs.
2. **$\Omega_m$ does not exist.** No block-bootstrap covariance of the empirical moments is
   computed anywhere. `observables_moments.py` produces point estimates only.
3. **$m(\theta)$ pipeline is stale.** `Programs baseline/second_moments.jl` still uses
   $\lambda=100{,}000$ and only $\{u,v,\theta,\text{labor\_prod},ls,z,\delta,w_R\}$ — it is
   missing $s$, $f$, and $N^e$, and uses the filter the paper explicitly rejects.
4. **$\beta(\theta)$ extraction does not exist.** No routine maps the state-space solution
   to a model IRF vector conforming to the LP object (1-SD $\delta$ shock, $u$ and $v$ in
   percentage points, quarterly, $h=0..20$, levels difference vs. $t-1$).
5. **No sampler.** `run_solution_core.jl:112` still reads `estimate = []  # filled in when
   SMM is wired up`; `priors = (;)`. `run_solution.jl` carries placeholder $\rho_s=0.90$,
   $\sigma_s=0.010$. The only `posterior_mode.mat` in the repo belongs to the **old
   Matlab/Dynare generation** of the model and must not be reused.

## Open specification decisions

- **Units/scaling of $\beta(\theta)$ vs $\hat\beta$.** $\hat\beta$ is pp response per
  1-SD of the *Bartik instrument* ($\hat\sigma=17.4$ pp), not per 1-SD of the structural
  $\delta$ shock. A scale factor must be either calibrated externally or estimated as a
  nuisance parameter; this is not addressed in the draft yet and it directly determines
  whether the IRF block identifies $\sigma_\delta$ or only the *shape* of the response.
  Cleanest fix: normalize both $\hat\beta$ and $\beta(\theta)$ by their own $h$-peak, or
  add a free scale parameter with a tight prior. **Decide before writing the sampler.**
- **Time FE in the LP vs. aggregate model IRF.** The LP has time fixed effects, so
  $\hat\beta$ is a *relative* (cross-state) response with the national component absorbed;
  the model IRF is an aggregate response. Same issue as above in a different guise — needs
  an explicit statement in §5.2 that the two are comparable up to scale.
- **Whether $s$-shock IRFs enter Block B.** Currently only $\delta\to u$ and $\delta\to v$.
  LD is not separately identified (r(δ,LD)=0.334; joint LP kills the LD coefficient), so
  the current answer is no — but that means $\rho_s,\sigma_s$ lean entirely on Block M.

See [`pending_tasks.md`](pending_tasks.md) for the ordered task list.
