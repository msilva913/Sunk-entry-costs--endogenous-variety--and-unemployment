# Pending Tasks
**Last updated:** September 5, 2026 (status re-verified against repo after 3-month pause;
previous update June 5, 2026)

## ⛔ Blockers for Section 5 (do these first)

**B1. Resolve the $\delta_e$ calibration contradiction.** ⛔ HIGH
The paper's stated departure from Coles-Kelishomi is that $\delta$ is set from BED Deaths
rather than from the aggregate separation rate — but the calibration does not implement it.
- Draft §5.1 says $\bar\delta \approx 1.1\%$/quarter (BED Deaths) → **annual 4.25%**,
  monthly 0.361%, $\delta_e/\tau = 0.116$.
- Draft §5.2 external block + `tab:calib_targets` say $\delta_e^{ann}=7.54\%$, monthly
  0.651%, from Jaimovich-Siu 21% × $\tau=3.1\%$ → $\delta_e/\tau = 0.210$.
- `Programs baseline/steady_state.jl:566` has `dest_ann = 0.0754` with the comment
  "[BED Deaths, emp-weighted]" — the value is the Jaimovich number, the label is BED.
This is a 1.8× difference in the steady-state level of the paper's central shock. It moves
the entry-cushion coefficient $\bar\delta_e/(r+\bar\delta_e)$ in Proposition 5 Part 3, the
$\delta_e/\tau\approx0.21$ quantitative claim at `Draft.tex:1481`, and Comparison D.
Note the direction: the BED value (0.116) makes Part 3's sufficient condition *easier* to
satisfy, so fixing this strengthens the proposition. Decide one source, then propagate to
§5.1, §5.2, `tab:calib_targets`, `steady_state.jl`, and every mechanism figure.

**B2. Wild cluster bootstrap must return the full $\Omega_\beta$.** ⛔ HIGH — `part5_wcrb.py` ❌
`part5_lp.py` runs each horizon as a separate regression and saves only a scalar clustered
`se` per $h$. `eq:ql_irf` needs the $42\times42$ covariance across $h=0..20$ **and** across
the two outcomes. Bootstrap SEs alone are not sufficient. See
[`estimation_design.md`](estimation_design.md).

**B3. Estimation code does not exist.** ⛔ HIGH
`run_solution_core.jl:112` still reads `estimate = []  # filled in when SMM is wired up`;
`priors = (;)`; `run_solution.jl` carries placeholder $\rho_s=0.90$, $\sigma_s=0.010$.
Missing: $\Omega_m$ (block bootstrap), a refreshed $m(\theta)$ simulator at $\lambda=1{,}600$
covering $\{u,v,s,f,\delta,N^e,z\}$, a $\beta(\theta)$ extractor matching the LP object, and
a sampler. The repo's `posterior_mode.mat` is from the **old Matlab/Dynare model** — do not
reuse. Also unresolved: the scale/units mapping between $\hat\beta$ (pp per 1-SD of a Bartik
instrument, time-FE-absorbed) and $\beta(\theta)$ (aggregate model IRF).

**B4. Draft cross-references are broken.** ⚠️ LOW effort, do opportunistically
From `Draft.log`: undefined `sec:conclusion`, `app:robustness`, `app:weighting_robustness`,
`eq:labor_C_N`; multiply-defined `eq:profit_share`. `app:weighting_robustness` is promised
by the §5.2 weighting discussion and needs to be written, not just relabelled.

## Status corrections found on re-verification (Sept 5, 2026)

- Section 4.3 **severity-placebo paragraph is written** (`Draft.tex:2176-2196`), contrary to
  the June entry below. Still missing from 4.3: the $\delta$-LD correlation paragraph and an
  in-text joint-LP sentence (currently only a bare pointer to `app:diagnostics`).
- `part7b_sloos.py` **exists and has been run**, but its outputs
  (`lp_irf_delta_sloos_*.csv`, `instr_sd = 11.16`) predate the May 14 switch to BED Deaths +
  v2 and are **stale**. Same for `lp_irf_delta_gfc.csv`. Both need a re-run against the
  current instrument ($\hat\sigma = 17.4$ pp) before being cited.
- `part7d_ld_gfc.py` does not exist; the GFC diagnostic has never been run as specified.

## Model Mechanism Tasks

See `Inspecting_mechanism_setup.md` for full design rationale. All three comparisons use
`b_ratio=0.9, x_v=0.5` as shared calibration settings (activated May 21). Figures are
4×2 portrait grids (rows: u/v, θ/w_int, N/N_e, asset values or exit diagnostic).

### Comparison A — Level of δ (CK timing)   ✅ COMPLETE
- Files: `run_solution_all_delta.jl`, `plot_mechanism_comparison.jl`
- Output: `irf_all_delta.jls`, `mechanism_A_z_shock.pdf`, `mechanism_A_delta_shock.pdf`
- Key finding: high-δ (δ_e=τ) shows larger N response to z shock (LOM multiplier 4.8×)
  and larger Q response to δ shock (short-duration asset amplification); baseline
  (low δ_e) shows more persistence in δ shock due to slow N recovery.

### Comparison B — Endogenous vs. Exogenous Exit   ✅ COMPLETE
- Files: `run_solution_endog_exit.jl`, `plot_endog_exit_comparison.jl`
- Output: `irf_endog_exit.jls`, `mechanism_B_z_shock.pdf`, `mechanism_B_delta_shock.pdf`
- Key finding: δ shock response is SMALLER under endog exit, primarily because
  dest_end_frac=0.5 halves δbar (δbar_endog=0.00326 vs δbar_exog=0.00651); the
  dynamic x_c amplification mechanism is dormant at calibrated ψ≈0.014.
- Calibrated ψ≈0.014 is far below Broer et al. (2025) preferred ψ=1; report implied
  dest_el as outcome alongside ψ for external validation.

### Comparison B (extension) — Re-calibrate targeting ψ ≈ 1.0   ❌ PENDING
- Motivation: at ψ≈0.014, the x_c continuation-cost margin is nearly dormant; the
  current Comparison B illustrates exposure effect (δbar halved) not dynamic buffering.
  Targeting ψ≈1.0 directly (as in Broer et al. IER 2025) would activate the x_c channel.
- Implementation: add `ψ_target=1.0` to `calibrate_shares` targets in `run_solution_endog_exit.jl`
  and re-run; or add a third variant to the existing comparison.
- Decision required: is the current exposure-effect result sufficient for the paper, or
  does the mechanism section need to show the dynamic x_c story?

### Comparison C — Variety Effects (N → ρ → w_int → JCC)   ✅ COMPLETE
- Files: `run_solution_no_variety.jl`, `plot_no_variety_comparison.jl`
- Output: `irf_no_variety.jls`, `mechanism_C_z_shock.pdf`, `mechanism_C_delta_shock.pdf`
- Implementation: set ρ≡1 by replacing f[13] in model equations post-gen; dedicated
  `calibrate_shares_no_variety()` with ρ=1 hardwired in Stage 4
- Targets: dest_end_frac=0.5, p_0=0.5, dest_elast_target=5.0, b_ratio=0.9, x_v=0.5
- Color: darkorange dotted for no-variety arm
- Draft: Section 4.3 Comparison C text + figure environments added (mechanism_C_z_shock, mechanism_C_delta_shock)
- Notation: survival quantile renamed \varsigma (was \zeta) throughout draft to avoid clash with variety taste parameter

### Comparison D — Role of Entry Elasticity (ξ_inv)   ✅ COMPLETE
- Files: `run_solution_entry_elasticity.jl`, `plot_xi_inv_comparison.jl`
- Output: `irf_xi_inv.jls`, `mechanism_D_z_shock.pdf`, `mechanism_D_delta_shock.pdf`
- Comparison: baseline ξ_inv=1.0 vs. near-free-entry ξ_inv=0.1 (both recalibrated; ϕ, κ, x_m adjust)
- Key findings: lower ξ_inv amplifies z-shock amplification (larger entry collapse via G(Q) channel)
  and attenuates δ-shock unemployment response (stronger entry cushion from duration shortening).
  At baseline δ_e/τ≈21%, negative u-v comovement for δ shocks is robust across full ξ range
  including free entry — empirically calibrated δ̄_e rules out CK Beveridge curve shift.
- New appendix section: app:additional_comparisons ("Impulse responses: additional comparisons"),
  subsection app:comparison_D. Placeholder ready for future comparisons.
- D5 extended: eq:K_decomp, eq:Q_ll_delta, eq:e_ll_delta (full K decomposition for δ shocks;
  corrected transmission chain showing K↑ → Q↑ → e↑ for δ shocks, not K↓)

### Free-entry comparison (steady_state_checks.jl)   ❌ PENDING (lower priority)
- Goal: side-by-side baseline (ξ_inv=1) vs. free-entry (ξ_inv=0) in steady_state_checks.jl
- Purpose: show ξ_inv has no handle on z-shock amplification (ϕ always recalibrates)
- Identification swap: DROP x_v; ADD b/w_int target (=baseline 0.8586); INHERIT κ; ϕ residual
- Functions written: `calibrate_shares_free_entry(targets, κ_fixed, b_w_int_target)` and
  `steady_state_free_entry(para)` — both in steady_state.jl
- Current failure: θ=1.80 found instead of ~0.51; JCC has two branches at ϕ=0.624;
  solver init log(0.51) may land near wrong branch
- Most promising fix: use bracketed solver in steady_state_free_entry,
  e.g. find_zero(jcc_res, (log(0.1), log(0.6))) to force low-θ branch
- Key insight: ϕ genuinely differs at free entry (0.624 vs 0.670) because K halves;
  this is correct behavior but shifts the JCC curve

## Empirical Code Tasks (priority order)

1. **Wild cluster bootstrap** — `part5_wcrb.py` ❌ — **see blocker B2 above**
   - n=50 state clusters is at the lower bound of asymptotic SE reliability
   - Use `wildboottest` or `linearmodels` bootstrap; supplement all main LP tables
   - NOT just SEs: must persist the full cross-horizon, cross-outcome covariance
     (42×42 for h=0..20 × {u,v}). That matrix is Ω_β in `eq:ql_irf` and nothing
     in the pipeline currently produces it.

2. **GFC diagnostic for LD unemployment** — `part7d_ld_gfc.py` ❌
   - LD unemployment IRF rises monotonically through h=20, inconsistent with ρ_LD≈0.3
   - Apply GFC_{t+h} outcome-quarter control dummy to test whether GFC drives the pattern

3. **SLOOS C&I interaction** — `part7b_sloos.py` ⚠️ WRITTEN, RUN, but OUTPUT STALE
   - Replace NFCI_risk with Senior Loan Officer Opinion Survey C&I net tightening (FRED)
   - Tests whether NFCI result reflects credit supply vs. recession severity
   - Directly relevant to industry-credit-conditions hypothesis for r(δ,LD)=0.334
   - Existing `lp_irf_delta_sloos_*.csv` were produced 2026-04-09 with instr_sd=11.16,
     i.e. the pre-BED-Deaths instrument. Re-run before citing.

4. **Pre-GFC sample restriction** ❌
   - Add `max_qt="2007Q4"` option to `run_lp()` in `part5_lp.py`
   - Tests whether monotonically rising s/LD IRF is GFC-driven

5. **Refresh BED cache to 2024Q4** ❌
   - Run `refresh_bed_cache.py` locally (blocked in sandbox by BLS API rate limit)
   - Current `ext_2024` sample in part11 truncates at 2021Q4

6. **build_report_html.py** ❌
   - Pandoc-based Markdown → HTML with base64-embedded figures
   - Requires `conda install -c conda-forge pandoc` locally

## Paper Writing Tasks (priority order)

1. ~~**Section 3.10 — δ vs. s asymmetry mechanism**~~ ✅ DONE — superseded by `prop:ds_asymmetry` (Proposition 5, complete June 5, 2026). Formal 3-part proposition + lem:vpre + full proof in app:proof_ds. See findings.md for full structure.

2. **Section 4.3 — Results** (mostly drafted) ⚠️
   - Numbers in the draft are the BED-Deaths baseline: δ→u plateau +1.71 pp at h=17–20
     (sig. from h=0); δ→v trough −0.58 pp at h=18 (sig. from h=2). The +1.50/−0.44
     figures listed here in June came from the superseded closings×π instrument.
   - ✅ severity placebo paragraph is written (`Draft.tex:2176–2196`)
   - Still needs: δ–LD correlation paragraph; an in-text joint-LP sentence (the section
     currently ends with a bare, malformed `\ref{app:diagnostics}`)

3. **Section 5 — Quantitative analysis** ⚠️ design written, results absent
   - ✅ Calibration written as the four-stage sequential approach (`tab:calib_targets`)
   - ✅ BSMM design + two-block weighting written (`eq:posterior`) — see
     [`estimation_design.md`](estimation_design.md)
   - ❌ §5.4 `sec:posterior` is an empty stub (three TODO comments only)
   - ❌ Blocked on B1 (which δ_e?), B2 (Ω_β), B3 (no estimation code exists)
   - ❌ `app:weighting_robustness` promised in §5.2 but never written
   - ❌ Two-model counterfactual (baseline vs. AGS) — not started

4. **Section 6 — Conclusion** ❌

## Completed Code Tasks

- ✅ `run_solution_core.jl` full audit (May 20): naming r/ρ, κ per match, u LOM total v_t, BFE at f[4], SS_symbolics x_c/C/Y corrected
- ✅ `part6b_shock_persistence_cyclical.py`: bivariate VAR(1) HP-1600 log(z)/log(δ), Cholesky z-first → ρ_z^m=0.902, ρ_δ^m=0.592
- ✅ `run_solution_entry_elasticity.jl` + `plot_xi_inv_comparison.jl` (May 31): Comparison D complete; see above
- ✅ `run_solution_all_delta.jl` (May 20–21): baseline vs. high-δ (δ_e=τ, pure exogenous); b_ratio=0.9, x_v=0.5; serializes `irf_all_delta.jls`
- ✅ `plot_mechanism_comparison.jl` (May 21): 4×2 portrait (u/v, θ/w_int, N/N_e, K/Q); K+Q row reveals short-duration asset amplification; interpreted comments added
- ✅ `run_solution_endog_exit.jl` (May 21): full baseline (endog exit) vs. exog exit; b_ratio=0.9, x_v=0.5; appends exit_flow=δ_e+N; serializes `irf_endog_exit.jls`
- ✅ `plot_endog_exit_comparison.jl` (May 21): 4×2 portrait with exit_flow replacing Q; comprehensive IRF interpretation comments added (δbar halving dominant, ψ≈0.014 dormant x_c)
- ✅ `Inspecting_mechanism_setup.md` (May 21): Comparison B results and two-factor interpretation documented; plot design table and implementation order updated

## Completed Writing Tasks

- ✅ Appendix D5 paragraph 9: K decomposition for δ shocks (eq:K_decomp, eq:Q_ll_delta, eq:e_ll_delta);
  corrected transmission chains; sign asymmetry ê^δ>0, ê^z<0 explained via duration shortening
- ✅ New appendix section: "Impulse responses: additional comparisons" (app:additional_comparisons);
  Comparison D with CK connection, ξ governs cushioning strength, δ_e/τ robustness result

- ✅ `prop:bgm_nest` (Nesting LR-BGM) + LR-BGM definition paragraph — Section 3 after Remark 1
- ✅ Proof of `prop:bgm_nest` — Appendix C (three-block: LOM eq:N_lom_eq, Euler eq:N_euler_eq, GDP eq:gdp)
- ✅ Section 4.2 LP specification: generic y_{s,t} formulation with label eq:lp
- ✅ Figure placement: fig:delta_distribution in §4.1; fig:delta_uv_irf in §4.2
- ✅ Vacancy convention corrected throughout: pp levels difference, not log (Principle 15)
- ✅ lp_irf_delta_uv_baseline.png now generated by part5_lp.py section 8g (_plt_diagnostics)
- ✅ instrument_delta_distribution.png generated by part3_resid_instruments.py standalone block
- ✅ Section 5.1 CK comparison compressed
- ✅ Section 5.2: filter switched to λ=1,600; Table 4 updated; moments paragraph rewritten; QL justification added; τ→s notation fixed
- ✅ Prop 6 (DS-CES equilibria) restated for DS-CES; detailed 5-step proof in new appendix section
- ✅ Prop 5 (prop:ds_asymmetry, δ-s asymmetry) updated with formal lem:vpre + 3-part proof in app:proof_ds (June 5, 2026); draft renamed to `Draft/Draft.tex`
- ✅ Appendix B new subsection `app:filter_robustness`: 3-filter dynamic correlation table with bold HP-100k outlier cells
- ✅ `observables.py`: BED Deaths via BLS API (dataclass 08), 1994 CPS correction, BAWBATOTALSAUS, plt.show() fix
- ✅ `labor_market_dyn_corr_by_filter.py`: standalone script producing dynamic correlation plots and CSV by filter

## Specification Decisions Still Open

- **Primary LP for paper**: joint (part5_joint_lp.py) or separate (part5_lp.py)?
  Current recommendation: separate LP as primary; joint as robustness for δ. δ vacancy IRF stable across both.
- **v3 residualization** (+ monetary policy sensitivity φ_j × ΔMP_t): coded but not yet run; appears in draft as appendix robustness.
- **Industry credit control** (Route B): would further reduce r(δ,LD) and potentially validate separate LPs. Sources: BofA/ICE OAS by sector (Bloomberg), or Compustat leverage. Not yet implemented.
