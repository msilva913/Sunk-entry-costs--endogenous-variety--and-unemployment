# Pending Tasks
**Last updated:** May 21, 2026

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

### Comparison C — Variety Effects (N → ρ → w_int → JCC)   ❌ PENDING
- Files to write: `run_solution_no_variety.jl`, `plot_variety_comparison.jl` (or extend
  `plot_mechanism_comparison.jl`)
- Template: `run_solution_general_CES.jl` (has ζ parameterization)
- Implementation: set ζ=0 so ρ=N^0=1 regardless of N; recalibrate to same targets
- Serialize `irf_no_variety.jls`
- Mechanism: after δ shock, N↓ → ρ↓ → w_int=ρ·z/μ↓ → JCC tightens → vacancies
  suppressed; shutting ζ=0 removes this amplification, directly answering "what can
  your model do that AGS cannot?"
- Color convention: Comparison C alternative uses a reserved color (TBD)

## Empirical Code Tasks (priority order)

1. **Wild cluster bootstrap** — `part5_wcrb.py` ❌
   - n=50 state clusters is at the lower bound of asymptotic SE reliability
   - Use `wildboottest` or `linearmodels` bootstrap; supplement all main LP tables

2. **GFC diagnostic for LD unemployment** — `part7d_ld_gfc.py` ❌
   - LD unemployment IRF rises monotonically through h=20, inconsistent with ρ_LD≈0.3
   - Apply GFC_{t+h} outcome-quarter control dummy to test whether GFC drives the pattern

3. **SLOOS C&I interaction** — `part7b_sloos.py` ❌
   - Replace NFCI_risk with Senior Loan Officer Opinion Survey C&I net tightening (FRED)
   - Tests whether NFCI result reflects credit supply vs. recession severity
   - Directly relevant to industry-credit-conditions hypothesis for r(δ,LD)=0.423

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

1. **Section 3.10 — δ vs. s asymmetry mechanism** ❌
   - Three-channel transmission: u_t (separation), v_{t-1} (predetermination), N_t (variety)
   - Formal proposition: local vs. global scope TBD pending quantitative analysis
   - Bridge sentence to Section 4 LP

2. **Section 4.3 — Results** (partially drafted) ⚠️
   - δ→u IRF numbers confirmed: peak +1.50 pp at h=20, significant from h=0
   - δ→v IRF numbers confirmed: trough −0.44 pp at h=14, significant from h=3
   - Needs: severity placebo paragraph, δ-LD correlation paragraph, joint LP robustness

3. **Section 5 — Quantitative analysis** ❌
   - Update calibration section to four-stage sequential approach
   - BSMM: unconditional moments (Shimer 2005/2012 + BED exit rate)
   - Bayesian IRF matching: δ→u and δ→v from Bartik LP
   - Two-model counterfactual: baseline vs. AGS

4. **Section 6 — Conclusion** ❌

## Completed Code Tasks

- ✅ `run_solution_core.jl` full audit (May 20): naming r/ρ, κ per match, u LOM total v_t, BFE at f[4], SS_symbolics x_c/C/Y corrected
- ✅ `part6b_shock_persistence_cyclical.py`: bivariate VAR(1) HP-1600 log(z)/log(δ), Cholesky z-first → ρ_z^m=0.902, ρ_δ^m=0.592
- ✅ `run_solution_all_delta.jl` (May 20–21): baseline vs. high-δ (δ_e=τ, pure exogenous); b_ratio=0.9, x_v=0.5; serializes `irf_all_delta.jls`
- ✅ `plot_mechanism_comparison.jl` (May 21): 4×2 portrait (u/v, θ/w_int, N/N_e, K/Q); K+Q row reveals short-duration asset amplification; interpreted comments added
- ✅ `run_solution_endog_exit.jl` (May 21): full baseline (endog exit) vs. exog exit; b_ratio=0.9, x_v=0.5; appends exit_flow=δ_e+N; serializes `irf_endog_exit.jls`
- ✅ `plot_endog_exit_comparison.jl` (May 21): 4×2 portrait with exit_flow replacing Q; comprehensive IRF interpretation comments added (δbar halving dominant, ψ≈0.014 dormant x_c)
- ✅ `Inspecting_mechanism_setup.md` (May 21): Comparison B results and two-factor interpretation documented; plot design table and implementation order updated

## Completed Writing Tasks

- ✅ `prop:bgm_nest` (Nesting LR-BGM) + LR-BGM definition paragraph — Section 3 after Remark 1
- ✅ Proof of `prop:bgm_nest` — Appendix C (three-block: LOM eq:N_lom_eq, Euler eq:N_euler_eq, GDP eq:gdp)
- ✅ Section 4.2 LP specification: generic y_{s,t} formulation with label eq:lp
- ✅ Figure placement: fig:delta_distribution in §4.1; fig:delta_uv_irf in §4.2
- ✅ Vacancy convention corrected throughout: pp levels difference, not log (Principle 15)
- ✅ lp_irf_delta_uv_baseline.png now generated by part5_lp.py section 8g (_plt_diagnostics)
- ✅ instrument_delta_distribution.png generated by part3_resid_instruments.py standalone block
- ✅ Section 5.1 CK comparison compressed
- ✅ Section 5.2: filter switched to λ=1,600; Table 4 updated; moments paragraph rewritten; QL justification added; τ→s notation fixed
- ✅ Prop 5 (DS-CES equilibria) restated for DS-CES; detailed 5-step proof in new appendix section
- ✅ Appendix B new subsection `app:filter_robustness`: 3-filter dynamic correlation table with bold HP-100k outlier cells
- ✅ `observables.py`: BED Deaths via BLS API (dataclass 08), 1994 CPS correction, BAWBATOTALSAUS, plt.show() fix
- ✅ `labor_market_dyn_corr_by_filter.py`: standalone script producing dynamic correlation plots and CSV by filter

## Specification Decisions Still Open

- **Primary LP for paper**: joint (part5_joint_lp.py) or separate (part5_lp.py)?
  Current recommendation: separate LP as primary; joint as robustness for δ. δ vacancy IRF stable across both.
- **v3 residualization** (+ monetary policy sensitivity φ_j × ΔMP_t): coded but not yet run; appears in draft as appendix robustness.
- **Industry credit control** (Route B): would further reduce r(δ,LD) and potentially validate separate LPs. Sources: BofA/ICE OAS by sector (Bloomberg), or Compustat leverage. Not yet implemented.
