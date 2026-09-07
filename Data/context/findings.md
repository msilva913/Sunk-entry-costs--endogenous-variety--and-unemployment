# Key Empirical Findings
**Last updated:** September 6, 2026. The empirical (LP) results below are unchanged since
June 5. Model-side results were re-run after the September 5–6 code fixes; see §D1/M5.

## Shock Persistence (part6, 2001Q1+ window)

| Series | ρ raw | ρ resid | corr(η^δ,η^LD) |
|--------|-------|---------|----------------|
| δ | 0.617 | 0.647 | — |
| LD | 0.489 | 0.305 | +0.488 (resid) |

Note: ρ_δ resid > ρ_δ raw because residualization removes industry-demand variation that is *less* persistent than the structural destruction process. Both bracket a narrow range (0.617–0.647); calibration is insensitive to choice. Model zero-correlation assumption (corr=0) is violated throughout.

VAR(1) panel calibration (part2d, industry×time FEs, 2005Q3–2021Q4):
- ρ_δ = 0.600, ρ_LD = 0.510 (within-estimator; lower than aggregate due to time FE absorption)
- β(LD←δ lag-1) = +0.228 — endogenous exit mechanism
- corr(u^δ, u^{LD}) = +0.229

## Baseline IRFs (BED Deaths, v2 residualized, 2001Q1–2019Q4, σ̂=17.4 pp)

**δ → unemployment:**
- Rises from +0.22 pp (h=0) to plateau +1.71 pp (h=17–20)
- Significant throughout (p<0.001 from h=0)
- 1997Q1+ extension: peak +1.07 pp — attenuation from pre-2001 expansion years (stale; needs rerun)

**δ → vacancy rate:**
- Falls from −0.11 pp (h=0) to trough −0.58 pp (h=18)
- Significant p<0.01 from h=2, sustained through h=20
- **Severity placebo passes decisively:** both u^nat level and change interactions insignificant at ALL h=0..20 — cleanest result in paper

**LD → unemployment:** Peak +1.68 pp, significant throughout. Rises monotonically through h=20 (anomaly; inconsistent with ρ_LD≈0.3).

**LD → vacancy rate:** Peak −0.45 pp at h=20. Persistently negative despite residualization (see §LD Vacancy Interpretation below).

**QU → unemployment/vacancy (PLACEBO — FAILS):** Peak ≈ +1.32/−0.44 pp, significant from h=4/h=7. Near-identical to δ and LD. QU is supposed to be a flat-zero placebo.

## Instrument Correlation

| Pair | v1 | v2 (closings×π) | v2 (Deaths) |
|------|----|----|-----|
| r(δ, LD) | 0.389 | 0.423 | **0.334** |
| r(δ, QU) | −0.004 | — | — |
| r(LD, QU) | −0.492 | — | — |

r(δ,LD) = 0.334 with BED Deaths + v2 residualization (down from 0.423 with closings×π). Switching to Deaths strips the temporary-shutdown component that was correlated with layoffs, reducing the residual correlation. Leading hypothesis for remaining correlation: **industry-specific financial conditions**. Consistent with NFCI state-dependence finding.

## QU Placebo Failure — Smoking Gun for Demand Contamination

Within-instrument corr(β_unemp_h, β_vac_h) for h=0..20:
- QU: **−0.971** — when unemployment rises x pp, vacancies fall proportionally at every horizon
- δ: **−0.948**

This near-perfect mechanical proportionality is the signature of instruments that identify demand-side contractions (which move u and v symmetrically along the Beveridge curve), not structurally distinct shocks with different vacancy implications. QU contamination is structural (industry-level quit dynamics), not fixable with macro controls.

## Joint LP Results (part5_joint_lp.py, BED Deaths + v2 residualized)

δ and LD simultaneously in same regression at each h:
- **β^δ vacancy**: qualitatively stable — same sign, timing, shape; trough −0.535 pp (joint) vs −0.582 pp (separate); confidence bands overlap at most h ✅
- **β^LD vacancy**: collapses to zero and insignificant at all h (0/21 horizons p<0.10)
- **β^δ unemployment**: shrinks by ~0.6–0.8 pp vs. separate; joint peak +1.15 pp, remains significant
- **β^LD unemployment**: significant at only 4/21 horizons (p<0.10) in joint spec

Interpretation: the zero LD coefficient reflects identification limits (r=0.334 collinearity), not structural confirmation of reposting. The δ vacancy IRF is qualitatively robust; LD is not separately identified. Note: the `<0.13 pp` stability claim from closings×π no longer holds exactly — joint-separate differences now reach ~0.48 pp for vacancy — but sign and shape are preserved.

## Recession-Severity Placebo (part7c)

Augmented LP adds B^k × ū^nat_t and B^k × Δu^nat_t interactions:

| Instrument-Outcome | β_h survives | Conclusion |
|---|---|---|
| δ → u | ✅ | Clean |
| δ → v | ✅ | **Cleanest result — BOTH interactions insig. ALL h** |
| LD → u | ✅ | Survives; anti-confound pattern at long horizons |
| LD → v | ⚠️ | Marginal; momentum confound at some horizons |
| QU → u/v | ❌ | Structural contamination, cannot be fixed |

## Granger Causality (part2c, residualized, F-statistics)

| Horizon | δ→LD | LD→δ |
|---------|------|------|
| h=0 | 4.89 | 4.68 |
| h=4 | 8.08 | 3.43 |
| h=8 | 8.89 | 1.82 |
| h=12 | 3.42 | 1.07 |

δ→LD persistent and significant through h=12. LD→δ fades below critical value by h=8 and collapses at h=12. Supports δ as more primitive shock; endogenous exit (δ drives future LD) but not vice versa at long horizons.

## NFCI Interaction (part5)

**δ × NFCI:** β_h insignificant beyond h=0 at average conditions; interaction δ_h highly significant throughout. Entire δ unemployment effect concentrates in tight financial conditions. Puzzle relative to model (no financial frictions).

**s × NFCI:** β_h grows and becomes significant from h=7 onward even at average conditions. Monotonically rising s IRF survives NFCI interaction — not explained by financial amplification.

> ⚠️ **The SLOOS companion results are stale.** `lp_irf_delta_sloos_*.csv` (2026-04-09) and
> `lp_irf_delta_gfc.csv` (2026-05-12) both predate the May 14 switch to BED Deaths + v2;
> their `instr_sd` is 11.16, not 17.4. Re-run before citing either in the draft.

## LD Vacancy Interpretation — Why Persistently Negative

JOLTS LD mixes three events: (a) exogenous match dissolution → reposting [model's s shock], (b) deliberate workforce reduction → no reposting, (c) partial establishment closure → vacancy withdrawn. Cases (b) and (c) are firm-level idiosyncratic, survive industry-VA residualization, and produce simultaneous LD + vacancy withdrawal. **Resolution requires firm-level data.** Decision: LD→vacancy moved to appendix; τ calibrated as external moment (3.1%/month) rather than identified from IRF.

## Motivating Descriptive Evidence

**Cross-recession scatter** (part9): Unemployment CUG scales monotonically with cumulative exit rate across 2001, 2008–09, 2020 recessions. CVS pattern less clean (2001 vacancy formation unusually slow).

**Cross-state scatter** (part10, 50 states):
- r_w(CUG, avg_delta_rate) = **+0.647** — headline result
- r_w(CVS, avg_delta_rate) = **+0.282**
- Recovery speed has wrong/null sign due to industry-composition confound (high-delta = Sun Belt; low-delta = Rust Belt). CUG is the right metric.

**JF Table 2 extension** (part11, ext_2019 preferred):
- Entry share of gross gains: 20.0%; exit share of gross losses: 19.5%
- Hamilton rel. std dev: C3=0.269, C4=0.239 (HP: C5=0.275, C6=0.223)
- C3/C4 fall in extended samples: GFC intensive margin dominates denominator
- Entry/exit margin accounts for ~10% of unemployment-driven flow variation → supports δ/τ ≈ 10% calibration

## SMM Moments Filter Decision (May 19, 2026)

**Baseline filter switched from λ=100,000 to λ=1,600.** Key findings:

- cor(δ,u): −0.041 at λ=100k → **+0.310 at λ=1,600** (sign flip)
- cor(δ,v): +0.127 at λ=100k → **−0.246 at λ=1,600** (sign flip)
- cor(δ,s): +0.492 → +0.681
- No trend contamination: cor(cyc_1600_δ, trend_100k) = 0.005

The λ=100k reversal is confirmed as an over-smoothing artifact: across HP-1600 and the Hamilton (2018) regression filter, cor(δ,u)>0 and cor(δ,v)<0 are robust at h=0..+3. HP-100k is the outlier at 7 of 18 cells in Appendix Table B (app:filter_robustness).

**New Table 4 values (λ=1,600, log-HP, pairwise-complete):**

| Series | σ(x) | cor(u) | cor(v) | cor(s) | cor(f) | cor(δ) | cor(N^e) | cor(z) | ρ₁ |
|--------|-------|--------|--------|--------|--------|--------|----------|--------|-----|
| u | 0.1496 | 1 | −0.804 | 0.577 | −0.791 | **+0.310** | −0.105 | −0.009 | 0.792 |
| v | 0.1419 | −0.804 | 1 | −0.422 | 0.620 | **−0.246** | 0.482 | 0.306 | 0.902 |
| s | 0.1019 | 0.577 | −0.422 | 1 | 0.020 | 0.681 | −0.345 | −0.228 | 0.250 |
| f | 0.0977 | −0.791 | 0.620 | 0.020 | 1 | 0.213 | −0.226 | −0.212 | 0.808 |
| δ | 0.0793 | +0.310 | −0.246 | 0.681 | 0.213 | 1 | −0.496 | −0.250 | 0.114 |
| N^e | 0.0813 | −0.105 | 0.482 | −0.345 | −0.226 | −0.496 | 1 | 0.381 | 0.545 |
| z | 0.0128 | −0.009 | 0.306 | −0.228 | −0.212 | −0.250 | 0.381 | 1 | 0.766 |

σ(θ)/σ(z) = 11.70


## Model Mechanism Findings (Comparisons A–D)

Shared settings unless noted: `b_ratio=0.9, x_v=0.5`; B and D additionally use PATH B
(`dest_elast_target=5.0`). Files and outputs are listed in [`pipeline.md`](pipeline.md)
Part 2. ⚠️ All four were produced at δ_e/τ = 0.210 and would need regenerating under
[D1](decisions.md).

**A — level of δ (CK timing).** High-δ (δ_e = τ, pure exogenous) shows a larger N response
to a z shock (LOM multiplier 4.8×) and a larger Q response to a δ shock (short-duration
asset amplification). The baseline (low δ_e) shows more persistence after a δ shock because
N recovers slowly.

**B — endogenous vs. exogenous exit.** ⚠️ **Redesigned May 22, 2026.** The current design
(`run_solution_endog_exit.jl`, and what the draft §5.3 reports) equates δ̄ *and* τ across the
two specifications so that the residual difference is purely the χ^c margin: the exogenous specification sets
δ_{e,exog} = δ̄_endog ≈ 0.326%/month, about half of δ_{e,endog} ≈ 0.653%, and s rises to hold
τ = 3.1% fixed.

Current finding: the endogenous-exit specification **recovers faster**, the mirror image of
Comparison A. From `eq:N_lom_ll` the entry coefficient is δ_e and the survival coefficient
(1−δ_e), so the higher steady-state destruction rate implies a higher *replacement* rate. The
exogenous specification's low-churn firm stock keeps N depressed longer, sustaining the
ρ↓ → w^int↓ → θ↓ chain and prolonging unemployment. The χ^c margin is a secondary channel
whose sign is ambiguous (a δ↑ impulse raises per-firm profit via market share but lowers ν_f
via ρ).

The **superseded** "exposure-effect" design — same `dest_ann`, same δ_e, different δ̄, giving
"δ-shock response is smaller under endogenous exit because δ̄ is halved" — is retained in the
runner's header as a possible appendix robustness check. Do not cite it as the Comparison B
result.

This also resolves the apparent ψ conflict: **ψ ≈ 0.014 came from PATH A** (the old
`Xc_Y`-targeted design, where ψ is an outcome), while **ψ ≈ 0.033 is the PATH B value**
(`dest_elast_target = 5.0`, `steady_state.jl:607` and draft §5.3). Both are correct for their
own path; only 0.033 applies to the current comparison. See [D2](decisions.md).

**C — variety effects (N → ρ → w^int → JCC).** Implemented by forcing ρ ≡ 1 (replacing f[13]
post-generation) with a dedicated `calibrate_shares_no_variety()`. Isolates the variety
externality operating through recruiter compensation.

**D — entry elasticity ξ_inv (baseline 1.0 vs. near-free-entry 0.1, both recalibrated).**
Lower ξ_inv *amplifies* z-shock amplification (larger entry collapse via the G(Q) channel)
and *attenuates* the δ-shock unemployment response (stronger entry cushion from duration
shortening). **Headline:** at the calibrated δ_e/τ, negative u–v comovement after a δ shock
is robust across the full ξ range including free entry — the empirically calibrated δ̄_e
rules out the CK Beveridge-curve shift. Appendix `app:comparison_D`; the K decomposition
(`eq:K_decomp`, `eq:Q_ll_delta`, `eq:e_ll_delta`) gives the corrected transmission chain
K↑ → Q↑ → e↑ for δ shocks, with sign asymmetry ê^δ > 0 vs. ê^z < 0.

## Where the draft status lives

Section-by-section draft status, proposition inventory, and proof structures moved to
[`draft_status.md`](draft_status.md) on September 5, 2026. This file is empirical results only.

---

## D1/M5 diagnostic — re-run September 6, 2026 (supersedes the retracted first run)

The September 5 run was invalid: it predated the `SS_numeric` fix (residuals 22/14/12), the
`hp_filter` fix (a first-difference smoother inflating σ(labor_prod) ~5x), the LOM convention
restoration, and the X-formula fix. Re-run after all of them, residuals now 3.6e-15 to 1.8e-13.

Θ_e are still hand-set (`b_ratio=0.9, x_v=0.5, ξ_inv=1.0, p_0=0.5, ω_δ=0.5,
dest_elast_target=5.0`), and ρ_s, σ_s remain placeholders. Read the ordering, not the levels.

| | BED (0.0320) | CODE (0.0754) | BGM (0.0963) | data |
|---|---|---|---|---|
| δ_e/τ | 0.087 | 0.210 | 0.271 | — |
| σ(θ)/σ(labor_prod) | 1.02 | 2.21 | 2.93 | **11.70** |
| δ→u peak (pp) | 0.012 | 0.027 | 0.035 | +1.71 |
| δ→u peak horizon (qtr) | 1 | 1 | 1 | **17–20** |
| δ→v trough horizon (qtr) | 4 | 3 | 3 | **18** |

**What changed from the retracted run.** Amplification roughly tripled (was 0.48/0.85/1.11),
almost entirely from the `hp_filter` fix removing spurious variance from the denominator.
σ(labor_prod) is now 0.0143 against 0.0128 in the data, so the "3x too volatile" symptom
behind **M7** is resolved. The degeneracy is not: cor(labor_prod, N^e) = 0.99, so model
measured productivity is still nearly a monotone function of entry.

**What did not change.** The δ→u peak sits at h=1 quarter for every δ_e in the range. The
hump-shape gap (**M8**) is invariant to δ_e and cannot be fixed by choosing it.

**What δ_e does buy.** Amplification is monotone and materially sensitive: 1.02 → 2.93 across
the range, a factor of ~2.9. So δ_e is not irrelevant to Block M, contrary to the retracted
run's conclusion. It is simply unable to close a gap to 11.70 on its own, falling 4x short
even at the BGM value.

## The placeholder s shock flips the Beveridge curve

Running the BED calibration with σ_s = 0 isolates it:

| | with placeholder s shock | s shock silenced | data |
|---|---|---|---|
| cor(u,v) | **+0.761** | **−0.911** | **−0.804** |
| σ(θ)/σ(labor_prod) | 1.02 | 0.75 | 11.70 |
| σ(u) | 0.0111 | 0.0044 | 0.1496 |

The model's unconditional Beveridge curve is fine, and slightly steeper than the data, which
is the right side to be on before a properly calibrated s shock pulls it back. The positive
correlation is an artifact of ρ_s = 0.90, σ_s = 0.010, which have no empirical basis
(`principles.md` §21). The s shock is also doing most of the volatility work: silencing it
cuts σ(u) by 60%.

**Implication for sequencing.** Block M moments cannot discriminate among δ_e values while the
s process is a placeholder, because it contaminates every correlation in the table and
supplies much of the volatility. Calibrating ρ_s and σ_s is now the binding task for Block M,
ahead of anything involving δ_e.

