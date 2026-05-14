# Key Empirical Findings
**Last updated:** May 12, 2026

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

## Baseline IRFs (v2 residualized, 2001Q1–2019Q4, σ̂=14.0 pp)

**δ → unemployment:**
- Rises from +0.18 pp (h=0) to plateau +1.50 pp (h=17–20)
- Significant throughout (p<0.001 from h=4)
- 1997Q1+ extension: peak +1.07 pp — attenuation from pre-2001 expansion years

**δ → vacancy rate:**
- Falls from −0.06 pp (h=0) to trough −0.44 pp (h=14)
- Significant p<0.01 from h=3, sustained through h=20
- **Severity placebo passes decisively:** both u^nat level and change interactions insignificant at ALL h=0..20 — cleanest result in paper

**LD → unemployment:** Peak +1.68 pp, significant throughout. Rises monotonically through h=20 (anomaly; inconsistent with ρ_LD≈0.3).

**LD → vacancy rate:** Peak −0.45 pp at h=20. Persistently negative despite residualization (see §LD Vacancy Interpretation below).

**QU → unemployment/vacancy (PLACEBO — FAILS):** Peak ≈ +1.32/−0.44 pp, significant from h=4/h=7. Near-identical to δ and LD. QU is supposed to be a flat-zero placebo.

## Instrument Correlation

| Pair | v1 | v2 |
|------|----|----|
| r(δ, LD) | 0.389 | 0.423 |
| r(δ, QU) | −0.004 | — |
| r(LD, QU) | −0.492 | — |

r(δ,LD) = 0.423 **survives v2 enriched residualization** (controls for aggregate productivity + lagged industry VA + lagged tightness). The shared variation is not industry demand or aggregate productivity cycles. Leading hypothesis: **industry-specific financial conditions** — when credit tightens in a sector, firms simultaneously exit (↑δ) and lay off at survivors (↑LD). Consistent with NFCI state-dependence finding.

## QU Placebo Failure — Smoking Gun for Demand Contamination

Within-instrument corr(β_unemp_h, β_vac_h) for h=0..20:
- QU: **−0.971** — when unemployment rises x pp, vacancies fall proportionally at every horizon
- δ: **−0.948**

This near-perfect mechanical proportionality is the signature of instruments that identify demand-side contractions (which move u and v symmetrically along the Beveridge curve), not structurally distinct shocks with different vacancy implications. QU contamination is structural (industry-level quit dynamics), not fixable with macro controls.

## Joint LP Results (part5_joint_lp.py, v2 residualized)

δ and LD simultaneously in same regression at each h:
- **β^δ vacancy**: stable vs. separate (differences <0.13 pp) — δ result not an LD artifact ✅
- **β^LD vacancy**: collapses to zero and insignificant at all h (replicates colleague's result)
- **β^δ unemployment**: shrinks by ~0.5–1.3 pp vs. separate — δ absorbs shared variation

Interpretation: the zero LD coefficient reflects identification limits (r=0.423 collinearity), not structural confirmation of reposting. The δ vacancy IRF is robust; LD is not separately identified.

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

## Paper Draft Status (as of May 14, 2026)

- **Introduction**: Complete
- **Section 2 Environment**: Complete
- **Section 3 Equilibrium**: Substantially complete (3.10 δ-vs-s mechanism pending)
  - `prop:bgm_nest` (Nesting LR-BGM) + LR-BGM definition paragraph: **Complete**
  - Proof of `prop:bgm_nest` in Appendix C: **Complete**
- **Section 4.1 Instrument construction**: Complete (`fig:delta_distribution` standalone panel added)
- **Section 4.2 LP specification**: Complete (eq:lp generic y_{s,t}; `fig:delta_uv_irf` placed here)
- **Section 4.3 Results**: Partially drafted (fig:delta_uv_irf prose; severity placebo and joint LP paragraphs pending)
- **Section 5 Quantitative**: Not revised
- **Section 6 Conclusion**: Not started
- **Appendix C proofs**: Props 1–4 complete (prop:bgm_nest added May 14)

Key notation: χ_t^c (not x_t^c), Λ_t ≡ F_χ(χ_t^c), eq:/tab:/fig:/sec: label prefixes throughout. Draft file: `Draft/Draft_11May2026.tex`.

## prop:bgm_nest — Proof Structure (Appendix C)

Three-block proof mirroring prop:ags style:
1. **Setup**: p_0=0 → F_χ(0)=0 → X_t^c=0, Λ_t=1, δ_{e,t}=δ_{t-1}
2. **Firm LOM** (eq:N_lom_eq): substituting Λ=1 + gross-output identity Y^c = z(1-u)ρ/μ → parenthesized term = N^e_{t-1} = BGM entry flow
3. **Business formation Euler** (eq:N_euler_eq): Λ_{t+1}=1, X^c=0 → μ-ratio form with ν^f = f_e ρ(N) = BGM equity asset pricing
4. **Resource constraint** (eq:gdp, NOT eq:rc): Y^c = C + X_t under p_0=0; net X_t → Y_t = C_t + ν^f N^e = BGM eq. (3)
5. **Closing remark**: isomorphism holds for any realization of {1-u_t, X_t}, not that equilibrium distributions coincide

Key: isomorphism at GDP level (eq:gdp), not gross output (eq:rc). X_t nets out as intermediate input on both sides.
