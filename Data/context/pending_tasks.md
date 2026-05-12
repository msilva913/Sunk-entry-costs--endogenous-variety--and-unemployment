# Pending Tasks
**Last updated:** May 12, 2026

## Code Tasks (priority order)

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

## Specification Decisions Still Open

- **Primary LP for paper**: joint (part5_joint_lp.py) or separate (part5_lp.py)?
  Current recommendation: separate LP as primary; joint as robustness for δ. δ vacancy IRF stable across both.
- **v3 residualization** (+ monetary policy sensitivity φ_j × ΔMP_t): coded but not yet run; appears in draft as appendix robustness.
- **Industry credit control** (Route B): would further reduce r(δ,LD) and potentially validate separate LPs. Sources: BofA/ICE OAS by sector (Bloomberg), or Compustat leverage. Not yet implemented.
