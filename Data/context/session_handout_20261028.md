# Session Handout — October 28, 2026 (cloud session)
**Branch:** `costly_vacancy_reposting` (commits `ebe921a`, `c70d760`, `69a70f2`, and this
context pass). **Previous handout:** [`session_handout_20261020.md`](session_handout_20261020.md).

> Written for a session with no conversation memory. This file is the entry point.

---

## 1. What was done

1. **Q2 shock processes rebuilt on the S14 windows**, all in
   `Bartek analysis/part6b_shock_persistence_cyclical.py` (no new files). It writes
   `shock_calibration.csv` and a self-documenting report, `shock_calibration.txt`.
2. **Two bugs in τ and jf fixed** in `observables.py` (`shimer_flows`):
   - the 1994 CPS redesign correction ran in the wrong direction;
   - jf and τ were misdated by one month.
3. **D14 settled** for the baseline: three independent AR(1)s, with s not purged.

## 2. The methodology (part6b; full text at the top of `shock_calibration.txt`)

**Model side.** Three independent AR(1)s with uncorrelated innovations. "Independent" describes
the model's processes, not the estimator.

**Windows (S13, S14).**
- Every series is cut at 2019Q4 *before* HP filtering (λ = 1,600, on logs).
- z, τ and s run from 1951Q1; δ runs from 1992Q3, when BED Deaths begin.
- Each series is filtered over its own window.

**Step 1, z.** Univariate AR(1) on 1951Q1–2019Q4. z is ordered first and its cross-lags are
insignificant, so this is the z equation of any VAR, run on the longer window.

**Step 2, δ.** A (z, δ) VAR(1) on 1992Q3–2019Q4, with z ordered first in the Cholesky step.
Keep only ρ_δ (the own lag) and σ_δ (the innovation SD after removing the part the z innovation
explains).
- The model generates z → χ^c → δ_e itself, so leaving that correlation in σ_δ would count the
  channel twice.
- The cross-lag terms are reported but not used: they are equilibrium outcomes.
- The VAR's own z equation is not used for ρ_z or σ_z.

**Step 3, s.** Univariate AR(1) on the spliced 1951–2019 series. σ_s is **not** purged.
- The model has no endogenous s, so z never moves s in the model. Purging would discard
  variation the model cannot regenerate (D14).
- ρ_s comes from the univariate AR(1), not from the diagonal of a VAR with its cross-lags
  dropped.

**Why s is not τ.** τ = δ_e + (1 − δ_e)s, so s = (τ − δ_e)/(1 − δ_e).

**The splice (s before 1992Q3).**
- Before 1992Q3: s = (τ − δ̄)/(1 − δ̄), with δ̄ = 0.00278, the mean of the first 10 BED quarters
  (1992Q3–1994Q4). From 1992Q3: the measured identity.
- In log-linear terms, the pre-1992 s cycle is k_pre × the τ cycle (k_pre = 1.084).
- The level of δ̄ barely matters for the cycle, because a constant shift and a slow drift both
  go into the HP trend.
- The one thing lost is δ's own cycle before 1992, which enters s with weight w = 0.093. No
  τ-based method recovers it.
- The early anchor cuts the step in ln s at 1992Q3 from −1.93 % (BED-window mean) to −0.36 %.
- A rejected alternative is to regress s on τ after 1992 and impute. In cycles that is k × the
  τ cycle, so it only rescales σ. In levels it extrapolates outside the post-1992 τ range.

**Checks in the report.**
- (a) the log-linear weights;
- (b) the step at the splice;
- (c) the pre-1992 construction tested on the BED window, where the true s is known;
- (d) sub-windows, including 1985–2019 (S14 iv);
- (e) the D14 independence rule;
- (f) Ljung-Box tests.

**τ and jf construction** (`observables.py`, `shimer_flows`, used by both the full build and the
flows-only refresh):
- **CPS 1994:** UEMPLT5 × 1.1 from 1994-01 (Shimer, working paper App. A, fn. 27). The old code
  multiplied *pre*-1994 values by 1.16, which doubled the break: ln τ fell 33.7 % at 1994Q1,
  against 2.8 % now.
- **Timing:** jf_t = 1 − (u_{t+1} − u^s_{t+1})/u_t (Shimer eq. 4), and τ uses u^s_{t+1}. The old
  `u_level[1:]` kept its date labels, so pandas aligned it with u_t and jf_t collapsed to
  u^s_t/u_t.
- **Rebuild option (b) (MS):** `flows_only = True` recomputes only `s` (τ) and `jf` in
  `raw_data.pkl`, from UNEMPLOY, CE16OV and UEMPLT5 via the public FRED CSV endpoint (no key).
  Every other column keeps its earlier vintage. Set it to `False` for a full rebuild.

## 3. Results (monthly; `shock_calibration.txt`)

| | ρ | σ | `SHOCKS_BASE` (code, still) |
|---|---|---|---|
| z | 0.9121 | 0.0128 | 0.902 / 0.0092 |
| δ | 0.5824 | 0.0663 | 0.592 / 0.0669 |
| s | 0.8054 | 0.0711 | 0.8741 / 0.0854 |

- **σ_z rises 39 %**, from the pre-1985 volatility (Great Moderation). See the sub-windows.
- **D14 rule passes narrowly.** z explains 9.2 % of the s innovation (threshold 10 %), corr −0.30,
  t = −4.38: 1951–84 13.3 %, 1985–2019 2.4 %. It read 19.1 % before the timing fix. Report the
  correlation; a z → s loading is a robustness check.
- **Ljung-Box** rejects the AR(1) for z and δ. For s, p = 0.29 at lag 4 and 0.014 at lag 8.
- **τ means:** 1951–93 0.0356 (Shimer: about 3.5 %); BED window 0.0258. δ_e/τ on the BED window
  is 0.0873.

## 4. Next, in order

1. **With MS's approval:** put the S14 values in `SHOCKS_BASE` (`run_solution_core.jl`) and `SHK`
   (`sigma_determinacy.jl`), then rerun the R8 and `Q2_*.jl` programs. Julia is not installed in
   the cloud container.
2. **Recompute the τ̄ target** (3.10 %/month) and every Block M moment that uses τ or jf
   (`observables_moments.py`) on the S14 windows. Then Q3.
3. **Align `part8_shimer_national.py`:** it scales post-1994 values by 1.16, the right direction
   but a different factor.
4. **Security:** `observables.py:24` hardcodes the FRED API key in a tracked file. Move it to
   `.env` and rotate the key.
5. Before citing Shimer in the draft, check the page and footnote numbers in the published *RED*
   version (doi:10.1016/j.red.2012.02.001).
6. **`part6c_shock_var_trivariate.py`** imports part6b, so its committed outputs are pre-S14.
   Rerun or retire it (MS).

## 5. Environment notes (cloud)

- No conda and no Julia in the container. Python 3.13 plus a venv with the `environment.yml`
  pins (numpy 2.4, pandas 3.0, scipy 1.18, statsmodels 0.14, matplotlib 3.11, requests,
  fredapi).
- With them, the old part6b reproduced the committed `shock_calibration.csv` to about 1e-12.
- Pushing to `costly_vacancy_reposting` works.
