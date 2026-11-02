# The data SVAR (Block B′ data side): setup, findings, implications for the model
**Written:** November 2, 2026 · **Branch:** `costly_vacancy_reposting`
**Program:** `Data/Bartek analysis/q3_data_svar.py`. All data numbers below come from its outputs:
`q3_data_svar.txt` (sections 1–6), `q3_data_svar_headline.txt/.csv`, `q3_data_svar_fevd.csv` and
`q3_data_svar_irf.csv/.png`.
**Model numbers:** `Data/Programs baseline/Q2_amplification_sweep.txt` (S14 shocks, Oct 29).
**Decisions behind it:** S13, S14, S15 and D14 in [`../context/decisions.md`](../context/decisions.md).
Shorter entries are in [`../context/findings.md`](../context/findings.md), "Q3: data SVAR".

> This note is a synthesis. If a number here and in a program output disagree, the program
> wins: rerun it and correct the note. Items marked *back-of-envelope* are arithmetic on program
> outputs, not program outputs, and must not enter the draft (N15). Citations marked
> *unverified* must go through LibKey before they enter the draft (N11–N14).

---

## 1. Setup and specifications

### 1.1 Purpose

D13 (Oct 14) replaced the Bartik local-projection block of the estimation with **Block B′**:
aggregate SVAR responses of u and v to the model's shocks, matched by indirect inference. This
program is the data side. The model side, the same SVAR run on simulated data, is not built yet.

### 1.2 Series and data availability

All series are quarterly, from `Data/raw_data.pkl` (built by `Data/observables.py`).

| Series | Source | Raw availability | Used |
|---|---|---|---|
| z, labor productivity | BLS nonfarm business output per worker, FRED `PRS85006163` | 1951Q1 on | 1951Q1–2019Q4 |
| δ, permanent exit rate | BLS Business Employment Dynamics, employment at dying establishments (`BDS0000000000000000110008LQ5`) / `PAYEMS`, divided by 3 to make it monthly | 1992Q3–2021Q4 | 1992Q3–2019Q4 |
| τ, total separation rate | Shimer method from CPS `UNEMPLOY`, `CE16OV`, `UEMPLT5`. Short-term unemployment × 1.1 from January 1994 (CPS redesign) | 1951Q1 on | input to s |
| s, match separation rate | s = (τ − δ)/(1 − δ) from 1992Q3. Before that, δ is held at δ̄ = 0.00278 (mean of the first 10 BED quarters) | built 1951Q1 on | 1951Q1–2019Q4 |
| u, unemployment rate | FRED `UNRATE` | 1951Q1 on | 1951Q1–2019Q4 |
| v, vacancy rate | Composite Help-Wanted Index to 2000-11, then JOLTS `JTSJOL`; divided by `CLF16OV` | 1951Q1 on | 1951Q1–2019Q4 |

- The composite index is Barnichon's (citation unverified; `observables.py` cites it as *AEJ Macro* 2010, which needs checking).
- **The binding constraint is δ.** BED begins in 1992Q3, so any SVAR that contains δ runs on 1992Q3–2019Q4 only: 110 quarters, all of it after the Great Moderation.

**Two data fixes this SVAR depends on (Oct 28).** Before them, τ had a spurious break of −34 % in log at 1994Q1, which fell inside A's window.
- The 1994 CPS correction now scales short-term unemployment *after* 1994, by Shimer's 1.1. The old code scaled the *pre*-1994 values by 1.16, which doubled the break.
- The job-finding and separation rates are now dated correctly (Shimer eq. 4). Before, pandas aligned u_{t+1} with u_t, misdating them by a month.

### 1.3 The series in levels

Figure: `Data/Bartek analysis/data/results/q3_data_svar_levels.png`, from `q3_data_svar.py`. Each series is shown in levels over 1951Q1–2019Q4, with its HP-1600 trend dashed.

Decade means (back-of-envelope from `raw_data.pkl`, percent):

| Decade | τ (monthly) | δ (monthly) | u | v | Productivity growth (per year) |
|---|---|---|---|---|---|
| 1950s | 3.07 | — | 4.43 | 2.98 | 1.95 |
| 1960s | 3.20 | — | 4.78 | 3.63 | 2.26 |
| 1970s | 3.89 | — | 6.22 | 3.99 | 1.37 |
| 1980s | 4.06 | — | 7.27 | 3.77 | 1.41 |
| 1990s | 3.15 | 0.28 | 5.76 | 3.21 | 2.16 |
| 2000s | 2.70 | 0.24 | 5.54 | 2.68 | 2.30 |
| 2010s | 2.14 | 0.17 | 6.22 | 3.21 | 1.19 |

- **Productivity trends up, with growth breaks:** the 1970s slowdown, the speed-up from the mid-1990s, and the slowdown after 2005. A linear trend misses these.
- **Separations rise and then fall.** τ and s climb into the early 1980s, then fall by about half by the 2010s. u follows a similar arc. The usual explanation is demographic (the baby boom raising turnover, then aging out; citations unverified), together with the secular decline in turnover.
- **δ falls by about 40 % over 1992–2019,** mostly as a level shift during 2001–2010. This is the decline in business dynamism.
- **u and v have no trend but long swings.** The HP-1600 trend follows part of the cycle: it rises through 2009–13, so some of the slow recovery after the Great Recession counts as trend.

**Implications:**
1. **Filtering is needed.** The model is stationary around a fixed steady state, and its shocks are business-cycle shocks. A levels VAR with linear trends would leave the τ and u humps in, and would attribute them to persistent s shocks.
2. **But HP-1600 may remove some persistence that matters,** such as slow recoveries, which is the δ mechanism. The choice of filter is under review (§3.3).

### 1.4 Transformation

- Logs, then an HP filter with λ = 1,600, the same as the moment targets and the model simulations.
- Every series is cut at 2019Q4 *before* filtering (S13), so the pandemic doesn't leak into the 2017–19 cycles.
- Each series is filtered over its own window.
- The SD of each cycle (×100) is: z 1.23, δ 5.96, s 5.77, u 12.48, v 13.82 (report §1).

### 1.5 The three SVARs

| SVAR | Variables | Window | T | Role |
|---|---|---|---|---|
| **A** | z, s, δ, u, v | 1992Q3–2019Q4 | 110 | The only one with δ. Gives δ's share and the full Block B′ target |
| **B** | z, s, u, v | 1951Q1–2019Q4 | 276 | z and s on the long S14 window |
| **C** | z, s, u, v | 1992Q3–2019Q4 | 110 | A without δ, on A's window. A against C measures what omitting δ does |

δ can't enter B, because before 1992Q3 a constant δ̄ has an HP cycle of exactly zero, which would distort every coefficient on δ. So in B, δ's contribution is spread over z, s and the residual. A against C measures that shift on 1992–2019. Assuming the same shift before 1992 can't be tested.

### 1.6 Estimation and identification

- **VAR:** reduced form with a constant, estimated equation by equation by OLS.
  - BIC chooses one lag in all three SVARs, over 1–4 on a common sample. Two lags are always reported too.
  - All are stable; the largest root is 0.910 (A), 0.884 (B) and 0.906 (C).
- **Identification:** recursive (Cholesky), with the shock block before u and v, and **z always first** (S15, Oct 29).
  - **z first.** As in much of the search literature, z is a broad aggregate productivity shock, not technology alone (citation unverified).
    - Ordering it first credits z with all same-quarter comovement between measured productivity and the other variables, including demand-driven movements in productivity; a richer model would separate those. The z share is therefore an upper bound on technology's share.
    - The ordering is conservative for δ. With z first and δ last (the preferred ordering), δ's share is at its minimum: about 12 % of var(u) and var(v).
    - The model generates z → δ_e within the period, so any exit that z drives counts toward z; the δ share is the δ shock's own contribution only.
    - The model-side SVAR must use the model's *measured* labor productivity, not the exogenous z state, with the same ordering.
  - **Preferred ordering for A: z, s, δ**, because it is model-consistent: s is purely exogenous in the model, while measured exit responds to the cycle within the period. **z, δ, s** is reported for completeness. The two differ only in how the same-quarter δ–s comovement is split; z, δ + s and the residual are identical.
  - B and C have a single ordering.
- **Residual.** The innovations to u and v left unexplained by the three shocks are reported as one "Residual". The order of u and v between themselves is arbitrary, so only their sum is identified.
- **Variance shares:** forecast-error variance decomposition at h = 8 quarters, the business-cycle horizon, and at h = 200, the unconditional share, which is the model-side object.
- **Impulse responses:** one-SD shocks over h = 0–12, from A in the preferred ordering. Bands come from a residual bootstrap: recursive design, initial values fixed, 1,000 draws, 68 % and 90 % percentiles. Bootstrap seeded inside the function (`SEED = 20261029`), so bands reproduce across runs and machines.

---

## 2. Key findings

### 2.1 Variance decomposition (percent; one lag)

**A, preferred ordering (z, s, δ)**

| | z | s | δ | Residual |
|---|---|---|---|---|
| var(u), 8q | 37.3 | 14.4 | 12.4 | 35.9 |
| var(u), long run | 35.2 | 12.9 | 11.7 | 40.2 |
| var(v), 8q | 33.6 | 9.5 | 12.4 | 44.5 |
| var(v), long run | 32.3 | 9.4 | 11.8 | 46.5 |

**A, alternative ordering (z, δ, s):** the long-run s / δ split is 9.5 / 15.1 for u and 6.5 / 14.6 for v. z and the residual are unchanged.

**B and C (z, s)**

| | B: z | B: s | B: Residual | C: z | C: s | C: Residual |
|---|---|---|---|---|---|---|
| var(u), long run | 32.8 | 25.0 | 42.2 | 42.1 | 19.1 | 38.8 |
| var(v), long run | 32.4 | 15.1 | 52.5 | 38.9 | 15.1 | 46.0 |

**With two lags (A, preferred, long run):** var(u) is z 26.4, s 11.3, δ 15.6, Residual 46.7 (`q3_data_svar_fevd.csv`).

**Findings:**
1. **z is the largest single source of u and v fluctuations**, about a third of each, *given* that z is ordered first. With z last (orderings no longer reported), δ would lead. So "z dominates" rests on the z-first assumption, which we adopt on prior grounds.
2. **δ and s are of similar size for u, and δ is larger for v.**
   - Preferred ordering: δ 11.7 against s 12.9 for u; δ 11.8 against s 9.4 for v.
   - The ordering moves only the split. δ + s is 24.6 % of var(u) either way.
   - δ's larger role for v fits the paper's mechanism: product-line exit destroys vacancies.
3. **The three shocks explain at most 60 % of var(u) and 53 % of var(v).** That's an upper bound: ordering them before u and v credits them with any same-quarter effect of demand or financial shocks on measured productivity and separations. The residual, 40–47 %, reflects forces outside the model, measurement error, and propagation a one-lag VAR misses.
4. **Shares barely change between 8 quarters and the long run** (at most 4.3 points). That is mostly mechanical: HP cycles are short-lived. So the "long run" here is a business-cycle share, not a long run in levels.
5. **Leaving δ out inflates s.** On the same window, s's share of var(u) is 12.9 % in A and 19.1 % in C. So B's 25.0 % for s is an upper bound on s's role over 1951–2019.
6. **z is not weaker after 1992.** Its share of var(u) is 32.8 % in B (1951–2019) and 42.1 % in C (1992–2019). Part of S14's rationale was that a post-1985 window makes z look weak; that came from z's correlation with u, not from its variance contribution. The rationale should be restated before the paper relies on it.

### 2.2 Impulse responses (A, one lag, preferred ordering, one-SD shocks, log points × 100)

68 % bootstrap bands in brackets.

| Shock → response | h = 0 | h = 2 | h = 4 | h = 8 | Peak |
|---|---|---|---|---|---|
| z → u | −0.99 [−1.17, −0.74] | −2.19 [−2.50, −1.60] | −2.42 [−2.83, −1.57] | −1.36 [−1.78, −0.49] | −2.42 at h = 4 |
| z → v | +0.78 [+0.46, +1.04] | +2.55 [+1.85, +2.93] | +2.50 [+1.55, +2.98] | +0.97 [+0.16, +1.40] | +2.64 at h = 3 |
| s → u | +1.12 [+0.87, +1.31] | +1.54 [+1.06, +1.80] | +1.39 [+0.79, +1.68] | +0.48 [−0.03, +0.77] | +1.54 at h = 2 |
| s → v | −0.36 [−0.59, −0.07] | −1.57 [−1.91, −1.01] | −1.28 [−1.65, −0.62] | −0.20 [−0.50, +0.26] | −1.57 at h = 2 |
| δ → u | +0.35 [+0.12, +0.56] | +1.28 [+0.84, +1.54] | +1.41 [+0.83, +1.73] | +0.78 [+0.21, +1.07] | +1.41 at h = 4 |
| δ → v | −0.22 [−0.52, +0.07] | −1.59 [−1.91, −1.02] | −1.48 [−1.84, −0.80] | −0.56 [−0.86, −0.02] | −1.60 at h = 3 |

**Each shock's size in its own units, and the peak responses per unit** (report §5, "Scale"):

| Shock | Impact on own series | Peak u per unit | Peak v per unit |
|---|---|---|---|
| z | +0.54 | −4.47 | +4.87 |
| s | +3.82 | +0.40 | −0.41 |
| δ | +5.24 | +0.27 | −0.30 |

**Findings:**
1. **All three shocks produce a negative u–v comovement** (a Beveridge-curve movement), with peaks after 2–4 quarters and responses gone by about h = 12.
2. **s acts fastest.** u peaks at h = 2 after s, against h = 4 after z and δ. s's effect is front-loaded and short-lived. It is significant at the 68 % level only through about h = 6.
3. **v falls after an s shock** (−1.57 at h = 2). In the model, costless reposting implies v *rises* after s, as the separated positions are reposted. The data argue for little reposting, a low Λ̄_r.
4. **v falls after a δ shock** (−1.60 at h = 3), and u rises persistently (significant through h = 8). That is the paper's mechanism: exit destroys vacancies as well as jobs.
5. **δ propagates well beyond its direct job loss** (*back-of-envelope*).
   - The shock raises δ by 5.24 % of its mean (0.00225 per month), which is about 0.033 % of the labor force in extra exits over the quarter.
   - u's peak rise is 1.41 % of its mean (0.058), about 0.082 % of the labor force.
   - That is 2.5 times the direct job loss. δ is not persistent (quarterly ρ ≈ 0.2), so even counting every displaced worker as still unemployed at the peak, direct job loss explains under half of the u response.
6. **s moves u much less than flow accounting suggests** (*back-of-envelope*).
   - In steady state, d ln u ≈ (1 − u)·d ln s ≈ 0.94 for a permanent change in s.
   - The data's peak is 0.40 per unit, and the response is short-lived.
   - Part of the gap is s's low persistence (quarterly ρ ≈ 0.5), but the data's s shock clearly raises unemployment less, and more briefly, than a pure inflow shock would.

### 2.3 Against the model (working point of `Q2_amplification_sweep.txt`: p_0 = 0.8, Λ̄_r = 0.25)

**Rerun Nov 2, 2026 at the matched (filter-consistent) shock processes** now in `SHOCKS_BASE`: z (0.9628, 0.0062), δ (0.5319, 0.0641), s (0.8048, 0.0437). The earlier version of this section used shocks whose monthly σ was inflated 1.5–2.5× by an inverted formula in part6b's `to_monthly` (Sept 21 – Nov 2). Its conclusion, that the model's total volatility was about right while s was over-propagated about 2×, was mostly an artifact of that bug. Details in findings.md.

The data column is √share × 12.48 (*back-of-envelope*), using A (HP), preferred ordering, long run, one lag. Ranges in the text add two lags and the alternative ordering.

| σ(u), percent | Data | Model | Gap |
|---|---|---|---|
| z | 7.40 | 2.09 | 3.5× too small (3.1–3.5×) |
| δ | 4.27 | 0.27 | 15.8× too small (15.8–18.3×) |
| s | 4.48 | 5.81 | 1.3× too large (1.3–1.5×) |
| Three shocks together | 9.65 (59.8 % of var(u)) | 6.23 | 0.65× (too low) |
| Shares of var(u): z / s / δ | 35 / 13 / 12 | 11 / 88 / 0.2 | — |

**Findings:**
- **s is about right.** With correct shock sizes, the model's s contribution is close to the data's (1.3–1.5× too large).
- **z is 3–3.5× too small.** Its contribution barely changed from the old run (2.09 against 2.15), even though σ_z halved, because the matched z is more persistent (ρ_m 0.963 against 0.912).
- **δ is 16–18× too small,** and is the central problem.
- **The model's total volatility from the three shocks is now too low,** about two thirds of the data's.
- **The Beveridge correlation is still far off.** cor(u,v) depends only on the reposting rate Λ̄_r. It reaches −0.664 at Λ̄_r = 0.10 against −0.925 in the data (`Q2_p0_floor.txt`, rerun).

**Caveat:** the model's numbers are true-shock shares, while the data's are Cholesky shares. The like-for-like comparison needs the model-side SVAR (§3.3).

---

## 3. Implications for the model

### 3.1 The problem is relative transmission, not total amplification

After the Nov 2 correction, s is close to its data contribution, while z is 3–3.5× and δ 16–18× too small. Amplifying all shocks equally (rigid wages of a general kind, a smaller fundamental surplus, a lower b/w gap) would push s over its data share. In the sweep (matched shocks), one lever at a time:
- **ϕ is not a lever** (the Hagedorn–Manovskii cancellation; findings.md "Q2").
- **Raising z's share to 27–31 %** via b/w = 0.95, x_v = 0.2 or ξ_inv = 0.25 still leaves z's σ(u) at 3.4–3.9, about half the data's, and δ at 0.2 %.
- **Costless reposting** brings s's σ(u) down to 3.44, slightly below the data range, but collapses z to 0.40.

So the model needs levers that act on z and δ but not on s. Damping s is now a second-order concern.

**The z gap is an upper bound on the technology-transmission gap**, because the data's z is broader than the model's (§1.6). The model's z is pure TFP, while the data's z also absorbs demand-driven movements in measured productivity. So part of the 3.1–3.5× z gap is definitional, and amplifying the model's z until it matches the data's would overshoot. δ (15.8–18.3× too small) and s (1.3–1.5× too large) are the cleaner targets. The δ gap is the central problem.

### 3.2 Candidate levers, by which shock they act on

These are proposals, not results. None has been implemented. Citations are unverified (N11–N14).

**Damp s only.** Less pressing after the Nov 2 correction, since s is now only 1.3–1.5× too large. The data say s's effect on u is small and brief (§2.2, items 2, 3 and 6).
- **(s1) Measure s without temporary layoffs.** This is the cheapest lever, with no model change.
  - The Shimer τ counts every employment-to-unemployment flow, including temporary layoffs with recall.
  - Using the CPS unemployed by reason (job losers on temporary layoff against permanent job losers) to strip out recalls would lower σ_s and leave s as permanent separations, closer to the model's s.
  - Test: rebuild s, rerun part6b and the sweep, and see how far s's model share falls.
- **(s2) Recall or rehiring in the model.** Give workers separated by s a high, separate job-finding rate. That matches both facts: u up briefly, and v down rather than up, since the worker returns to the same job and no vacancy is posted.
  - The reposting margin is the firm-side analog.
  - Literature on recall and unemployment: Fujita and Moscarini (citation unverified).

**Amplify z only.**
- **(z1) Wage rigidity.** Rigid wages amplify mainly through θ, the channel z works through, while s raises u mostly through inflows, which don't depend on the wage. So rigidity raises z's share much more than s's. This is the standard Shimer-puzzle response (Hall; Hall and Milgrom; citations unverified). It is a model change, and a well-trodden one; check that it doesn't blur the departure from Coles and Kelishomi.
- **(z2) Check the business block first (diagnostic 4 in `pending_tasks.md`).** Compare the model's z-elasticity of θ with the closed form, with the business block on and frozen. If the variety and entry block absorbs z, the fix belongs in that block, not in the labor market.

**Amplify δ only.** δ is the largest gap (§2.3), and the data show δ propagating well beyond its direct job loss (§2.2, item 5).
- **(d1) The δ-propagation decomposition, a diagnostic.** Split the model's δ → u response into shock scale, jobs and vacancies destroyed per exit, and the variety and entry channel. It shows which of the levers below is worth building.
- **(d2) Slower product-line recovery.** Stickier entry makes the loss of product lines, and of their vacancies, more persistent.
  - Risk: it also slows entry's response to z, so it trades δ against z. The decomposition will show whether the net effect is positive.
- **(d3) Worse job prospects for displaced workers.** Workers from closing establishments plausibly have lower job-finding than other separators. This is specific to δ and has empirical support in the displaced-worker literature (citations unverified).
  - Implementation: a second unemployment pool, or a temporary fall in matching efficiency after δ shocks.
- **(d5) The persistence of δ: data-based and probably first in line.**
  - The HP-matched δ process has a monthly half-life of about one month (ρ_m 0.53). Under the Hamilton filter, the matched process has a half-life of about six months (ρ_m 0.89; `shock_calibration.txt` §8). An AR(1) can't fit both, which suggests measured δ is transitory BED noise on top of a slow-moving component.
  - Modeling measured δ as an AR(1) plus measurement noise, and feeding the model only the persistent part, should fit under both filters. It would give the model a more persistent δ shock, which should propagate more.
  - It is a data-side fix, not a model change.
- **(d4) The scale of δ.** BED establishment deaths count only whole-establishment exits. The model's δ is product-line destruction, which also happens inside surviving firms. If BED understates it, σ_δ is too small for what δ means in the model. This is a measurement question, not a model change.

### 3.3 Before choosing a lever

1. **Build the model-side SVAR.** Run the same specification on simulated data (data-length samples, the same HP filter, the same ordering), so the composition gap is measured like-for-like.
2. **Resolve the residual problem in the moment system.**
   - Block M currently targets total σ(u), σ(v) and cor(u,v), which include the 40–47 % the SVAR assigns outside the model.
   - The three-shock model also has a stochastic-singularity problem: its own SVAR would find almost no residual.
   - Options:
     - (a) measurement error in u and v in the simulations;
     - (b) restrict the targets to what the three shocks explain: their impulse responses and the composition of their share;
     - (c) a fourth, hiring-side shock.

     (b) with (a) for robustness looks most consistent with the indirect-inference design (D13). It would need a new decision entry.
3. **Revisit the filter and the shock processes.**
   - The SVAR's HP filter is two-sided, which blurs the recursive timing assumption. A one-sided Hamilton filter avoids this, at the cost of 11 quarters per series (A falls to about 99).
   - The shock processes in part6b take the AR(1) of the *filtered* data as the shock process. The model's simulations are then time-averaged and filtered again, which changes persistence and volatility.
   - A scratch simulation suggests the effect matters mainly for z. After averaging and HP filtering, the model's z has ρ_Q of about 0.68 against 0.76 in the data, and keeps about 77 % of its SD. For s and δ the two biases roughly cancel.
   - The fix is to choose each monthly (ρ, σ) so that the averaged and filtered simulated series matches the data's filtered AR(1) coefficient and SD. That would make part of the z gap disappear mechanically.
4. **Then test the cheap levers first:** (s1), then diagnostics (z2) and (d1). Move to structural changes, (s2), (z1), (d2) or (d3), only if those fall short. Structural changes should go on their own branch, since they change `run_solution_core.jl`.

---

## 4. Caveats

- **Short window for δ.** A has 110 quarters, and the bands are wide; the δ → v response at h = 0 includes zero.
- **BED Deaths are noisy at quarterly frequency.** δ's own equation fits poorly (R² ≈ 0.09 in part6b). Noise inflates δ's variance but not its share of var(u).
- **Recursive identification.** The shares of z, s and δ are upper bounds against same-quarter feedback from demand. The opposite bound (u and v ordered first) is not computed.
- **The splice.** Before 1992Q3, s is the τ cycle scaled by k_pre = 1.084, so in B it carries δ's pre-1992 cycle, with weight w ≈ 0.09.
- **A single VAR over 1951–2019 (B)** averages across the Great Moderation.
- **These are data statistics.** The comparison with the model is valid only after the model-side SVAR (§3.3).
