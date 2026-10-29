# Pipeline — Empirical and Model Code
**Last updated:** October 12, 2026 (Oct 11 theory programs added; `run_prop5_*` relabeled to the new `sec:mechanism` results; X^c exact in the core code). Earlier: September 23, 2026 (added the ξ sweeps and E8; corrected the
`run_solution_delta_target.jl` status. Full inventory re-verification September 5, 2026;
previous update May 12, 2026, which was missing 13 scripts and the entire model side)

Two independent pipelines, plus one standalone utility ([§Paper library](#paper-library-utility)):

- **Empirical (Python)** — `Data/Bartek analysis/`, subdirs `data/cache/`,
  `data/instruments/`, `data/results/`. Produces the Bartik instruments and LP IRFs.
- **Model (Julia)** — `Data/Programs baseline/`. Produces the steady state, the perturbation
  solution, and the mechanism IRFs.

They meet only at estimation, which does not exist yet — see
[`estimation_design.md`](estimation_design.md).

---

## Part 1 — Empirical pipeline (Python)

### Core path (run in this order)

| Script | Function | Status | Key output |
|--------|----------|--------|------------|
| `part1_shares.py` | 2006 QCEW state×supersector shares ω_{s,j} | ✅ | `shares_base2006.parquet` |
| `part2_shock_rates.py` | BED deaths → LOO national g^δ | ✅ | `shock_rates_1992Q3_2024Q2.parquet` |
| `part2_shock_rates_s.py` | JOLTS separations → LOO g^LD, g^QU, g^TS | ✅ | `shock_rates_{ld,qu,s}_*.parquet` |
| `part2b_extend_va.py` | Extends BEA sector VA back from 2005Q1 via constrained Chow-Lin | ✅ | `bea_va_quarterly_12ind_extended.parquet` |
| `part2b_residualize_shocks.py` | **v2** residualization on Δlog p_t + Δlog VA_{j,t−1} | ✅ | `shock_rates_{delta,ld,qu}_resid.parquet` |
| `part3_resid_instruments.py` | Bartik aggregation for ν^δ, ν^LD, ν^QU | ✅ | `{delta,ld,qu}_instrument_resid_base2006.csv` |
| `part4_outcomes.py` | LAUS unemployment + JOLTS state vacancies | ✅ | `laus_quarterly.parquet` |
| `part5_lp.py` | All panel LPs, both outcomes; §8g emits the paper figure; **§11 is the E8 instrument-persistence diagnostic** (within-state autocorrelation of B̃^δ + augmented LP with p lagged instruments) | ✅ | IRF CSVs + plots, `irf_delta_augmented_comparison.png` |

### Supporting and diagnostic scripts

| Script | Function | Status |
|--------|----------|--------|
| `construct_delta_instrument.py` | Library: QCEW shares, BED rates, δ Bartik aggregation | ✅ |
| `construct_s_instrument.py` | Library: JOLTS rates, s Bartik aggregation | ✅ |
| `part2b_residualize_shocks_v3.py` | **v3**: adds monetary-policy sensitivity φ_j × ΔMP_t | ⚠️ written, **never run** — task E6 in [`pending_tasks.md`](pending_tasks.md) |
| `part2c_granger_lp.py` | Panel LP Granger causality, δ vs. LD | ✅ |
| `part2d_var_calibration.py` | Bivariate panel VAR(1) on (ν^δ, ν^LD) | ✅ |
| `part3_instrument.py` / `part3_instrument_s.py` | Raw (non-residualized) Bartik aggregation | ✅ appendix only |
| `part5_joint_lp.py` | Joint LP: δ and LD in the same regression | ✅ robustness, not structural ([D6](decisions.md)) |
| `part5_joint_diagnostics.py` | Partial-R²/first-stage test of Story 1 vs. Story 2 for the zero β^LD | ✅ |
| `part6_shock_persistence.py` | AR(1) persistence of the raw and residualized shocks | ✅ |
| `part6b_shock_persistence_cyclical.py` | **Model shock calibration (S14, Oct 28):** (1) z AR(1) on 1951Q1–2019Q4; (2) (z, δ) VAR(1) on 1992Q3–2019Q4, Cholesky z-first, keeping ρ_δ and σ_δ ⊥ z; (3) s AR(1) on the spliced 1951–2019 series (before 1992Q3, s = (τ − δ̄)/(1 − δ̄) with δ̄ the first-10-BED-quarter mean), σ_s not purged. Series cut at 2019Q4 before filtering (S13), each filtered over its own window. Report includes the method, splice checks, sub-windows, the D14 independence test and Ljung-Box. Run code under `if __name__ == "__main__":` | ✅ Oct 28 → `shock_calibration.csv/.txt` (ρ, σ for z, δ, s; τ as diagnostic) |
| `q3_data_svar.py` | **Q3 Block B′ data side (Oct 29):** recursive SVARs A (z, δ, s, u, v) 1992Q3–2019Q4, B (z, s, u, v) 1951–2019, C (z, s, u, v) 1992Q3–2019Q4; cycles from part6b plus u, v; BIC lag (p = 1, 2 also reported); FEVD at h = 4, 8, 20, 200 with min–max over orderings of the shock block; IRFs of u, v at h = 0–12 with residual-bootstrap 68 %/90 % bands (1,000 draws). part6c (trivariate VAR) deleted Oct 29 | ✅ Oct 29 → `q3_data_svar.txt`, `_fevd.csv`, `_irf.csv`, `_irf.png` |
| `part7_nfci.py` | Chicago Fed NFCI → quarterly | ✅ |
| `part7b_sloos.py` | SLOOS C&I tightening in place of NFCI_risk | ⚠️ run, but **output stale** (pre-Deaths instrument, `instr_sd`=11.16) |
| `part7c_recession_placebo.py` | Dual-severity placebo (u^nat level + change) | ✅ headline validation |
| `part8_shimer_national.py` | Diagnostic only: national Shimer flows; does not feed the LPs | ✅ |
| `part9_recession_scatter.py` | Cross-recession motivating scatter | ✅ |
| `part10_state_scatter.py` / `part10b_..._sensitivity.py` | Cross-state motivating scatter (50 states) + x-axis sensitivity | ✅ |
| `part11_jf_table2.py` | Jaimovich-Floetotto (2008) Table 2 replication + extension | ✅ |
| `refresh_bed_cache.py` | Refetch BED from the BLS API | ❌ not run — cache truncates at 2021Q4 |
| `run_v2_locally.py` | Orchestrates part2b → part3_resid → part5 (needs `FRED_API_KEY`) | ✅ |
| `build_report_html.py` | Pandoc Markdown→HTML with base64 figures | ❌ needs local pandoc |
| `BEA_diagnostic_script.py`, `diagnostics.py`, `jolts_coding_test.py`, `calibration_shocks.py` | Ad-hoc diagnostics; not part of any run path | — |

### Not yet written

| Script | Purpose | Blocking |
|---|---|---|
| `part5_wcrb.py` | Wild cluster bootstrap → **full 42×42 Ω_β** (not just SEs) | 🔴 estimation |
| `moments_bootstrap.py` | Block bootstrap → Ω_m, block length 8 | 🔴 estimation |
| `part7d_ld_gfc.py` | GFC dummy diagnostic for the monotone LD unemployment IRF | 🟡 |
| **E9 model-side LP** (name TBD; Julia, or Julia simulation → Python LP) | Simulate a 50-"state" panel and run the baseline LP on it, to settle [D10](decisions.md) | 🔴 estimation — **highest priority** |

> ⚠️ **Known bug.** `overlay_plot_irf` in `part5_lp.py` §8 crashes with
> `OSError: Invalid argument` when writing `lp_irf_delta_ld_qu_overlay.png` — a PIL or
> matplotlib version problem in the `econ2315` conda environment, not a logic error. The
> call is wrapped in `try/except` so the rest of §8 still runs. Low priority.

### Python environment (Oct 20, 2026)

**Spec:** [`../environment.yml`](../environment.yml) (env name `econ2315`, conda-forge, minor
versions pinned). This is the cross-machine contract, the Python counterpart of the Julia
`Project.toml`/`Manifest.toml`. On a new machine: `conda env create -f Data/environment.yml`.
After a spec change: `conda env update -n econ2315 -f Data/environment.yml --prune`.

**Check after any environment change:** run `part6b_shock_persistence_cyclical.py` and
`git diff --stat data/results/shock_calibration.csv`. Only last-digit noise is allowed (1e-12
on Oct 20). Then `git checkout` the regenerated outputs.

**Per machine.**
- **Oct 20 machine (m.silva):** Anaconda 2026.07 at `%LOCALAPPDATA%\anaconda3`. Its **base**
  env has every package and passes the check, so it is used directly; there is no `econ2315`
  env here yet. Not on PATH: call `"$env:LOCALAPPDATA\anaconda3\python.exe"` from PowerShell.
- **Other machine (msilv):** the original `econ2315` env, versions unrecorded. Bring it in line
  with the update command above, then run the check.
- Julia's private Conda Python (`~/.julia/conda/3`) is for PyCall/SymPy only. Never use it
  for the empirical pipeline.

`raw_data.pkl` (built by `observables.py`, tracked in git) is a pandas pickle: keep pandas on the
pinned major version, or rebuild it.

**`observables.py` (Oct 28).** τ (`s`) and jf come from `shimer_flows`: UEMPLT5 × 1.1 from 1994-01
(CPS redesign; Shimer working paper App. A, fn. 27), jf_t = 1 − (u_{t+1} − u^s_{t+1})/u_t, and τ from
u^s_{t+1}. `flows_only = True` (the default since Oct 28) recomputes only `s` and `jf` in the existing
`raw_data.pkl` from the FRED CSV endpoint (no key); `False` rebuilds every series at the current vintage.
⚠ The FRED key is hardcoded at line 24; move it to `.env`.

### Run order for a fresh start

```
part1 → part2 → part2b_extend_va → part2b_residualize_shocks → part3_resid → part4 → part5 → part6
```

With v2 residualization (requires a FRED key):

```powershell
$env:FRED_API_KEY = "<key>"
python run_v2_locally.py   # orchestrates part2b → part3_resid → part5
```

After the first run, `data/cache/bea_va_quarterly_12ind.parquet` is cached and the key is no
longer needed.

### Technical conventions

- **BED Deaths (dataclass 08)** — zero third-month employment in each of the four
  consecutive quarters following the last quarter with positive employment, i.e. permanent
  exits. The δ numerator. Deaths ⊆ Closings; weighted mean ratio 0.90 (2001–2019).
  **The death is dated to the quarter of closure, not to the quarter of confirmation**, so
  the series is correctly timed; the rule induces a publication lag and vintage revisions
  only. See [D10](decisions.md) — this is why the recognition rule cannot explain the
  δ→u peak at h = 17–20.
- **LOO** — national shock rates always exclude the state being instrumented (denominator
  only; state-level deaths by supersector are not public).
- **Base year** — 2006 QCEW shares, `agglvl=54`, private sector.
- **v2 residualization** — `log g^k_{j,t} = α_j + γ Δlog p_t + λ Δlog VA_{j,t−1} + ν^k_{j,t}`.
  VA from FRED, extended to 1992Q1 by Chow-Lin + Denton.
- **State FIPS** — zero-padded 2-digit strings throughout.
- **Aggregation** — vacancies are a stock, so quarterly = *average* of monthly; separations
  are a flow, so quarterly = *sum*.
- **COVID cap** — 2019Q4 for all LP outcomes; instruments may extend further.

---

## Part 2 — Model pipeline (Julia)

Everything is run from `Data/Programs baseline/`.

### Environment (this machine, configured September 5, 2026)

The Julia model code needs SymPy via PyCall. Julia's private Conda.jl environment
(`~/.julia/conda/3/x86_64`, a Miniforge3 install) is **broken** — no `conda-meta/history`,
so `Conda.add` fails with `NoBaseEnvironmentError`, and its libmamba solver fails to load.
PyCall was therefore rebuilt against the user's miniconda instead:

```julia
ENV["PYTHON"] = raw"C:\Users\msilv\miniconda3\python.exe"
using Pkg; Pkg.build("PyCall")
```

with `sympy` installed there via `python -m pip install sympy`. To revert to the Conda.jl
Python: `ENV["PYTHON"] = ""; Pkg.build("PyCall")` (but fix the base env first).

Full package set required by the model code: BenchmarkTools, CSV, DataFrames, Distributions,
GLM, KernelDensity, LaTeXStrings, LeastSquaresOptim, MKL, MappedArrays, NLsolve, Optim,
PGFPlotsX, Parameters, Plots, PrettyPrinting, PrettyTables, PyCall, Roots, ShiftedArrays,
StatsBase, SymPy, TexTables, TypedTables.

### Core

| File | Role |
|---|---|
| `solution_functions.jl` | Perturbation toolkit (Salazar-Perez & Seoane 2023). Not modified. |
| `steady_state.jl` | **Canonical.** `calibrate_shares` (4 stages, PATH A/B — see [D2](decisions.md)), `steady_state`, and the free-entry variants. Holds `const TARGETS`. |
| `run_solution_core.jl` | **The 41-equation model** (since R8, Oct 14–19, 2026): draft timing, partial reposting, states [Ñ, Ẽ, Ṽ^u, Ṽ^s, z, δ, δl, s]; equation map in [`model_equations.md`](model_equations.md) "Code status". Helpers: `repost_par(cal)` → (α, rp), `build_PAR(cal, SHK)`, `shock_col(state, sd)`, `const SHOCKS_BASE`. `SS_numeric(PAR, targets)` calls `steady_state(cal)` and errors unless `PAR == build_PAR(calibrate_shares(targets), …)`. `solution_interface` requires that numeric SS and errors on a failed BK count, citing ε* = 2.40 and σ* = 0.011. `SS_symbolics` is a placeholder (zeros): never a linearization point. The steady-state residual check hard-errors — never downgrade it to a warning. `estimate = []`, `priors = (;)` (no estimation code yet) |
| `run_solution.jl` | Thin baseline runner (`build_PAR`, `SHOCKS_BASE`, numeric SS); serializes `model_output.jls` |
| `time_series_fun.jl` | `simulate_model`, `monthly_to_quarterly`, `hp_filter`, `hamilton_filter`, `moments`. ⚠️ `hp_filter` was **fixed Sept 5, 2026** — it had been a first-difference smoother, not HP. `moments(...)` requires `lags=2` (two autocorrelation labels are hard-coded) |
| `impulse_response_plots.jl` | Shared IRF plotting helpers |

### Mechanism comparisons (all complete; figures in the draft)

| Runner | Plotter | Comparison | Output |
|---|---|---|---|
| `run_solution_all_delta.jl` | `plot_mechanism_comparison.jl` | **A** — level of δ (CK timing) | `irf_all_delta.jls`, `mechanism_A_*.pdf` |
| `run_solution_endog_exit.jl` | `plot_endog_exit_comparison.jl` | **B** — endogenous vs. exogenous exit | `irf_endog_exit.jls`, `mechanism_B_*.pdf` |
| `run_solution_no_variety.jl` | `plot_no_variety_comparison.jl` | **C** — variety effects (ρ ≡ 1) | `irf_no_variety.jls`, `mechanism_C_*.pdf` |
| `run_solution_entry_elasticity.jl` | `plot_xi_inv_comparison.jl` | **D** — entry elasticity ξ_inv | `irf_xi_inv.jls`, `mechanism_D_*.pdf` |

### Draft-bound statistics

| Runner | Purpose | Output |
|---|---|---|
| `mechanism_stats.jl` | **Reproducibility (`principles.md` N15).** Reads the four `irf_*.jls` and emits every scalar quoted in §5.3 and the Comparison D appendix: half-lives, peak ratios, entry costs, the exit-flow rate/stock decomposition, the w^int impact gap. 66 quantities. Run after the four mechanism runners | `mechanism_stats.txt` (readable), `mechanism_stats.tex` (`\newcommand` macros) |

The `.tex` file is ASCII-clean and meant to be `\input` by the manuscript so the prose stops
carrying hard-coded figures. MS decided to keep literals in the prose; `mechanism_stats.txt`
is the authority to diff against.

| `beveridge_free_entry.jl` | **Added Oct 6, 2026.** Calibration, tightness coefficient, IRF and sign-flip helpers for the free-entry GS economy. Imported by `gs_free_entry_beveridge.jl` | `beveridge_free_entry.txt` |
| `gs_free_entry_beveridge.jl` | **Added Oct 6, 2026; extended Oct 7.** The R15 proposition's numbers: log-linear weights (including the log-shock weights δ̄/(r+δ̄), f̄δ̄/τ̄), IRF tables at iid and estimated persistence, path slopes, sign-flip boundaries in rho and delta_e, the `Kfixed` flow-cost comparator that isolates the user-cost channel, and (Oct 7) a table of ω, R, ω+R and ρ(1−λ²)/(1+ρλ). **Every closed form is asserted against a simulated IRF of the same linear system**; Oct 7 adds assertions that the note's effective-rate form of a (`a_note`) and the R form of the slope and both conditions match the original algebra at ρ ∈ {0, 0.3, 0.592, 0.9}. Internal formulas for the δ shock use level shocks; the note's log-shock coefficient is δ̄ × the level one. **s-shock section (Oct 7):** `gs_s(c; ρ)` gives b, ℓ_s, R_s, the slope and the part-3 bound in log units; `irf_s` simulates the linear s-shock system; `ρ_s_cal = 0.874` from `parameters.md`. Asserts b at ρ = 1 against the steady-state finite difference d ln θ/d ln s, ℓ_s against a one-step nonlinear response of ln u, the slope formula against the simulated path at ρ ∈ {0.3, 0.592, 0.874, 0.99}, v̂_t > 0 for all t ≥ 1 whenever the bound holds, and θ̂ ≡ 0 for iid s. Prints the s table for D1, 7.54 %/yr and 10 %/yr (δ_e = τ has no s margin), R_δ for comparison, and the D1 s-shock IRF at ρ = 0.874. Includes `beveridge_free_entry.jl`; Base Julia only (no packages) | `gs_free_entry_beveridge.txt` |
| `predetermined_vacancy_share.jl` | **Added Oct 5, 2026.** ⚠ The T threshold was withdrawn Oct 6 (per unit level, not per 1 %); the predetermined share stands. Predetermined share of the vacancy stock (v̄^pre/v̄ = 1 − ē/v̄) and the threshold entry elasticity T = κ_v/(δ̄_e/τ̄) − 1 at four destruction calibrations including δ_e = τ (CK). **Self-contained** — no solve, no `steady_state.jl` include; mirrors its closed-form algebra and asserts all three closed forms against each other. Runs in seconds | `predetermined_vacancy_share.txt` |
| `gs_finite_xi_entry_offset.jl` | **Added Oct 11, 2026.** Finite-ξ GS economy of `prop:fe` (κ = 0, states u, v_pre): perfect-foresight stacked linear solve, sparse coloring Jacobian, T = 2000, BK check via `eigvals(Bm, Am)` (asserts 2 stable roots). h = 1 entry offset on a δ̄ × ξ grid and the threshold δ̄*(ξ). Source of `rem:finite_xi` (3.8%/2.7%/1.6% per year at ξ = 1/2/10; offset 0.85 at D1, ξ = 1) and of the BK claim in `prop:ags` | `gs_finite_xi_entry_offset.txt` |
| `double_limit_check.jl` | **Added Oct 11, 2026.** Steady-state N under f_e = f̄_e/ε as ε grows converges to z̄L̄(1−δ̄)/(f̄_e(r+δ̄)) (`prop:double_limit` (i)) | `double_limit_check.txt` |
| `equilibria_convexity_check.jl` | **Added Oct 11, 2026.** `prop:equilibria`: Φ(θ), condition (C), crossing counts of N_JC and N_RC under capped Cobb–Douglas, DHRW and SMCD. Loads `MatchingSpec` from the companion repo, which must be checked out as a sibling of this repo (`GitHub/Soft-Min-Cobb-Douglas-Matching-Function`; path hard-coded in `SMCD_LIB`) | `equilibria_convexity_check.txt` |

### Calibration diagnostics (not paper comparisons)

| Runner | Purpose | Output |
|---|---|---|
| `run_solution_delta_target.jl` | **D1/M5 diagnostic** — dest_ann ∈ {0.0320 BED, 0.0754 code, 0.0963 BGM}; reports steady state, Block M moments at HP λ=1,600, and δ→u/δ→v IRF persistence in quarters. Decided that δ_e can be fixed. ✅ **run September 6, 2026** (after the `SS_numeric`/`hp_filter` fixes) and **re-run September 21, 2026** with the calibrated s process. Results in [`findings.md`](findings.md) §D1/M5 | `irf_delta_target.jls` |
| `run_prop5_s_check.jl` | **`lem:accounting` Part 1** (old Prop. 5 Part 1) at the p_0 = 0 benchmark: s-shock responses at h = 0, 1 (Δe_t, ΔM_t, Δe_{t+1}, Δu, Δv) checked against the exact accounting, and conditions (i)–(ii) as shares of the inflow d, over dest_ann ∈ {0.0320, 0.0754}, ρ_s ∈ {0, 0.5, 0.874, 0.95}, ξ_inv ∈ {0.5, 1, 2}. Source of the 2.6% / 0.7% bounds in the draft. ✅ run Oct 5, 2026 | `prop5_s_check.csv` |
| `run_prop5_parts23_check.jl` | **`lem:accounting` δ part and `prop:exit_s`** (old Prop. 5 Parts 3 and 2). δ part (δ shock, p_0 = 0, exact): first-month conditions, monthly paths, quarterly averages; source of the 81% offset and the quarterly-sign claims. Exit part (p_0 = 0.5 baseline): g_N, g_u, F′(χ^c), det J of the 3×3 exit system — **indicative** until R8 fixes the code's exit dating. ✅ run Oct 5, 2026; redo the exit part after R8 | `prop5_parts23_check.csv`, `prop5_part3_quarterly.csv` |

| `run_xi_sweep.jl` | ξ_inv ∈ {0.5, 1, 2, 3, 5, 8} at the BED spec — does entry-cost convexity move the δ→u peak horizon? | ✅ run September 21, 2026. Peak never passes h=2 |
| `run_xi_eps_sweep.jl` | Joint ξ_inv × ε sweep — same question plus the variety channel; also maps the Blanchard-Kahn ceiling | ✅ run September 21, 2026. ε < 3 violates BK |

The two sweeps were written as scratch and committed in `dd0fcdf`. They establish the
**structural impossibility of a hump-shaped δ→u IRF** recorded in
[`findings.md`](findings.md), which is what motivated [D10](decisions.md) and E8. Keep them:
they are the evidence for that claim and will be re-run if the model gains a propagation
mechanism.

`run_solution_delta_target.jl` also doubles as the **prototype m(θ) simulator** for Block M:
it is the first code in the repo to build model moments at λ=1,600 over the full series set
{u, v, θ, τ, f, δ_e, N^e, labor_prod}. Promote it to a standalone module when wiring
the estimation (see `estimation_design.md` build step 4). Its `simulated_moments` function
(line 154) and `simulate_model` (`solution_functions.jl:783`) are also the starting point
for **E9**, the model-side LP test.


Comparisons B and D use **PATH B** (`dest_elast_target = 5.0`) plus `b_ratio = 0.9,
x_v = 0.5`, which are *not* the `TARGETS` defaults. See [D2](decisions.md).

### Q1 (R8) programs and helpers — Oct 14–18, 2026 (branch `r8_timing_reposting`)

Since R8 the state vector is [Ñ, Ẽ, Ṽ^u, Ṽ^s, z, δ, δl, s] (41 equations). Runners must build
PAR with `build_PAR(cal, SHK)`, single-shock columns with `shock_col(state, sd)`, and take shock
processes from `SHOCKS_BASE` (S14 part6b values since Oct 29). Never assign numbers to the
symbolic ρ_z, σ_z, … (`gen_model_equations` reads them). `solution_interface` errors with the
qzflag case (0 no bounded equilibrium, 2 indeterminate, 3 invertibility). Environment:
`Project.toml`/`Manifest.toml` + `setup_environment.jl`; run with `--project=.` after
`Pkg.instantiate()`, or in the global environment if that already has the packages.

| Program | Purpose | Output |
|---|---|---|
| `R8_snapshot_pre.jl` | Pre-R8 IRF snapshot (regression reference) | `R8_snapshot_pre.jls/.txt` |
| `R8_test_timing.jl` | Tests 1–2: draft timing; p_0 = 0 reproduces the snapshot to 1.8e-12 | `.txt/.jls` |
| `R8_test_steady_state.jl` | Test 3: steady state with partial reposting, 14 identities | `.txt` |
| `R8_test_reposting.jl` | Tests 4–5: reposting dynamics; costless = step 2; positions conserved | `.txt/.jls` |
| `R8_step5_defaults.jl` | Step 5: defaults converge; δ_e/τ = 0.0873; calibrated vector | `.txt` |
| `R8_step7_baseline.jl` | Step 7: IRFs and simulated cor(u,v) over Λ̄_r (exploratory) | `.txt` |
| `R8_test6_free_entry.jl` | Test 6: `prop:fe` closed forms under AGS emulation; documents the user-cost offset (R19) | `.txt` |
| `sigma_determinacy.jl` | σ* by bisection under four specs (own `SHK` copy) | `.txt` |
| ~~`Q2_shock_processes_check.jl`~~ | Deleted Oct 29 with part6c (it compared part6b against part6c processes); result kept in findings.md "Q2" | — |
| `Q2_p0_floor.jl` | cor(u,v) over (p_0, Λ̄_r): only Λ̄_r matters; p_0 sets the floor | `.txt` |
| `Q2_amplification_sweep.jl` | One-at-a-time levers (b/w, x_v, X_Y, ξ_inv, exit elasticity, ω_δ): ϕ, fundamental surplus, z/δ/s shares of var(u). Never sweeps η_L or ε | `.txt` |
| `bk_epsilon_threshold.jl` | BK threshold ε* by bisection under four specs (`tab:eps_threshold`) | `.txt` |
| `lrbgm_determinacy.jl` | LR-BGM block determinacy over ε, σ (App. K) | `.txt` |
| `cutoff_profit_share.jl` | Flow-profit share s_R of the survival cutoff (App. K) | `.txt` |

### Comparison / legacy variants

`run_solution_CK.jl` and `steady_state_CK.jl` (Coles-Kelishomi benchmark);
`run_solution_general_CES.jl` (general CES — **loads the old Matlab `posterior_mode.mat`**,
which belongs to a previous generation of the model and must not be reused);
`GS.jl` (Gabrovski-Silva JEDC benchmark); `run_solution_{alt,elastic,high_b,kappa,phi,
risk_neutral,simplified}.jl` (single-parameter sensitivity runs). ⚠ **Unmaintained (Oct 18,
2026):** these call `calibrate_labor_share` (defined only in `steady_state_CK.jl` /
`steady_state_general_CES.jl`, which the core does not include) and build PAR by hand, so
they were broken before R8 and were not converted in Q1 step 6. Rebuild from
`run_solution.jl` (`build_PAR`, `shock_col`, `SHOCKS_BASE`) if one is needed.

### Diagnostics and output

`steady_state_checks.jl` (calibration inspection; free-entry branch bug — task M2),
`steady_state_checks_*.jl`, `accuracy_checks.jl`, `comparative_statics.jl`,
`export_endog_exit_irfs.jl` (IRFs + SS → CSV), `generate_tables.jl`,
`second_moments.jl` (⚠️ **stale and invalid**: λ=100,000, only 8 series, AND it predates the `hp_filter` fix of Sept 5, 2026 — every moment it ever produced is wrong; superseded by `simulated_moments` inside `run_solution_delta_target.jl`), `Hill estimator/` (tail-index estimation for the cost distribution).

### Model code conventions (audited May 20, 2026; updated for R8 Oct 19)

- **Naming** — `r` is the discount rate; `ρ` is the relative price N^(1/(ε−1)). Consistent
  with `steady_state.jl`.
- **κ** — matching cost paid *per match* (Pissarides 2009). The resource constraint uses
  `κ·q·v`, not `κ·v`.
- **u LOM** — uses total `v_t = v_pret + e_t`; entrants participate in matching within
  period t (draft `eq:v_lom`).
- **States** (since R8) — `x = [Nt, Et, Vu, Vs, z, δ, δl, s]`: pre-exit stocks, the three
  shocks, and δl = δ_{t−1} (bookkeeping, for the realized δ_e only); 33 controls; monthly.
- **Equation order** — f[1] χ^c, f[2] δ_e, f[3] JCC, f[4] business-formation Euler, f[5] K,
  f[6] w^int, f[7] Nash wage, f[8] vacancy creation, f[9]–f[21] static conditions,
  f[22]–f[24] same-period stocks, f[25]–f[29] state LOMs, f[30]–f[35] observables,
  f[36]–f[38] shocks, f[39]–f[41] reposting (Λ_r, Q^rep, H).
- **Full equation reference with interpretations** — [`model_equations.md`](model_equations.md):
  Block 1 (baseline) and Block 2 (reposting margin, in code since R8; see [D11](decisions.md)).
  Its "Code status" section maps the equations to the code's f[k].
- **Defaults (`TARGETS`, Q1 step 5)** — D1 `dest_ann = 0.0320` + PATH B
  `dest_elast_target = 5.0`, `b_ratio = 0.9`, `x_v = 0.5`, `Λr_target = nothing` (costless
  reposting). A runner changes a target with `merge(TARGETS, (...))`, never by editing the const.
- **PAR and shocks** — always `build_PAR(calibrate_shares(tg), SHK)`, with `SHK = SHOCKS_BASE`
  or a NamedTuple with the same fields; never a hand-built vector (`SS_numeric` errors on a
  mismatch). Shock columns via `shock_col(:δ, SHK.σδ)`.
- **SS_symbolics** — a placeholder since R8 (zeros). Linearize only at `SS_numeric`.
- **Known toolkit bug** — `eval_SS` iterates `for ip in npar` (once) instead of `1:npar`,
  so it substitutes only the last parameter. Every runner passes the numeric SS, so it is
  latent.

### `steady_state.jl` architecture (as of May 31, 2026)

- `calibrate_shares(targets)` → 18-parameter NamedTuple. **PATH A** (default; targets
  `Xc_Y`) or **PATH B** (triggered by supplying `dest_elast_target`; targets the exit
  elasticity, making `Xc_Y` an outcome). ϕ is a residual from Nash in Stage 5.
- `steady_state(para; init=0.51)` → 2×2 solver in (log θ, log χ^c).
- `steady_state_free_entry(para)` → 1D solver for ξ_inv = 0; outer root-find over θ (JCC),
  inner over δ_e. **Currently converges to the wrong (high-θ) branch** — task M2.
- `calibrate_shares_free_entry(targets, κ_fixed, b_w_int_target)` → free-entry calibration:
  drops `x_v`, inherits κ, targets b/w_int, leaves ϕ residual.

---

## Part 3 — Paper library utility

`Data/Key papers/pdf_to_markdown.py` mirrors every PDF in `Key papers/` as plain text under
`Key papers/markdown/`. It exists to make [`principles.md`](principles.md) **N11–N14**
enforceable: a claim about the literature needs a verbatim quote, and checking one by
scrolling a PDF is slow enough that it gets skipped.

```
cd "Data/Key papers" && python pdf_to_markdown.py     # new or changed PDFs only
grep -rn "product destruction" markdown/              # check a claim across the library
```

| Point | Detail |
|---|---|
| Requires | the `pdftotext` executable. Ships with Git for Windows (`mingw64/bin`), already on PATH here |
| Output | `Key papers/markdown/<slug>.md`, **gitignored** — derived from the committed PDFs, regenerates in seconds |
| Page markers | `<!-- page N -->` between pages, so a quote can be cited by page. This is the point of the format, not decoration |
| Ligatures | always expanded. BED-era PDFs store `ﬁ` as one glyph, so without this `grep finding` silently misses the word — the exact failure that makes you conclude a paper does not say something it says |
| Default is quote-safe | `-layout` keeps line breaks and columns, so a phrase spanning a line break will **not** grep. `--reflow` fixes phrase search, `--dehyphenate` joins split words. Both stamp a warning in the file header: **do not quote from output generated with either flag** |

Re-run after adding a paper. The PDFs remain the source of record; the Markdown is
disposable.
