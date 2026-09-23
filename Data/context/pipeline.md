# Pipeline — Empirical and Model Code
**Last updated:** September 23, 2026 (added the ξ sweeps and E8; corrected the
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
| `part2b_residualize_shocks_v3.py` | **v3**: adds monetary-policy sensitivity φ_j × ΔMP_t | ⚠️ written, **never run** — see [D8](decisions.md) |
| `part2c_granger_lp.py` | Panel LP Granger causality, δ vs. LD | ✅ |
| `part2d_var_calibration.py` | Bivariate panel VAR(1) on (ν^δ, ν^LD) | ✅ |
| `part3_instrument.py` / `part3_instrument_s.py` | Raw (non-residualized) Bartik aggregation | ✅ appendix only |
| `part5_joint_lp.py` | Joint LP: δ and LD in the same regression | ✅ robustness, not structural ([D6](decisions.md)) |
| `part5_joint_diagnostics.py` | Partial-R²/first-stage test of Story 1 vs. Story 2 for the zero β^LD | ✅ |
| `part6_shock_persistence.py` | AR(1) persistence of the raw and residualized shocks | ✅ |
| `part6b_shock_persistence_cyclical.py` | **Model shock calibration**: bivariate VAR(1) on HP-1600 log z, log δ, Cholesky z-first | ✅ → ρ_z^m, ρ_δ^m, σ_z^m, σ_δ^m |
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
| `run_solution_core.jl` | The 33-equation model system, symbolic parameters/variables, and `solution_interface`. Audited May 20, 2026. **`:97`** hard-errors when the steady-state residual exceeds tolerance — never downgrade it to a warning. **`:411` `SS_numeric`** now calls `steady_state(cal)` and asserts `PAR` consistency at `:455`. ⚠️ **`:498` `SS_symbolics` still carries the original gross-Lerner-share defect** — latent only because every runner passes `SS_precomputed`. Never use its output as a linearization point. `:131` still reads `estimate = []` |
| `run_solution.jl` | Thin baseline runner; serializes `model_output.jls`. Carries **placeholder** ρ_s = 0.90, σ_s = 0.010. |
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
carrying hard-coded figures. Wiring that up is the remaining half of M9.

### Calibration diagnostics (not paper comparisons)

| Runner | Purpose | Output |
|---|---|---|
| `run_solution_delta_target.jl` | **D1/M5 diagnostic** — dest_ann ∈ {0.0320 BED, 0.0754 code, 0.0963 BGM}; reports steady state, Block M moments at HP λ=1,600, and δ→u/δ→v IRF persistence in quarters. Decided that δ_e can be fixed. ✅ **run September 6, 2026** (after the `SS_numeric`/`hp_filter` fixes) and **re-run September 21, 2026** with the calibrated s process. Results in [`findings.md`](findings.md) §D1/M5 | `irf_delta_target.jls` |

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

### Comparison / legacy variants

`run_solution_CK.jl` and `steady_state_CK.jl` (Coles-Kelishomi benchmark);
`run_solution_general_CES.jl` (general CES — **loads the old Matlab `posterior_mode.mat`**,
which belongs to a previous generation of the model and must not be reused);
`GS.jl` (Gabrovski-Silva JEDC benchmark); `run_solution_{alt,elastic,high_b,kappa,phi,
risk_neutral,simplified}.jl` (single-parameter sensitivity runs).

### Diagnostics and output

`steady_state_checks.jl` (calibration inspection; free-entry branch bug — [D9](decisions.md)),
`steady_state_checks_*.jl`, `accuracy_checks.jl`, `comparative_statics.jl`,
`export_endog_exit_irfs.jl` (IRFs + SS → CSV), `generate_tables.jl`,
`second_moments.jl` (⚠️ **stale and invalid**: λ=100,000, only 8 series, AND it predates the `hp_filter` fix of Sept 5, 2026 — every moment it ever produced is wrong; superseded by `simulated_moments` inside `run_solution_delta_target.jl`), `Hill estimator/` (tail-index estimation for the cost distribution).

### Model code conventions (audited May 20, 2026)

- **Naming** — `r` is the discount rate; `ρ` is the relative price N^(1/(ε−1)). Consistent
  with `steady_state.jl`.
- **κ** — matching cost paid *per match* (Pissarides 2009). The resource constraint uses
  `κ·q·v`, not `κ·v`.
- **u LOM** — uses total `v_t = v_pret + e_t`; entrants participate in matching within
  period t (draft `eq:v_lom`).
- **States** — `x = [u, N, v_pret, z, δ, s]`; 31 controls; monthly frequency.
- **Equation order** — f[1] exit threshold χ^c, f[2] δ_e, f[3] JCC, f[4] business-formation
  Euler (BGM), f[5] K value, f[6] MRP, f[7] Nash wage, f[8] vacancy creation, f[9]–f[30]
  static conditions / LOMs / observables, f[31]–f[33] AR(1) shocks.
- **SS_symbolics** — `x_c_ss = Y_c·(μ−1)/(μ·N) + ν_f` (consistent with f[1]); `C_ss`
  includes `X_c`; `Y_ss = C + ν_f·N_e`.
- **Known toolkit bug** — `eval_SS` iterates `for ip in npar` (once) instead of `1:npar`,
  so it substitutes only the last parameter. Callers must pass a pre-computed Float64 SS to
  `solution_interface`; `run_solution.jl` uses the slow symbolic fallback instead.

### `steady_state.jl` architecture (as of May 31, 2026)

- `calibrate_shares(targets)` → 18-parameter NamedTuple. **PATH A** (default; targets
  `Xc_Y`) or **PATH B** (triggered by supplying `dest_elast_target`; targets the exit
  elasticity, making `Xc_Y` an outcome). ϕ is a residual from Nash in Stage 5.
- `steady_state(para; init=0.51)` → 2×2 solver in (log θ, log χ^c).
- `steady_state_free_entry(para)` → 1D solver for ξ_inv = 0; outer root-find over θ (JCC),
  inner over δ_e. **Currently converges to the wrong (high-θ) branch** — [D9](decisions.md).
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
