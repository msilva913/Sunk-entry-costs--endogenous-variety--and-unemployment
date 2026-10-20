# Session Handout — October 20, 2026 (for the cloud session)
**Branch:** `costly_vacancy_reposting`. Work only here; the cloud default checkout is `main`,
which is far behind. `r8_timing_reposting` is an older pointer to the same line.
**Previous handout:** [`session_handout_20261019.md`](session_handout_20261019.md). Read its
§4–5 (Oct 20) for the full results. This file is the entry point.

> Written for a Claude Code cloud session with no conversation memory. It replaces the local
> Windows session of Oct 15–20.

---

## 0. Before any work

1. `git checkout costly_vacancy_reposting && git pull`.
2. Check that pushing works: `git push --dry-run origin HEAD`. A 403 means the Claude GitHub
   App lacks write access to `msilva913/Sunk-entry-costs--endogenous-variety--and-unemployment`.
   Stop and tell MS; otherwise the work is lost when the container is reclaimed.
3. Read `Data/CLAUDE.md`, then `context/principles.md` (N1–N15), then this file.

## 1. One-time environment setup (fresh Linux machine)

- **Python:** `conda env create -f Data/environment.yml` (env `econ2315`), then
  `conda activate econ2315`. Without conda, pip-install the same packages at the pinned minor
  versions.
- **Julia:** install Julia, then from `Data/Programs baseline/` run
  `julia --project=. setup_environment.jl` and `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
  The model needs SymPy via PyCall, so point PyCall at a Python that has `sympy`
  (`ENV["PYTHON"] = <that python>; Pkg.build("PyCall")`).
- **Not in git:** `Data/Bartek analysis/.env` (the FRED key). It is only needed to re-download
  data; `raw_data.pkl` and `data/cache/` are tracked. If needed, MS adds `FRED_API_KEY` as an
  environment secret.
- **Draft builds:** `pdflatex` only, three passes, **never `bibtex`**.
- The Windows paths in `pipeline.md` (`%LOCALAPPDATA%\anaconda3`) apply to the old machine only.

**Reproduction checks after setup** (do these before anything new):
- `python part6b_shock_persistence_cyclical.py` in `Data/Bartek analysis/`. Then
  `git diff --stat data/results/shock_calibration.csv` should show last-digit noise at most.
  Restore the regenerated outputs with `git checkout`.
- In `Data/Programs baseline/`, run R8_test_timing, R8_test_reposting, R8_test_steady_state,
  R8_step5_defaults, R8_step7_baseline, R8_test6_free_entry, bk_epsilon_threshold and
  sigma_determinacy (`julia --project=. <name>.jl`). Each `.txt` should match the committed
  one (the steady-state test may differ at 1e-15).

## 2. Where things stand

- **Q1 (R8 code restructure) is closed.** 41 equations; states [Ñ, Ẽ, Ṽ^u, Ṽ^s, z, δ, δl, s].
  δl is a bookkeeping state for the realized exit rate. Λ̄_r (`Λr_target`) replaces α in Θ_e.
  See `model_equations.md` "Code status".
- **Determinacy:** σ > σ* = 0.011 and ε > ε* = 2.40.
- **Q2 (shock processes) is in progress.** Results so far (findings.md "Q2"):
  - **Trivariate VAR:** the s innovation is nearly orthogonal to z and δ (< 2 %), so D14 is low
    stakes.
  - **Beveridge correlation:** cor(u,v) depends only on Λ̄_r. Even Λ̄_r = 0.10 gives −0.63, so
    the reposting margin alone cannot hit the data.
  - **Shock contributions:** in the model, s supplies 93–99 % of var(u); in the data SVAR, at
    most about 25 %. z is under-amplified 3–6×, and **δ barely moves u (> 20× short)**.
  - **ϕ is not an amplification lever** (Hagedorn–Manovskii cancellation).
- **Decisions on Oct 20 (MS), in `decisions.md`:**
  - **S13:** end every series at 2019Q4, cut *before* filtering.
  - **S14:** 1951–2019 for z, s and the labor-market moments; 1992–2019 for δ.
  - Never sweep η_L or ε; both are data-pinned.
- **Known stale items:** `SHOCKS_BASE` (in `run_solution_core.jl`) and `Q2_shock_processes_check.txt`
  predate S13. The old targets σ(u) 14.96 / cor(u,v) −0.804 included 2020–24; pre-pandemic
  (1951–2019) they are 12.48 / −0.925.

## 3. Next, in order

1. **Build the S14 shock processes** (new script in `Data/Bartek analysis/`, reusing part6b's
   functions):
   - (z, s) VAR(1) on HP-1600 log cycles over 1951Q1–2019Q4, Cholesky z → s.
   - Before 1992Q3, s = (τ − δ̄)/(1 − δ̄), with δ̄ the 1992–2019 BED mean; from 1992Q3 use the
     measured identity.
   - The δ equation on 1992Q3–2019Q4 with δ ordered last.
   - Off-diagonal terms excluded; convert to monthly as in part6b.
   - Output a `.txt` and `.csv` (N15).
2. **With MS's approval:** refresh `SHOCKS_BASE`, and the `SHK` copy in `sigma_determinacy.jl`,
   then rerun the R8 and `Q2_*.jl` programs.
3. **Regenerate the Block M targets on the S14 windows.** Formalize the data SVAR on
   (z, δ, s, u, v) as the data side of Q3 Block B′. Report FEVD bounds across recursive
   orderings and impulse responses with bootstrap bands.
4. **z-elasticity diagnostic:** compare the model's z-only d ln θ/d ln z with the closed form at
   our persistence, with the business block on and frozen. Candidates: low z persistence, or
   the variety block absorbing z.
5. **The δ-propagation gap:** δ alone moves u by about 0.3 %, against about 6–7 % implied by the
   data. Break it down into shock scale, jobs destroyed per exit, and the variety channel.
6. **R19 (theory, MS's decision)**, as before.

## 4. Standing rules (short list; full list in `principles.md` and `CLAUDE.md`)

- No `Draft.tex` edits without MS's go-ahead.
- N15: every draft number must come from a program.
- Never write bare "free entry" for vacancies; say "vacancy free entry".
- Never assign numbers to symbolic model names (ρ_z, s, b, δ, …), loop variables included.
- Write Julia and LaTeX edit scripts to files with the Write tool.
- Don't install packages beyond the environment specs without asking.
- Citations go through LibKey (DOIs) before they enter the draft.
- Commit often on `costly_vacancy_reposting` and push. The cloud container is not persistent.
