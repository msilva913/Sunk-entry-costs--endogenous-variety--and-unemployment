# Session Handout — October 19, 2026
**Branch to resume on:** `costly_vacancy_reposting` (fast-forwarded to `r8_timing_reposting`;
both point at the same commit). **Last commit:** see `git log -1`.
**Previous handout:** [`session_handout_20261018.md`](session_handout_20261018.md)

> Written to be resumed on another machine with no conversation memory.
> ⚠ The agent shell has no GitHub credentials. **Push `costly_vacancy_reposting` and
> `r8_timing_reposting` from your Git client before leaving this machine.**

---

## 0. Branches

| Branch | State |
|---|---|
| `costly_vacancy_reposting` | **Work here.** Everything: draft, Q1 code, context files |
| `r8_timing_reposting` | Same commit as above. Can be deleted after pushing |
| `r8_oct15_local` | Local backup of the Oct 15 work before reconciliation (`a9b253c`). Delete once satisfied |

## 1. What was done (Oct 15–19)

**Reconciliation (Oct 19, `853b053`).** The Oct 15 work on this machine and the Oct 18 work on
the other machine had diverged. Local r8 was reset to origin and the Oct 15 work re-applied on top:
- `run_solution_core.jl`: cleaned layout on the remote equation numbering. Comments keep the
  intuition. δl = δ_{t−1} is documented as a **bookkeeping state**: it enters only f[2], the
  realized measured exit rate δ_{e,t} = 1 − (1−δ_{t−1})Λ_t. Forward-dating δ_e to drop it gives
  the expected rate and breaks Test 5 at first order. `rp` is the reposting switch (1 partial,
  0 costless; α = 0 would divide by zero, Λ̄_r = 1 is a kink).
- `R8_test_timing.jl`: built PAR by hand (23 entries, stale since α and rp were added); now
  `build_PAR`. Test 1 passes (1.9e-12).
- `Project.toml`/`Manifest.toml`: JLD2 restored (needed by `run_solution_endog_exit`).
- **Verification:** all 8 R8 programs reproduce their committed outputs (steady-state test to
  1e-15 residual noise); all 10 runners pass.

**New result: determinacy in σ** (`sigma_determinacy.jl/.txt`; findings.md "Determinacy in σ").

| Configuration | σ* |
|---|---|
| Baseline (p_0 = 0.5, variety) | **0.0108** (no bounded equilibrium below) |
| p_0 = 0, variety | 0.0058 |
| p_0 = 0, ρ ≡ 1 | 0.0029 (indeterminate below) |
| as above, business block frozen | none: determinate at every σ, including 0 |

Source: product-line entry under near-linear utility. Free entry fixes ν_f and the entry Euler
pins only the expected dividend, so the timing of entry is a free direction. DMP and CK are
stabilized by congestion q(θ) and convex G; constant-cost entry has no such curvature. σ > 0
adds consumption smoothing, which removes the free direction. Same mechanism as R19 (the entry
Euler pins m(1−δ_t)). **Use:** a justification for σ > 0, and a footnote on `prop:independence`.

**Q1 closed (`f0694e4`), merged into `costly_vacancy_reposting`.**
- `model_equations.md`: "Code status" resolved; states, controls, and the map from the 41
  equations to the code's f[k].
- `parameters.md`: **Λ̄_r (`Λr_target`) replaces α in Θ_e**; α backed out (Θ_d).
- `pipeline.md`: the conventions section updated for R8 (states, equation order, `TARGETS`,
  `build_PAR` only).
- `decisions.md` D11: implementation settled; go/no-go open until Q2 and R2/Q4.

**Q2 scoped; D14 opened** (no code yet). See [`decisions.md`](decisions.md) D14 and the Q2 row
of [`pending_tasks.md`](pending_tasks.md). What is missing from the part6b approach
((z, δ) VAR + univariate AR(1) on s = (τ−δ)/(1−δ)):
1. σ_s is never orthogonalized against the z and δ innovations, yet the model treats them as
   independent. An overstated independent σ_s pushes cor(u,v) toward +1.
2. The double-counting argument used for z → δ does not apply to s: the model has no endogenous
   s. Orthogonalizing therefore drops real countercyclical separations (D14, open).
3. BED noise enters s with weight ≈ δ/τ ≈ 0.09, which builds a spurious negative s–δ correlation.
4. ρ_s absorbs any lagged z, δ → s feedback; the cross terms are never reported.
5. The shock processes must share the identification of Block B′ (D13).

## 2. Next steps, in order

1. **Q2: trivariate (z, δ, s) VAR.** Extend `Bartek analysis/part6b_shock_persistence_cyclical.py`:
   VAR(1), z → δ → s, HP-1600, 1992Q3–2019Q4. Report A₁ with t-stats, the residual correlations,
   σ_s raw vs orthogonalized, the share of s variance explained by z and δ, and monthly (ρ, σ).
   Write a `.txt` for §5.2 (N15). ⚠ Needs the **`econ2315` conda environment** (pandas,
   statsmodels). The Oct 19 machine has no Python with them: PyCall's Conda Python lacks pandas,
   and nothing was installed.
2. **D14 (MS).** Decide on orthogonalizing s using the Q2 numbers. Recommendation: orthogonalize
   in the baseline, correlated innovations as robustness, with the σ(τ) shortfall reported.
3. **Model check.** Update `SHOCKS_BASE` in `run_solution_core.jl` (and `SHK` in
   `sigma_determinacy.jl`, which keeps its own copy), then rerun `R8_step7_baseline.jl`: the
   cor(u,v) grid over Λ̄_r plus σ(τ). Before: +0.995 / +0.968 / +0.812 / +0.461 at costless /
   0.95 / 0.85 / 0.70, against −0.804.
4. **R19 (theory, MS).** Pricing replacement offsets `prop:fe`'s user-cost channel in the
   baseline. Double-limit check, then a remark after `prop:fe`, then the framing. The σ* result
   fits here as a footnote.
5. Then Q3 (targets), as in `pending_tasks.md`.

**Optional:** a convex entry-cost test of the determinacy intuition. With G convex in entry,
σ* should fall toward 0.

## 3. Standing constraints (unchanged)

- `pdflatex` only, three passes, never `bibtex`. Do not edit `Draft.tex` without MS's go-ahead.
- N15: every draft number comes from a program.
- Never write bare "free entry" for vacancies: "vacancy free entry" or "infinitely elastic vacancy creation".
- Never assign numbers to symbolic model names (ρ_z, s, b, δ, …), including loop variables.
- Write Julia/LaTeX edit scripts to files; the Bash tool collapses backslashes.
- Don't install packages without asking.
