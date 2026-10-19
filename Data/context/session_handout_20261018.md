# Session Handout — October 18, 2026
**Branch to resume on:** `r8_timing_reposting` (now also carries the latest draft:
`costly_vacancy_reposting` was merged in at `1a25965`). **Last commit:** see `git log -1`.
**Previous handout:** [`session_handout_20261014.md`](session_handout_20261014.md)

> Written to be resumed on another machine with no conversation memory.
> ⚠ The agent shell has no GitHub credentials. **Push `r8_timing_reposting` from your Git
> client before leaving this machine** (it was 4+ commits ahead of origin when this was written).

---

## 0. Branches

| Branch | State |
|---|---|
| `r8_timing_reposting` | Q1 code (steps 1–7), the BK programs, **and** the current draft (merge `1a25965`). Work here |
| `costly_vacancy_reposting` | Draft line through `ac78fd0` (Remark 3 + App. K). Behind r8 by the code commits |

Close-out of Q1 ends with merging `r8_timing_reposting` → `costly_vacancy_reposting`. Since the
reverse merge already happened, this is a fast-forward or trivial merge.

## 1. What was done (Oct 18)

**Q1 steps 5–7, all tests passing** (plan: [`../Notes/Q1_R8_D1_PATHB_plan.md`](../Notes/Q1_R8_D1_PATHB_plan.md)):
- **Step 5 (defaults).** `TARGETS` in `steady_state.jl`:
  - `dest_ann = 0.0320` (D1) and `dest_elast_target = 5.0` (PATH B, a placeholder);
  - `b_ratio = 0.9`, `x_v = 0.5`, `Λr_target = nothing` (costless).
  - Result: δ_e/τ = 0.0873; SS residual < 1e-12 for costless, the Λ̄_r grid and p_0 = 0 (`R8_step5_defaults.txt`).
  - `steady_state_checks.jl` reads `TARGETS` and checks the reposting identities. Its free-entry tail still fails (M2, pre-existing).
  - **D2 settled: PATH B.**
- **Step 6 (runners).**
  - Core helpers: `SHOCKS_BASE` (part6b shock processes until Q2), `shock_col(state, sd)`, and a qzflag guard in `solution_interface`.
  - Ten runners now use `build_PAR` and named shock columns, and all run. The Prop. 5 checkers reproduce their committed CSVs to 1e-11.
  - Comparison A's CK spec uses δ_e = τ(1−1e-6), because Ṽ^s is now a log state.
  - Comparison B's re-specification is documented for Q5.
  - Five legacy sensitivity runners are marked unmaintained.
- **Step 7 (baseline over Λ̄_r)**, `R8_step7_baseline.txt`, exploratory:

  | Λ̄_r | cor(u,v) | s only | δ only | z only |
  |---|---|---|---|---|
  | costless | +0.995 | +1.000 | −0.876 | −0.977 |
  | 0.95 | +0.968 | +0.974 | −0.869 | −0.961 |
  | 0.85 | +0.812 | +0.823 | −0.858 | −0.944 |
  | 0.70 | +0.461 | +0.478 | −0.841 | −0.932 |

  Data: −0.804. The s shock sets the sign. The s process is still the AR(1) placeholder, so
  this is a reason to do Q2 before drawing conclusions about the Λ̄_r required.
- **Test 6** (`R8_test6_free_entry.txt`): passes under **AGS emulation** (ρ ≡ 1, clone
  replacement via Ñ_{t+1} = 1 in place of the entry Euler, σ = κ = p_0 = 0). The gap to
  c_δ, c_s falls tenfold per decade of ξ_inv, to 1.9e-4 at 1e-8.

**Key finding: R19 (needs MS).** With product-line entry, the entry Euler pins the
survival-adjusted discount factor m(1−δ_t). The vacancy user cost K_t = Q_t − m(1−δ_t)Q_{t+1}
uses the same factor, so a δ shock barely moves K or θ: θ̂/δ̃ ≈ −0.003, against c_δ = −0.72
in AGS. `prop:fe` is correct for AGS, but its mechanism is nearly switched off once business
formation is priced. The baseline's δ-sign must come from the finite-ξ predetermined stock and
variety. Details: [`findings.md`](findings.md) "Q1 regeneration results";
[`pending_tasks.md`](pending_tasks.md) R19.

**Determinacy in ε**:
- Programs: `bk_epsilon_threshold.jl`, `lrbgm_determinacy.jl`, `cutoff_profit_share.jl`.
- No bounded equilibrium below **ε* = 2.40** at the defaults. The threshold is 1.40 with exit
  off, 2.12 at b/w = 0.71, and there is none without variety. The LR-BGM block is determinate
  for ε ∈ [1.05, 8].
- Mechanism: costs paid in goods are not neutralized by ρ(N). The survival cutoff has
  elasticity 1/(ε−1) − s_R in N, with s_R = 0.013.
- Corrected records: `findings.md` (was "ε < 3, indeterminacy"), S9, and
  `Notes/Baseline_Blanchard_Kahn.md` (its "BK for any ε > 1" holds for LR-BGM only).

**Draft**:
- `rem:dilution_variety` (Remark 3, after `rem:bgm_nest`), shortened by MS.
- New **Appendix K** `app:determinacy` with `tab:eps_threshold`; Notation is now Appendix L.
- The g_N paragraph after `prop:exit_s` now points to the remark.
- Correction made: the log-utility cancellation needs X_t = 0.

## 2. Next steps (in order)

1. **Close Q1** (pending_tasks Q1, items i–v):
   - `model_equations.md` "Code status" → resolved, with the new state vector and 41 equations;
   - `parameters.md`: Λ̄_r replaces α in Θ_e;
   - `pipeline.md` conventions (the Q1 program table is already there);
   - final D11 record;
   - merge into `costly_vacancy_reposting`.
2. **R19, MS to decide the order.** Recommended:
   - (a) a double-limit check (does the offset survive ε → ∞ with f_e = f̄_e/ε, where the entry
     Euler degenerates?);
   - (b) a short remark after `prop:fe`;
   - (c) revisit how `sec:mechanism` and D12 frame `prop:fe`.
3. **Q2:** the trivariate (z, δ, s) VAR in `part6b_shock_persistence_cyclical.py` (Python). Then
   Q3 → Q4 → Q5 (agree the simple-model-first threshold before Q5) → Q6.

## 3. Environment and gotchas

- **Julia:**
  - Run from `Data/Programs baseline` with `julia --startup-file=no <file>.jl`.
  - The pinned `--project=.` environment needs `Pkg.instantiate()` first; on this machine it
    lacked `Preferences`, and the global environment works.
  - Long runs (sweeps, bisection, step 7) take 10–30 minutes; run them in the background.
- **Never assign numbers to the symbolic model names** (`ρ_z`, `s`, `b`, …), even inside a
  REPL loop: VS Code's soft scope rebinds the global. That was the Oct 18 `DataFrame * Sym`
  error at `f[7]`, fixed in the R8 tests.
- **Shell:**
  - The Bash tool collapses backslashes in heredocs and inline `-e`. Write edit scripts with
    the Write tool, using `raw""` strings.
  - `@sprintf` needs a single literal format string.
  - `include` inside a scratch script resolves relative to the script's location, so use
    `joinpath(pwd(), …)`.
  - Windows path length: git worktrees need a short path such as `GitHub/sec-wt`.
- **LaTeX:**
  - Three `pdflatex` passes from `Data/Draft/`; **never `bibtex`**.
  - Close TeXstudio before scripted edits, since a held file makes git and writes fail.
  - Expect two undefined references: `app:robustness`, `app:weighting_robustness`.
- **Owed (MS):** add `silva_2026` and `den_haan_2000` to `references.bib`.

## 4. Reading order
1. This handout.
2. `../CLAUDE.md` "Priority tasks right now".
3. [`pending_tasks.md`](pending_tasks.md): Q1 row, R19, Q2.
4. [`findings.md`](findings.md): the last section.
5. [`../Notes/Q1_R8_D1_PATHB_plan.md`](../Notes/Q1_R8_D1_PATHB_plan.md), for the close-out list.
