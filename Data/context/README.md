# Context Folder — Index
**Last reorganized:** September 5, 2026 (branch `Organize_Project_State_Estimation`)
**Last consistency pass:** October 7, 2026 (branch `costly_vacancy_reposting`; R15 free-entry proposition restated in R form, `.tex` note authoritative, `.md` exposition deleted). Previous: October 5 (intro positioning, `app:loglin` ¶9 rewrite, predetermined-vacancy threshold; R12b/R12c added)

Nine files, each with one job, plus the dated session handouts. If you are about to write
something here, check this table first — most duplication in the past came from appending
status notes to whichever file was open rather than the file that owns the topic.

| File | Owns | Does **not** own |
|---|---|---|
| [`principles.md`](principles.md) | Standing rules (1–21) **and the non-negotiable rules N1–N15** that constrain every choice | Anything provisional — that is a decision |
| [`decisions.md`](decisions.md) | Open decisions D2, D3, D4, D10, D11; settled D1, D5–D7, S1–S12 (S10 timing, S11 X^r timing, S12 Prop. 5 Part 1), each with reasons | Tasks. A decision is a choice; a task is work |
| [`pending_tasks.md`](pending_tasks.md) | Actionable work, ordered; the critical path | Decisions, results, draft status |
| [`estimation_design.md`](estimation_design.md) | How Blocks M and B combine; what is missing to run it; build order | Empirical results |
| [`findings.md`](findings.md) | All empirical and model-mechanism results with numbers | Draft status, tasks |
| [`draft_status.md`](draft_status.md) | Section-by-section draft state; proposition inventory and proof structures | Results, tasks |
| [`data_and_files.md`](data_and_files.md) | Data sources, FRED IDs, file locations, empirical targets | How the code works |
| [`parameters.md`](parameters.md) | All structural parameters: economic meaning, classification, identifying moments | Code implementation details |
| [`model_equations.md`](model_equations.md) | All 33 model equations, dated in the draft's timing, with the code's deviations and the R8 restructuring spec (Block 1 = full reposting; Block 2 = partial reposting) | Calibration values, results |
| [`pipeline.md`](pipeline.md) | Both pipelines: script inventory, run order, code conventions | What the results were |

## Reading order when picking the project back up

1. **`session_handout_20261001.md`** — latest handout. Read its Oct 5 update first; §0 (uncommitted work) and the "push local commits" note are **resolved** as of Oct 5, 2026, and R13 is **done**
2. `../CLAUDE.md` — one-paragraph framing and current priorities
3. **`pending_tasks.md`** — the critical path (as of Oct 5: **E9**, then the single regeneration pass **R8 + D1 + PATH B**; R13 complete)
4. **`decisions.md`** — what is unsettled, and the settled timing (S10–S12)
5. `model_equations.md` — the model spec, if the task touches the code
6. `estimation_design.md` — if the task touches estimation
7. The rest as needed. Older dated `session_handout_*.md` files are history

## Conventions

- **Dates are absolute.** "Recently" and "currently" rot; "May 14, 2026" does not.
- **Cross-reference rather than duplicate.** If a number appears in two files, one of them
  should be a link.
- **Flag staleness inline**, next to the number, not in an appended note at the end of the
  file. A correction at the bottom does not stop anyone reading the wrong number at the top.
- **Draft labels in backticks** (`prop:ds_asymmetry`, `eq:posterior`, `tab:calib_targets`)
  so they can be grepped against `Draft.tex`.

## Related documents outside this folder

| Path | Contents |
|---|---|
| `../Notes/delta_calibration_and_the_reposting_margin.md` | **Sept 23, 2026.** Why δ̄_e = 3.2 % against BGM/GS/Shao-Silos at 10 %; why raising δ̄_e cannot fix the Beveridge curve; the reposting rate λ and whether to estimate it. Recommends keeping [D1](decisions.md) settled |
| `../Notes/beveridge_free_entry_GS.tex` (+ `.pdf`) | **Oct 6–7, 2026. Authoritative source for the R15 proposition** (the Oct 6 `.md` exposition was deleted Oct 7). Draft layout (elsarticle, review). AGS with σ = κ = 0 (= GS by `rem:nesting`), free entry, log shock δ_t = δ̄·exp(δ̃_t). Nonlinear system in GS's form and dating, steady state as GS Def. 2, every log-linear equation derived in full with its weights interpreted. `lem:theta` (tightness rule in the shock alone, with the determinacy condition ϕ(1−δ̄)f̄ < η_L(2+r−τ̄)), `rem:recursive`, `prop:fe` (conditions ρ > ω + R and ρ(1−λ²)/(1+ρλ) > ω + R), `cor:allh`. Not in the draft |
| `../Notes/predetermined_vacancies_and_the_CK_special_case.md` | **Rewritten Oct 6, 2026** (the Oct 5 threshold T was per unit level, not per 1 %, and is withdrawn). The Beveridge response to a δ shock: the accounting identity in consistent units, the effective stock elasticity ξ·ē/v̄ at finite ξ, the user-cost channel at free entry (`../Programs baseline/beveridge_free_entry.jl`), and what CK actually show. Oct 5 text: Derives the threshold entry elasticity T = κ_v/(δ̄_e/τ̄) − 1 from `eq:v_ll`; shows δ_e = τ is the boundary of the parameter space and **minimizes** T. Numbers from `../Programs baseline/predetermined_vacancy_share.jl`. Records two draft bugs (fixed) and the Prop. 5 proof bug (R12c) |
| `../Notes/Baseline_Blanchard_Kahn.md` | BK conditions, full baseline (DS-CES, p_0 = 0). BK holds for ε > 1; endogenous exit strengthens it |
| `../Notes/AGS_Blanchard_Kahn.md` | BK conditions for the AGS σ = 0 case |
| `../Inspecting_mechanism_setup.md` | Design rationale for mechanism Comparisons A–D |
| `../Story_1_vs_2.md` | Partial-R² test of the two readings of the zero β^LD vacancy coefficient |
| `../Bartek analysis/shock_specific_irf_procedure.md` | Original LP implementation procedure |
| `../Bartek analysis/report.md` | Standalone empirical framework report |
| `../Programs baseline/CALIBRATION_IMPROVEMENTS.md` | Record of the `calibrate_shares` refactor |
| `../project_brief_empirical_lp.md` | March 2026 handoff brief — **largely superseded** by `../CLAUDE.md` plus this folder; keep for provenance |
