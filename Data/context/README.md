# Context Folder — Index
**Last reorganized:** September 5, 2026 (branch `Organize_Project_State_Estimation`)
**Last consistency pass:** October 18, 2026 (branch `r8_timing_reposting`, with `costly_vacancy_reposting` merged in: Q1 steps 5–7, BK ε-threshold and `rem:dilution_variety` + App. K, R19 user-cost offset). Handout `session_handout_20261018.md`. Earlier: Oct 14 (Q0, D13, Q1 steps 1–4); Oct 13 (P8); Oct 11–12 (theory rebuild, X^c). Handouts before Oct 9 are in `archive/`.
spliced (`prop:fe` leads `sec:mechanism`), every theory proposition audited and corrected, R16
repaired, X^c made exact in draft and code, `ass:matching` tightened. Handout
`session_handout_20261012.md`. Handouts before Oct 9 moved to `archive/`.

Eleven files, each with one job, plus the dated session handouts. If you are about to write
something here, check this table first — most duplication in the past came from appending
status notes to whichever file was open rather than the file that owns the topic.

| File | Owns | Does **not** own |
|---|---|---|
| [`principles.md`](principles.md) | Standing rules (1–21) **and the non-negotiable rules N1–N15** that constrain every choice | Anything provisional — that is a decision |
| [`decisions.md`](decisions.md) | Open decisions D2, D3, D4, D10, D11, D12 (venue / theory-first pivot); settled D1, D5–D7, S1–S12 (S10 timing, S11 X^r timing, S12 structure of `sec:mechanism`), each with reasons | Tasks. A decision is a choice; a task is work |
| [`pending_tasks.md`](pending_tasks.md) | Actionable work, ordered; the critical path | Decisions, results, draft status |
| [`estimation_design.md`](estimation_design.md) | How Blocks M and B combine; what is missing to run it; build order | Empirical results |
| [`findings.md`](findings.md) | All empirical and model-mechanism results with numbers | Draft status, tasks |
| [`draft_status.md`](draft_status.md) | Section-by-section draft state; proposition inventory and proof structures; build instructions | Results, tasks |
| [`data_and_files.md`](data_and_files.md) | Data sources, FRED IDs, file locations, empirical targets | How the code works |
| [`parameters.md`](parameters.md) | All structural parameters: economic meaning, classification, identifying moments | Code implementation details |
| [`model_equations.md`](model_equations.md) | All 33 model equations, dated in the draft's timing, with the code's deviations and the R8 restructuring spec (Block 1 = full reposting; Block 2 = partial reposting) | Calibration values, results |
| [`nesting_map.md`](nesting_map.md) | How the baseline relates to DMP, CK, GS and BGM: restrictions, type of equivalence (exact, limit, allocations only), what survives and is lost, draft label, proof location, program check; cross-cutting caveats (e_t ≥ 0, σ = 0, CK parameterization vs CK model) | Proof status and audit history (`draft_status.md`), numbers (`findings.md`) |
| [`pipeline.md`](pipeline.md) | Both pipelines: script inventory, run order, code conventions | What the results were |

## Reading order when picking the project back up

1. **`session_handout_20261018.md`** — latest handout (Oct 18: Q1 steps 5–7, Test 6 and the
   user-cost offset R19, BK ε-threshold, Remark 3 + App. K; resume with the Q1 close-out, then
   R19 and Q2). Before it: `session_handout_20261014.md` (Q0, D13, Q1 steps 1–4) and
   `session_handout_20261012.md` (theory rebuild, X^c fix)
2. `../CLAUDE.md` — one-paragraph framing, quick facts, current priorities
3. **`pending_tasks.md`** — the critical path Q0–Q6 (D13): Q1 close-out, then Q2 (trivariate VAR)
   → Q3 → Q4 → Q5 → Q6; theory task R19
4. **`decisions.md`** — what is unsettled (D12 and D2 are the ones MS must decide), and the
   settled timing (S10–S12)
5. `draft_status.md` — before editing `Draft.tex` (build rules, proposition inventory)
6. `model_equations.md` — the model spec, if the task touches the code
7. `estimation_design.md` — if the task touches estimation
8. The rest as needed. `archive/` holds the Sept 21–Oct 1 handouts (history only)

## Conventions

- **Dates are absolute.** "Recently" and "currently" rot; "May 14, 2026" does not.
- **Cross-reference rather than duplicate.** If a number appears in two files, one of them
  should be a link.
- **Flag staleness inline**, next to the number, not in an appended note at the end of the
  file. A correction at the bottom does not stop anyone reading the wrong number at the top.
- **Draft labels in backticks** (`prop:fe`, `eq:posterior`, `tab:calib_targets`)
  so they can be grepped against `Draft.tex` (use `grep -a`: the file contains NUL bytes).

## Related documents outside this folder

| Path | Contents |
|---|---|
| `../Notes/delta_calibration_and_the_reposting_margin.md` | **Sept 23, 2026.** Why δ̄_e = 3.2 % against BGM/GS/Shao-Silos at 10 %; why raising δ̄_e cannot fix the Beveridge curve; the reposting rate λ and whether to estimate it. Recommends keeping [D1](decisions.md) settled |
| `../Notes/beveridge_free_entry_GS.tex` (+ `.pdf`) | **Oct 6–7, 2026. Superseded Oct 11 by the draft** (`lem:theta`, `prop:fe`, `app:proof_fe`, with translated notation). Keep for the full derivation of every log-linear equation in GS's dating. Numbers: `../Programs baseline/gs_free_entry_beveridge.jl` |
| `../Notes/predetermined_vacancies_and_the_CK_special_case.md` | **Rewritten Oct 6, 2026** (the Oct 5 threshold T was per unit level, not per 1 %, and is withdrawn). The δ-shock Beveridge response: the accounting identity in consistent units, the stock elasticity ξ·ē/v̄ at finite ξ, the user-cost channel at free entry (`../Programs baseline/beveridge_free_entry.jl`), and what CK actually show |
| `../Notes/role_of_ingredients.md` | Ingredient ablation design (R14). ⚠ §2 is wrong (variety is not what makes the δ sign survive free entry); rewrite before R14 |
| `../Notes/Baseline_Blanchard_Kahn.md` | BK conditions, full baseline (DS-CES, p_0 = 0). BK holds for ε > 1; endogenous exit strengthens it |
| `../Notes/AGS_Blanchard_Kahn.md` | BK conditions for the AGS σ = 0 case |
| `../../../Soft-Min-Cobb-Douglas-Matching-Function/` | MS's companion SMCD repo (separate git repo). `equilibria_convexity_check.jl` loads its `MatchingSpec` library. Paper: Silva (2026), SSRN 6956640, cited as a working paper |
| `../Inspecting_mechanism_setup.md` | Design rationale for mechanism Comparisons A–D |
| `../Story_1_vs_2.md` | Partial-R² test of the two readings of the zero β^LD vacancy coefficient |
| `../Bartek analysis/shock_specific_irf_procedure.md` | Original LP implementation procedure |
| `../Bartek analysis/report.md` | Standalone empirical framework report |
| `../Programs baseline/CALIBRATION_IMPROVEMENTS.md` | Record of the `calibrate_shares` refactor |
| `../project_brief_empirical_lp.md` | March 2026 handoff brief — **largely superseded** by `../CLAUDE.md` plus this folder; keep for provenance |
