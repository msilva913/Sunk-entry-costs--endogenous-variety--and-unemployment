# Context Folder — Index
**Last reorganized:** September 5, 2026 (branch `Organize_Project_State_Estimation`)
**Last consistency pass:** September 23, 2026 (branch `instruments_LP_coefficient`)

Nine files, each with one job, plus the dated session handouts. If you are about to write
something here, check this table first — most duplication in the past came from appending
status notes to whichever file was open rather than the file that owns the topic.

| File | Owns | Does **not** own |
|---|---|---|
| [`principles.md`](principles.md) | Standing rules (1–21) **and the non-negotiable rules N1–N15** that constrain every choice | Anything provisional — that is a decision |
| [`decisions.md`](decisions.md) | Decisions D1–D11 (D1 settled; D2, D3, D11 open; D10 partial) + settled S1–S9, each with reasons | Tasks. A decision is a choice; a task is work |
| [`pending_tasks.md`](pending_tasks.md) | Actionable work, ordered; the critical path | Decisions, results, draft status |
| [`estimation_design.md`](estimation_design.md) | How Blocks M and B combine; what is missing to run it; build order | Empirical results |
| [`findings.md`](findings.md) | All empirical and model-mechanism results with numbers | Draft status, tasks |
| [`draft_status.md`](draft_status.md) | Section-by-section draft state; proposition inventory and proof structures | Results, tasks |
| [`data_and_files.md`](data_and_files.md) | Data sources, FRED IDs, file locations, empirical targets | How the code works |
| [`parameters.md`](parameters.md) | All structural parameters: economic meaning, classification, identifying moments | Code implementation details |
| [`pipeline.md`](pipeline.md) | Both pipelines: script inventory, run order, code conventions | What the results were |

## Reading order when picking the project back up

1. **`session_handout_20260923.md`** — latest session handout (E8 persistence test results)
2. `../CLAUDE.md` — one-paragraph framing and current priorities
3. **`decisions.md`** — what is unsettled; this is where the branch's work is
4. `pending_tasks.md` — the critical path
5. `estimation_design.md` — if the task touches estimation
6. The rest as needed

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
| `../Notes/Baseline_Blanchard_Kahn.md` | BK conditions, full baseline (DS-CES, p_0 = 0). BK holds for ε > 1; endogenous exit strengthens it |
| `../Notes/AGS_Blanchard_Kahn.md` | BK conditions for the AGS σ = 0 case |
| `../Inspecting_mechanism_setup.md` | Design rationale for mechanism Comparisons A–D |
| `../Story_1_vs_2.md` | Partial-R² test of the two readings of the zero β^LD vacancy coefficient |
| `../Bartek analysis/shock_specific_irf_procedure.md` | Original LP implementation procedure |
| `../Bartek analysis/report.md` | Standalone empirical framework report |
| `../Programs baseline/CALIBRATION_IMPROVEMENTS.md` | Record of the `calibrate_shares` refactor |
| `../project_brief_empirical_lp.md` | March 2026 handoff brief — **largely superseded** by `../CLAUDE.md` plus this folder; keep for provenance |
