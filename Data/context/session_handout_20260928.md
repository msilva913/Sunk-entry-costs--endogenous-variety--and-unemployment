# Session Handout — September 28, 2026
**Branch:** `costly_vacancy_reposting` (branched from `instruments_LP_coefficient`).

This session did two things: (1) finished the E8 Bartik-persistence robustness (p=16), and
(2) developed and settled the **formulation** of the reposting margin (D11), the structural
change needed to fit the Beveridge curve. Whether the margin enters the paper (D11) is still
the open decision; the mechanism is now pinned down.

## What to read first, in order

1. This handout.
2. [`decisions.md`](decisions.md) **D11** — the Sept-28 banner at the top of D11 has the settled formulation.
3. [`model_equations.md`](model_equations.md) — Block 1 (current model, all 33 eqs) and Block 2 (reposting modifications). **New file this session.**
4. [`../Notes/delta_calibration_and_the_reposting_margin.md`](../Notes/delta_calibration_and_the_reposting_margin.md) — full argument, §4.3a (inescapability), §7.3–§7.3b (formulation, rationale, one-time-vs-recurring verdict).
5. [`pending_tasks.md`](pending_tasks.md) R-series — the reposting task list (R7 is the next analytical step).

## The reposting margin — settled formulation

**Problem.** Model cor(u,v) = +0.995 vs data −0.804. The s shock under costless reposting is the
culprit; raising δ_e helps only weakly (+0.995 → +0.887 over the whole defensible range).

**Inescapable (verified, note §4.3a).** The reposting inflow ∂v_pre/∂s = (1−δ_e)(1−u) > 0 is an
accounting identity (Prop. 5 Part 1, Channel 1). With s data-calibrated and δ_e small, no
parameter fits the Beveridge curve under costless reposting. Some separations must fail to repost.

**Mechanism.** On separation, each vacated position draws a reactivation cost **α·χ, χ ~ F** (the
same distribution as firm exit — no new distribution). One-time keep-or-retire decision;
reactivate iff **(Q − K) ≈ Q ≥ α·χ**, giving reposting rate

```
Λ_r = F((Q − K)/α) ≈ F(Q/α) = min[(1 − p_0) + p_0·(Q/(α·f_m))^ψ, 1]
```

Key settled choices, each with a reason:
- **Per-position draws, not per-firm** → preserves the single firm size (per-firm draws break DS-CES symmetry; per-position + LLN keeps all firms symmetric). Note §7.3.
- **Draw from F, not a new distribution G** → F is the cost of *operating* capacity, G of *creating* it; reposting operates existing capacity. One new parameter α only. Note §7.3a.
- **Threshold = stock value Q − K ≈ Q, not flow K** → one-time reactivation with permanent retirement, not recurring maintenance. K would be degenerate and temporary. Note §7.3b.
- **Notation Λ_r** (capital), not λ (which is marginal utility of consumption, f[14]).

**Dynamic payoffs (to size by simulation, R9), beyond the recalibrated level:** amplification
(procyclical Λ_r + reposting option in the JCC) and **persistence** via stock depletion (retired
positions rebuild at full cost Q) — potentially relevant to the δ→u peak-horizon gap (M8).

## Equations that change (model_equations.md Block 2)

- **f[22]** vacancy LOM: `v_pret' = (1 − δ_e)·[(1 − q)v + Λ_r·s·sbar·(1 − u)]` (factors through δ_e because reposting ⊂ survival).
- **f[16]** resource constraint: add `X_r` (aggregate reactivation costs paid). ⚠ exact conditional-mean coefficient to match X_c convention.
- **f[new]** reposting rate `Λ_r = F(Q/α)`.
- **f[3] JCC** ⚠ **to be derived (R7)** — separation branch gains `Λ_r'·(vacancy value) − expected reactivation cost`. **This is the next analytical step.**
- Unchanged: f[2] δ_e, f[23] u LOM, f[24] N LOM (reposting destroys positions, not workers' flows or variety).

## Next steps (priority)

1. **R7 — derive the JCC** with the reposting option (Bellman rebuild). Blocks coding. See model_equations.md f[3].
2. **R1 — decide D11** (does the margin enter the paper).
3. **R2 — identification check** on simulated data: α vs σ_s vs δ_e vs f_m.
4. **R8 — implement** α, Λ_r, f[22], f[16], f[3] in `steady_state.jl` and `run_solution_core.jl`.
5. **R9 — simulate** to size amplification and persistence; ties to E9 (run the LP on model-simulated data) and M8.

Separately, still owed from before: **D1 cascade** (set `dest_ann = 0.0320` in `steady_state.jl:613`, re-run Comparisons A–D, regenerate mechanism PDFs, update §5.3/§5.2 and the four `δ_e/τ ≈ 0.21` claims). And **D10/E9** (model-side LP) for the Block B IRF target.

## Files changed this session

- **New:** `context/model_equations.md`, `context/session_handout_20260928.md`, memory `project-reposting-margin.md`.
- **Edited:** `Notes/delta_calibration_and_the_reposting_margin.md` (§4.3a added; §7 rewritten to per-position, Λ_r, §7.3b added), `context/decisions.md` (D11 Sept-28 banner), `context/parameters.md` (α, Λ_r), `context/pending_tasks.md` (R-series updated, R7–R9 added), `context/findings.md` (reposting mechanism section), `context/estimation_design.md` (α in Θ_e), `context/pipeline.md` (model_equations pointer), `context/README.md`, `CLAUDE.md`, `Draft/Draft.tex` (graphicspath fixes), `Bartek analysis/part5_lp.py` (E8 p=16).

## Key numbers

| Quantity | Value |
|---|---|
| model cor(u,v) at BED δ | +0.995 |
| data cor(u,v) | −0.804 |
| cor(u,v) with s silenced | −0.911 |
| cor(u,v) at δ_e = τ (extreme) | +0.887 |
| Λ_r bounds | [1 − p_0, 1] |
| reposting bites for | α > Q̄/f_m |
| K / Q (monthly) | ≈ 0.6% (so Q − K ≈ Q) |
