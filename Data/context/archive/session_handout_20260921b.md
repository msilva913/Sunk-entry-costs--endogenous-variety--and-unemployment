# Session Handout — September 21, 2026 (session B)
**Branch:** `Organize_Project_State_Estimation` · **Last commit:** 07dcd61

## What happened this session

### 1. ξ_inv sweep (entry cost convexity → δ→u peak horizon)

Swept ξ_inv ∈ {0.5, 1.0, 2.0, 3.0, 5.0, 8.0} at BED spec (dest_ann=0.0320) with
calibrated shock processes (ρ_s=0.874, σ_s=0.0854). **Result: peak never moves past h=2.**
Higher ξ_inv makes the IRF flatter (more persistent) but never hump-shaped. The u_h20/u_peak
ratio rises from 0.52 → 0.90 but the peak is always at impact.

### 2. Joint ξ_inv × ε sweep

Added ε ∈ {4.3, 3.0, 2.0, 1.5}. Key findings:
- ε = 3.0, ξ_inv = 8.0: peak at h=2, u_h20/u_peak = 0.98, max|eig(hx)| = 0.9997
  (near unit root, but still monotone decay)
- **ε < 3 causes Blanchard-Kahn violations** (indeterminacy) — the variety channel has a
  hard determinacy ceiling
- **Structural impossibility confirmed**: no feasible parameter combination within the
  current model generates a hump-shaped δ→u response

### 3. News shock idea — considered and rejected

Idea: if product line obsolescence is foreseen, a news shock to δ would produce hump-shaped
u. **Rejected** because the Bartik LP identifies the response to a *contemporaneous* shock
at h=0. The hump in the data is a propagation phenomenon, not an anticipation phenomenon.

### 4. Two paths forward identified

Either:
- **(a) New propagation mechanism** — financial accelerator, input-output networks, multi-period
  firm wind-down. Substantial model extension.
- **(b) Skepticism of the Bartik LP peak horizon** — the instrument may be persistent in the
  cross-section (industry composition changes slowly, industry-level death rates are serially
  correlated), so the LP IRF is a Wold representation that conflates propagation with
  cumulative serially-correlated shocks. **Testable**: add lagged B^δ_{i,t-p} controls.

Path (b) is testable with existing data and should be done first.

## What was updated

| File | What changed |
|------|-------------|
| `context/findings.md` | Added "ξ_inv × ε sweep" section with full tables, structural impossibility finding, and news shock rejection |
| `context/decisions.md` | Added **D10** (🔴): Is the LP peak horizon a propagation fact or a Bartik persistence artifact? PMW 2021 point, lagged-instrument test |
| `context/pending_tasks.md` | Updated M8 (diagnosis complete, gates on E8); added **E8** (🔴): Bartik instrument persistence test; updated s-shock calibration status |
| `context/parameters.md` | Updated ξ_inv description with sweep result; added ε determinacy ceiling to identification tensions |
| Memory | Added `project-hump-shape-structural.md` |

## Uncommitted files

- `run_xi_sweep.jl` — ξ_inv sweep script (scratch, can be deleted)
- `run_xi_eps_sweep.jl` — joint ξ_inv × ε sweep script (scratch, can be deleted)
- Context file updates listed above

## What to do next (priority order)

1. **E8: Bartik instrument persistence test.** In `part5_lp.py`, add lagged B^δ_{i,t-1},...,
   B^δ_{i,t-p} (p=4–8) as controls. Re-estimate δ→u and δ→v IRFs. Also report
   within-state autocorrelation of B^δ after time FEs. This is the highest-priority empirical
   task — it determines whether M8 is a real model problem or a measurement artifact.

2. **Step 0 cascade** (still pending from prior session). Set `dest_ann = 0.0320` in
   `steady_state.jl`, re-run Comparisons A–D, regenerate mechanism PDFs, update §5.3 and §5.2.

3. **Settle D2** (PATH A vs B) and **D3** (β̂ ↔ β(θ) scaling) — both gate the estimation
   build.

## Key numbers to remember

| Quantity | Value | Source |
|----------|-------|--------|
| ρ_s, σ_s (calibrated) | 0.874, 0.0854 | part6b AR(1) on s, HP-1600, 1992Q3–2019Q4 |
| ρ_z, σ_z | 0.902, 0.0092 | part6b Cholesky z-first |
| ρ_δ, σ_δ | 0.592, 0.0669 | part6b Cholesky z-first |
| dest_ann (D1) | 0.0320 | BED Deaths emp-weighted 1993–2019 |
| cor(u,v) at BED | +0.995 | D1/M5 diagnostic with calibrated s |
| σ(θ)/σ(LP) at BED | 5.50 | D1/M5 diagnostic with calibrated s |
| δ→u peak horizon | h=1–2 (model, any params) | ξ_inv × ε sweep |
| δ→u peak horizon | h=17–20 (data) | Bartik LP |
| BK ceiling for ε | ε ≥ 3 (approx) | Joint sweep |

## Open decisions

| # | Status | Summary |
|---|--------|---------|
| D1 | ✅ Settled | BED 0.0320, but not yet implemented in steady_state.jl |
| D2 | 🔴 Open | PATH A vs B — gates Θ_e contents |
| D3 | 🔴 Open | β̂ ↔ β(θ) scaling — gates Block B |
| D10 | 🔴 **NEW** | Is LP peak horizon propagation or persistence artifact? Gates M8 |
