# Session Handout — September 22–23, 2026
**Branch:** `instruments_LP_coefficient` · **Previous branch:** `Organize_Project_State_Estimation`

## What happened this session

### E8: Bartik instrument persistence test (completed)

Implemented and ran the augmented LP with lagged instruments in `part5_lp.py` section [11].
This is the test prescribed by D10 to determine whether the baseline LP's δ→u peak at
h=17–20 is genuine propagation or a Wold-representation artifact from instrument persistence.

**Within-state autocorrelation of B̃^δ (after time FEs):**
- Median ρ = 0.91 (lag 1), 0.83 (lag 4), 0.79 (lag 8), 0.76 (lag 12)
- 90–100% of states significant at every lag — the instrument is near-unit-root

**Augmented LP results (p = 4, 8, 12, 16 lagged instruments):**

δ→u:

| Spec | Peak h | Peak β (pp) | p |
|------|--------|-------------|---|
| baseline | 17 | 1.71 | 0.001 |
| p=4 | 10 | 0.95 | 0.004 |
| p=8 | 10 | 1.07 | 0.021 |
| p=12 | 10 | 1.35 | 0.015 |
| p=16 | 13 | 1.86 | 0.011 |

δ→v:

| Spec | Peak h | Peak β (pp) | p |
|------|--------|-------------|---|
| baseline | 18 | −0.58 | 0.001 |
| p=4 | 5 | −0.55 | 0.001 |
| p=8 | 18 | −0.57 | 0.001 |
| p=12 | 5 | −0.57 | 0.007 |
| p=16 | 5 | −0.65 | 0.001 |

### Key finding: the peak horizon is not well-identified

The δ→u peak drifts from h=10 to h=13 and the magnitude from 0.95 to 1.86 as p increases.
At p=16 the coefficient exceeds the baseline. This is symptomatic of extracting innovations
from a near-unit-root process — identifying variation shrinks, estimates drift.

**Stable across all specs:** δ→u at h=0–2, δ→v trough at h=4–6.
**Not stable:** peak horizon, peak magnitude, anything at h>8.

### Implication for estimation design

The original Block B design targets the full IRF path β(h), h=0,...,20 via Bayesian impulse
response matching. The instability means this path is not well-identified at long horizons.

**Three options identified (D10):**
- **(a)** Short-horizon targeting: match only h=0–2 (u) and h=4–6 (v). Conservative.
- **(b)** Baseline-to-baseline: run the same LP on model-simulated data. If the model's
  shock persistence produces comparable Wold contamination, the comparison is
  apples-to-apples and the full path is valid. **Must be tested (→ E9).**
- **(c)** Hybrid: match full path, reweight toward short horizons.

Option (b) is most favorable — it preserves the original estimation design — but requires
one model simulation to verify. This is now task E9.

## What was updated

| File | What changed |
|------|-------------|
| `Bartek analysis/part5_lp.py` | Added `attach_lagged_instruments()`, extended `run_lp_horizon`/`run_lp` with `lagged_instr_cols`, added section [11] with autocorrelation diagnostic + augmented LP (p=4,8,12,16) + comparison table + comparison plot. Wrapped section [8] plotting in try/except for pre-existing PIL bug |
| `context/findings.md` | Added "E8: Bartik instrument persistence test" section with full tables |
| `context/decisions.md` | D10 rewritten: now 🟡, four options (a–d), recommendation to test option (b) first |
| `context/pending_tasks.md` | E8 → ✅; added E9 (model-side LP test); M8 updated to gate on E9 |
| `context/README.md` | Reading order updated to this handout |
| Memory | `project-hump-shape-structural.md` updated with E8 results |

## Uncommitted changes

- `Data/Bartek analysis/part5_lp.py` — the augmented LP implementation (section [11])
- Context file updates listed above
- This handout
- Plot: `data/results/irf_delta_augmented_comparison.png`

## What to do next (priority order)

1. **E9: Model-side LP test.** Simulate a panel of 50 "states" from the model with
   calibrated shocks (ρ_δ=0.592, σ_δ=0.0669), heterogeneous industry shares, and run the
   baseline LP specification. If the model's LP also shows monotone buildup, the full IRF
   path is a valid Block B target. If it peaks at h=1–2, fall back to short-horizon targeting.
   **This gates the Block B estimation design.**

2. **Step 0 cascade** (still pending). Set `dest_ann = 0.0320` in `steady_state.jl`, re-run
   Comparisons A–D, regenerate mechanism PDFs, update §5.3 and §5.2.

3. **Settle D2** (PATH A vs B) and **D3** (β̂ ↔ β(θ) scaling).

4. **Section [8] plotting bug.** `overlay_plot_irf` crashes with `OSError: Invalid argument`
   on `lp_irf_delta_ld_qu_overlay.png` — likely a PIL/matplotlib version issue in econ2315
   conda env. Currently wrapped in try/except. Low priority.

## Key numbers

| Quantity | Value | Source |
|----------|-------|--------|
| Within-state ρ(B̃^δ) at lag 1 | 0.91 (median) | E8 section 11a |
| Within-state ρ(B̃^δ) at lag 8 | 0.79 (median) | E8 section 11a |
| δ→u stable features | h=0–2, β≈0.27–0.65 pp | E8, all specs |
| δ→v stable features | h=4–6, β≈−0.55 to −0.65 pp | E8, all specs |
| δ→u peak horizon (augmented) | h=10–13 (not stable) | E8 |
| δ→u peak magnitude (augmented) | 0.95–1.86 pp (not stable) | E8 |

## Open decisions

| # | Status | Summary |
|---|--------|---------|
| D1 | ✅ Settled | BED 0.0320, but not yet implemented in steady_state.jl |
| D2 | 🔴 Open | PATH A vs B — gates Θ_e contents |
| D3 | 🔴 Open | β̂ ↔ β(θ) scaling — gates Block B |
| D10 | 🟡 Partial | What IRF target for Block B? Test option (b) via E9 |
