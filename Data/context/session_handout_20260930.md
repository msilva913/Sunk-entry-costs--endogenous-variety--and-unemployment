# Session Handout — September 30, 2026
**Branch:** `costly_vacancy_reposting` (branched from `instruments_LP_coefficient`).

Resume point: **Chunk C** of the draft reposting-margin revision. Chunks A and B are written and
evaluated; two small fixes are pending in them (listed below) before/as you continue.

## Context in one paragraph

We are incorporating the **costly partial reposting margin** into `Data/Draft/Draft.tex`, chunk by
chunk. The user writes every revision; the assistant supplies structure, target equation forms,
and evaluation. Model derivation and the settled forms live in
[`model_equations.md`](model_equations.md) Block 2 (source of truth),
[`../Notes/delta_calibration_and_the_reposting_margin.md`](../Notes/delta_calibration_and_the_reposting_margin.md),
and [`../Notes/role_segmentation.md`](../Notes/role_segmentation.md). The approved plan is at
`~/.claude/plans/gentle-humming-snowflake.md`; its chunk structure is folded in below. Scope of
this draft pass: **model equations only** (recruiter block, vacancy LOM, resource constraint),
bottom-up. Deferred: Proposition 5, intro/mechanism prose, and the new parameter α in the
calibration section (see Deferred, bottom).

## Settled representation decisions (apply consistently in all chunks)

- **Notation.** Reposting rate `Λ_{r,t} = F(Q_t/α)`. Reactivation option value **`Q_t^{rep}`**
  `= Q_tΛ_{r,t} − α∫_0^{Q_t/α}χ dF` (draft `eq:Qrep`). The shortfall is written **inline** as
  `Q − Q^{rep}` (no symbol `D`). The mean cost integral is written inline (no symbol `M`). These
  renames avoid clashes: `Ψ` vs the shape index `ψ`; `D` vs the intermediary dividend `D^{int}`
  and duration multiplier `\mathcal{D}`; `M` vs the matching function `m`/SDF. Do not reintroduce
  `Ψ`, `D`, or `M`.
- **Avoid the nested expectation.** Keep the recruiter **surplus `𝒮 = J − Q` as a tracked
  variable**; write the JCC as `κ + K_t/q_t = E_t m_{t+1}(1−δ_t)F(χ_{t+1}^c)(J_{t+1}−Q_{t+1})` and
  let `eq:value_recruiter_surplus` carry the reposting term (a single `E_t`). Do **not** fully
  substitute the surplus into a closed-form JCC — that composes two one-period lookaheads into a
  `t+2`/nested term. This is representational, not economic; the model is first-order Markov if
  `𝒮` (or `J`) is a jump variable. Same move for the code (R8).
- **Threshold is `Q`** (the recovered vacancy value), not `Q−K`. **Nesting/limits:** the costless
  current model is recovered as `α→0`; the LOM *quantity* saturates at `Λ_r=1` (all repost, but
  costs still paid unless `α→0`).

## Draft state by chunk

**Chunk A — recruiter block (§sec:recruiters) — WRITTEN; one fix pending.**
- Done: `eq:value_recruiter_matched` separation branch → `s_t Q_{t+1}^{rep}`; new `eq:Qrep`
  defines `Q^{rep}`; costless-reposting sentence rewritten; `eq:value_recruiter_surplus` gains the
  `− s_t E_t m S (Q_{t+1}−Q_{t+1}^{rep})` term; notation clean.
- ⚠ **Pending:** `eq:jcc` (~L943) and `eq:jcc_wage` (~L985) are still the **old** fully-substituted
  forms, **missing** the reposting term. Fix per the "avoid nested expectation" decision: write the
  JCC referencing the surplus (`κ+K/q = E m S (J'−Q')`), and do the wage substitution **inside the
  surplus equation**, not in a re-expanded JCC. (Do not add the reposting term by substitution —
  that nests.)

**Chunk B — vacancy law of motion (§Environment, `eq:v_lom` ~L646) — WRITTEN; fixes pending.**
- Done: `Λ_r` inserted on the reposting inflow; `(1−δ_{e,t})` on the bracket; forward pointer to
  `sec:recruiters`.
- ⚠ **Timing bug:** `eq:v_lom` uses `Λ_{r,t-1}` but must be **`Λ_{r,t}`**. A position vacated at
  end of `t−1` is reactivated into a period-`t` vacancy worth `Q_t`, so the rate is
  `Λ_{r,t}=F(Q_t/α)`, pairing with the `t`-dated survival `(1−δ_{e,t})` and matching the recruiter
  Bellman. Your own prose (~L658) already says `Λ_{r,t}`/`Q_t`, so the equation contradicts the
  prose — fix the equation index to `t`.
- Minor: "survived both margins and the activation cost" → "…and chose to reactivate the vacancy";
  `\ref{sec:recruiters}` → `Section~\ref{sec:recruiters}`; optional sentence noting the complement
  `1−Λ_{r,t}` is permanently retired (the "missing middle").

**Chunk C — resource constraint / aggregation (§sec:aggregation, ~L1103) — NOT STARTED. Resume here.**

## Chunk C — structure to write

**C1. Add the aggregate reactivation cost `X_r`.** Reposted positions this period × per-position
reactivation cost. Reposted positions `= (1−δ_{e,t})·Λ_{r,t}·s(1−u)` (align timing with `eq:v_lom`
— note the `Λ_{r,t}` fix); per-position cost `= α·(mean χ among reposters)`.

⚠ **Convention reconciliation (the one open coefficient).** The draft writes the aggregate
continuation cost as `X_c = N_t·p_0·ψ_c·χ_t^c` (`eq:agg_fixed_costs`, ~L1109), which drops the
`(χ^c/χ_m)^ψ` factor that the exact integral `∫_0^{χ^c}χ dF = p_0 ψ_c χ^c (χ^c/χ_m)^ψ` carries.
**Write `X_r` using the same convention the draft uses for `X_c`**, so the two aggregate-cost
objects match. The exact integral form (from `model_equations.md` f[16]) is
`X_r = (1−δ_e)·s·(1−u)·ψ_c·Q·(Λ_r−1+p_0)` (interior regime `Λ_r<1`). Decide: adopt the draft's
simplified `X_c` convention for `X_r` too, or upgrade both to the exact integral. Note the choice.

**C2. Thread `X_r` into `eq:X_total` (~L1121), `eq:rc` (~L1185), `eq:gdp` (~L1134).** Either add
`X_r` as a third recruiting-cost addend (`X_t = X_v + X_f + X_r`) or carry it separately
(`Y_c = C + X + X_c + X_r`, `Y = Y_gross − X − X_c − X_r`). Recommend keeping `X_r` visible as its
own term.

**Unchanged — state explicitly (per `role_segmentation.md`, confirmed against `eq:firm_bellman`):**
the retailer block — firm value `ν_f`, exit cutoff `χ^c` (`eq:cutoff`), dividend `d_f`, free entry
(`eq:free_entry`), business-formation Euler (`eq:euler_bf`). The recruiter/retailer segmentation
quarantines the reposting option in the recruiter's `J`; these carry no reposting term and move
only through general equilibrium.

## Source of truth — settled forms (draft notation)

Conventions: SDF `m_{t+1}`; survival `(1−δ_t)F(χ_{t+1}^c)`; `Λ_t≡F(χ_t^c)` (`eq:Lambda`);
`F(χ)=(1−p_0)+p_0(χ/χ_m)^ψ` (`eq:F_chi`); `ψ_c≡ψ/(ψ+1)`.
```
Λ_{r,t} = F(Q_t/α) = (1−p_0)+p_0(Q_t/(α χ_m))^ψ          reposting rate
Q_t^{rep} = Q_tΛ_{r,t} − α∫_0^{Q_t/α}χ dF                reactivation option value (eq:Qrep)
shortfall = Q_t − Q_t^{rep} = Q_t(1−Λ_{r,t}) + α·M_t     (written inline; M_t = the integral)
α·M_t = ψ_c·Q_t·(Λ_{r,t}−1+p_0)                          interior closed form (Λ_r<1)
X_r = (1−δ_e)·s·(1−u)·α·M                                aggregate reactivation cost (see C1 convention)
```

## Immediate next steps (resume order)

1. Apply the two pending fixes: Chunk B `Λ_{r,t-1}→Λ_{r,t}` (+ minors); Chunk A `eq:jcc`/`eq:jcc_wage`
   restructuring to the surplus-reference (non-nested) form.
2. Write **Chunk C** (C1, C2 above); resolve the `X_r`/`X_c` convention.
3. Compile; run the nesting check (`Λ_r=1`, `X_r=0`, `α→0` recover the current equations).

## Deferred (out of scope for this draft pass; nothing lost)

- **Proposition 5** (`prop:ds_asymmetry`, ~L1423; `lem:vpre`, ~L4366; proof `app:proof_ds`,
  ~L4358): Part 1 becomes conditional on `Λ_r`. Separate re-derivation (task R3).
- **Intro / mechanism prose** asserting costless reposting (~L1480–1482 and the intro overview) —
  qualify once the equations are in.
- **New parameter α** in §sec:calib (~L2249): estimated block, `tab:calib_targets` (~L2489), priors.
  Entangled with **D2** (PATH A/B) and the **D1** cascade — recall `dest_ann=0.0320` only converges
  under PATH B (`dest_elast_target`), not the default PATH A. Also fix the superseded `7.54%` δ_e
  table row vs the updated 3.2% prose during that pass.

## Other open threads (unchanged from prior handout)

- **R8** implement the reposting margin in `steady_state.jl`/`run_solution_core.jl` (keep the
  recruiter surplus as a tracked jump — same non-nesting move). **R9** simulate for
  amplification/persistence; ties to **E9/D10** (model-side LP for the Block B IRF target).
- **D1 cascade** still owed (blocked on D2). **D2, D3** open.
