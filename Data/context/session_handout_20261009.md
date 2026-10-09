# Session Handout — October 9, 2026
**Branch:** `costly_vacancy_reposting` · **Last code/note commit:** `dea96eb` (Oct 7, "R15: add the
s-shock part; one free-entry proposition for both shocks")

> Written to be resumed in a fresh session, possibly on another machine, with no conversation
> memory. Everything needed is in this file or named by path.

> **Sync check.** As of Oct 9 the branch was level with `origin/costly_vacancy_reposting`
> (`03df797` and `dea96eb` pushed). This handout and the context updates are a separate commit.
> On the new machine, run `git pull` and confirm `git log -3` shows it.

---

## What this session was

Two strategy discussions, Oct 8–9, 2026. **No model runs, no draft edits, no code.** The outputs
are this handout, decision **D12**, tasks **R16** and **R17**, and annotations to R12, R15, D10 and
the proposition inventory.

1. What role is left for Prop. 5 (`prop:ds_asymmetry`) now that the R15 free-entry proposition
   (`prop:fe`) exists.
2. Should the paper pivot to theory-first and target the *Journal of Economic Theory*, moving the
   Bartik/LP work to a separate empirical paper? **Verdict: no-go for JET** (D12, open, MS).

---

## 1. Prop. 5 vs `prop:fe` — division of labor

`prop:fe` (in [`../Notes/beveridge_free_entry_GS.tex`](../Notes/beveridge_free_entry_GS.tex)) is
sharper wherever both apply. It gives iff conditions in parameters, covers the whole path rather
than only h = 1, and quantifies the δ–s asymmetry (R_δ = 0.045 vs R_s = 9.56 at D1, ρ = 0.592).
It also shows Prop. 5 Part 3's "δ̄_e/τ̄ small" fails at ξ → ∞ for an iid shock.

Prop. 5 still owns four things `prop:fe` cannot cover:
- **Finite ξ, the calibrated case.** The closed form exists only because ξ = ∞ (and ξ = 0) makes
  the system block recursive. At finite ξ, v_pre is a genuine state. There "δ̄_e/τ̄ small" is right,
  because it makes the entrant weight ē/v̄ = δ̄(v̄+1−ū)/v̄ small.
- **Endogenous exit** (Part 2, p_0 > 0), which nothing else covers.
- **Variety** (Part 1 needs only p_0 = 0, so σ, κ > 0 are allowed).
- **Partial reposting** as the second source of negative comovement in the draft's framing.

The two bracket ξ. At finite ξ a low δ̄/τ̄ delivers the negative slope; at ξ → ∞ only persistence
does, with the threshold in closed form.

**Recommendation (MS to decide, R15):** option (b).
- Lead with `prop:fe` as the δ–s theorem.
- Keep Prop. 5 Part 2 as a standalone exit result.
- Demote Parts 1 and 3 to a short general-ξ lemma or remark, with the R12c proof rewritten
  through ē/v̄ and `eq:K_jcc_ll`.

Option (a), keeping both whole, leaves two overlapping claims, one partly overturned by the other.

---

## 2. JET assessment (theory-first pivot)

### Verdict
No-go for JET as the paper stands. The theory reads as applied theory, and clearing JET's bar
would take months of uncertain new theory. A JET version also contradicts the stated objective
(JME / AEJ: Macro / JPE, `../CLAUDE.md`). Adopt the theory-forward *structure*, but aim at a
field macro journal.

### The theory inventory against a theory bar
- **Nesting and isomorphism, expository rather than new economics:** `prop:independence`,
  `prop:double_limit`, `prop:ags`, `prop:bgm_nest`. `prop:bgm_nest` is unaudited, and its proof's
  bracket is misstated (R12).
- **`prop:curves`** "follows directly from the expressions" — remark-level.
- **`prop:equilibria` — proof wrong as written (new task R16).**
  - Step 4 goes from "even number of crossings" to "generically exactly two". Even includes zero
    and four, and the proof concedes Δ > 0 everywhere is possible, so existence is not shown.
  - A tangency is codimension one, not two.
  - Step 5 takes ε → ∞ with f_e fixed, which is the degenerate limit (N → 0) the main text rules
    out before `prop:double_limit`.
  - Probable fix: N_JC = (μΦ/z)^{ε−1} is convex when ε > 2 and Φ is convex, while N_RC is concave,
    so there are at most two crossings. Add an explicit existence condition and redo Step 5 with
    f_e = f̄_e/ε. Verify the convexity of Φ(θ) first.
- **Prop. 5.**
  - Its conditions are written in endogenous objects (∂e/∂s, ∂e/∂δ).
  - Part 2 holds to first order in F′, with continuity in ρ_s.
  - The Part 3 proof has the R12c bug, and its sufficient condition fails at ξ → ∞.
- **`prop:fe` is the best theory in the project, but narrow.**
  - It is log-linear, with one shock at a time.
  - Variety, exit and reposting are off (σ = κ = 0, p_0 = 0).
  - It holds at the ξ → ∞ corner.
  - Its thresholds run through R = ℓ/|a|.

### Proximity to the closest papers (DOIs from `../Key papers/markdown/`)
- **Gabrovski–Silva**, JEDC 2025, DOI 10.1016/j.jedc.2025.105205. The δ/s split, the sunk vacancy
  cost and predetermined vacancies are GS's. GS has productivity shocks with constant δ and s, so
  R15 is an increment on the author's own paper.
- **Coles–Kelishomi**, AEJ: Macro 2018, DOI 10.1257/mac.20150040. Their abstract: vacancy creation
  "is less than infinitely elastic." R15 says that, with durable sunk vacancies and a separate δ
  margin, free entry suffices if δ is persistent. That is a real refinement, but CK's free-entry
  comparator is flow-cost DMP.
- **Bilbiie–Ghironi–Melitz**, JPE 2012, DOI 10.1086/665825. Variety enters the theory only
  through nesting and steady-state counting. **Every dynamic theorem switches the business
  formation block off.**
- What JET publishes in search and macro theory was **not checked against specific papers**:
  citation unverified; source missing.

### Weakest links, in order
1. The `prop:equilibria` proof.
2. The headline theorem is a linearized special case with variety off.
3. Prop. 5's conditions and the R12c bug.
4. Numerics standing in for proofs (ê^δ > 0, R12b; the Part 2 numerics).

### What would raise the odds
| Addition | Effort | Risk |
|---|---|---|
| R16: repair `prop:equilibria` | 1–2 weeks | Low |
| R12c + R12 audits | 2–4 weeks | Low–medium |
| R17: nonlinear `prop:fe`. At ξ → ∞ the JCC has no states, so with Markov (δ, s) the exact policy θ_t = Θ(δ_t, s_t) should exist globally. Prove existence and uniqueness (contraction), then monotonicity, then the u and v signs | 4–8 weeks | Medium |
| Finite-ξ results: corners plus continuity or monotonicity in ξ, not a closed form | 2–4 months | High |
| Variety inside the dynamic theorem (ζ > 0 breaks recursivity) | 2–3 months | High |
| Efficiency / Hosios with sunk vacancies and variety | 1–2 months | Medium |
| Drop the quantitative sections | ~1 week | Abandons the objective |

Only the two high-risk rows would change the paper's category.

### The quantitative part
- **Correction to the pivot's premise.** The SVAR is bivariate (z, δ), Cholesky with z first; s
  comes from a univariate AR(1) on s = (τ−δ_e)/(1−δ_e) (part6b; ρ_s = 0.874, σ_s = 0.0854).
- **At JET** the section would be at most a one-figure illustration, and full-information
  estimation should go.
- **At a field journal** keep the section. Calibration plus moment matching is enough;
  full-information estimation is optional and deferrable.
- **Live problems either way:**
  - model cor(u,v) is +0.995 against −0.804 in the data until D11 and R8;
  - there is no δ→u hump;
  - all §5.3 numbers await regeneration.

### Venues if not JET (ranked by fit)
1. **RED.** Moderate rewrite.
2. **JME.** The stated objective; needs R8 and D11.
3. **AEJ: Macro.** CK's venue; needs some empirical motivation kept.
4. **JEDC.** GS's venue; safest, below the objective.

TE is not a fallback (same bar as JET); IER is generic. Stripping the paper for JET then going to
JME or AEJ: Macro after a rejection would cost months of rebuilding.

### Split-paper risk
- **The CK departure narrows** to GS's δ/s split, the user-cost channel and BED measurement
  (δ_e/τ = 0.087). The severity placebo, the "headline validation" (principle 15), leaves with
  the Bartik work.
- **The empirical paper is not viable alone:**
  - QU placebo fails (corr −0.971);
  - LD is not separately identified (r = 0.334);
  - the δ→u peak is unstable (E8, instrument ρ = 0.91).
- **What survives:** δ→u at h = 0–2, the δ→v trough at h = 4–6, and the severity placebo. Keep
  these in this paper as motivating evidence (D10 option (a)) instead of splitting.

---

## 3. Recommended path and next tasks

Theory-forward applied-theory paper for JME, AEJ: Macro or RED:
- lead with `prop:fe` (R15 option (b));
- keep the SVAR-calibrated mechanism;
- cut the Bartik section to robust short-horizon evidence;
- defer full-information estimation.

1. **R16** — repair the `prop:equilibria` proof.
2. **R15** — MS decides (a)/(b); then splice. **R17** — test the nonlinear `prop:fe`.
3. **R12c, R12** — rewrite Prop. 5 Part 3's proof; recast Prop. 5 as the general-ξ statement plus
   Part 2.
4. **R8 + D1 + PATH B** regeneration and the **D11** go/no-go — the quantitative Beveridge sign must
   be fixed for any venue.
5. **D10 → option (a)** under D12's recommended structure; D3 and D4 then largely become moot.

## Open decisions for MS
- **D12** (venue / theory-first pivot) — new, recommendation above.
- **R15 (a) vs (b)** — recommendation (b).
- **D10** — recommendation (a) if IRF matching leaves the paper.
- Standing constraint: **do not edit `Draft.tex`** until these are settled.

## Read first on resumption
1. This file.
2. [`../CLAUDE.md`](../CLAUDE.md).
3. [`decisions.md`](decisions.md) D12, D10, D11.
4. [`pending_tasks.md`](pending_tasks.md) R15, R16, R17, R12c.
5. [`../Notes/beveridge_free_entry_GS.tex`](../Notes/beveridge_free_entry_GS.tex) (PDF alongside).
6. [`draft_status.md`](draft_status.md) — the proposition inventory.
