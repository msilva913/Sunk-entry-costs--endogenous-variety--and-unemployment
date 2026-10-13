# Session Handout — September 30, 2026 (session B)
**Branch:** `costly_vacancy_reposting` · **Last commit:** `ac9920a`
**Predecessor handout:** [`session_handout_20260930.md`](session_handout_20260930.md) (session A,
which set up the chunked draft revision)

> **Written to be resumed on a different machine with no session memory.** Everything needed
> is either in this file or named by path. Nothing relies on conversation context.

---

## 1. Where the draft pass stands

The chunked revision that session A planned is **finished**. All three chunks are in
`Draft/Draft.tex` and the document compiles clean: **0 errors, 0 missing graphics, 4
pre-existing undefined references** (`sec:conclusion`, `app:robustness`,
`app:weighting_robustness`, `eq:labor_C_N`) and one pre-existing multiply-defined label
(`eq:profit_share`). Those five predate the reposting work — do not treat them as regressions.

| Chunk | Content | State |
|---|---|---|
| A | Recruiter block: `eq:Qrep`, `eq:value_recruiter_surplus`, `eq:jcc` restructured to reference the surplus, new `eq:surplus_wage` | ✅ done (MS wrote it) |
| B | `eq:v_lom` with `Λ_{r,t}` on the reposting inflow | ✅ done, timing corrected |
| C | `eq:agg_repost_costs` (with intermediate algebra), threaded into `eq:rc` and `eq:gdp` | ✅ done |
| — | Equilibrium definition `def:equilibrium` updated | ✅ done |
| — | Propositions 1 and 3 audited | ✅ done |

`X_t` was deliberately **not** folded into: `eq:X_total` remains `X_v + X_f` (vacancy
*creation* costs, from the entry distribution G), and `X_t^r` is carried as its own term
because it is reactivation drawn from F. `eq:rc` is `Y_t^c = C_t + X_t + X_t^c + X_t^r`
(line ~1269); `eq:gdp` subtracts all three.

---

## 2. Resume here

**Priority 1 — R3, and it is more urgent than session A assumed.** `app:proof_ds` contains a
step that the reposting margin **breaks**, not merely qualifies. The proof argues that at
ρ_s = 0 the right-hand side of the JCC at t+1 "is invariant to `s_t` at leading order." That
was true when the surplus had no reposting term. `eq:surplus_wage` now carries
`−s_t E_t m(1−δ)Λ(Q_{t+1}−Q_{t+1}^{rep})`, which depends on `s_t` directly. The dangling
citation there was repointed so the draft compiles, but **the argument itself is not fixed and
should not be treated as sound.** Find it by searching `Draft.tex` for
`invariant to $s_t$`. `lem:vpre` and the Part 2 Jacobian touch the same channel.

**Priority 2 — two conventions that must be settled before R8 (code implementation).** Both
are recorded inline in [`model_equations.md`](model_equations.md); neither is decidable from
the equations alone.

1. **X_r timing.** As written the draft pays reactivation costs on a *lagged* flow valued at
   the current threshold: `X_{r,t} = (1−δ_{e,t})·s_{t−1}(1−u_{t−1})·α∫_0^{Q_t/α}χ dF`. The
   draft's LOM already carries lagged flows so this is natural there, but the code's f[22]
   uses unprimed `s·sbar·(1−u)`, so a t-dated `X_r` needs a lag the code does not track.
   Alternative: **pay at separation**, `X_r = (1−δ_e)·s·sbar·(1−u)·α·M'`, no new state. That is
   a different model, not a notational variant. See f[16].
2. **δ_e dating.** The draft writes `v_t = (1−δ_{e,t})[…]`, dating survival at the LHS period.
   The code's f[22] writes `(1−δ_e)` unprimed with `v_pret'` on the left, dating it at the RHS
   period. Under the draft's convention `Λ_{r,t}` and `(1−δ_{e,t})` pair; under the code's they
   do not. One is wrong or they differ deliberately. See f[22] and
   [`../Notes/LOM_timing_consistency.md`](../Notes/LOM_timing_consistency.md).

**Priority 3 — deferred from session A, still deferred.** α in `sec:calib`: it is **not yet in
the calibration section at all** (verified). Needs the estimated block, `tab:calib_targets`,
and priors. Entangled with **D2** (PATH A/B) and the **D1** cascade; recall `dest_ann = 0.0320`
converges only under PATH B. Also fix the superseded 7.54% δ_e table row against the updated
3.2% prose during that pass.

**Also outstanding:** `Λ_r` has **no labelled defining equation** — it appears inline as
`Λ_{r,t}=F(Q_t/α)` in the equilibrium definition's auxiliary list.
[`model_equations.md`](model_equations.md) f[new] recommends giving it its own equation as a
tracked control, which would make its IRF reportable for R9 and let the definition cite a label.

---

## 3. Settled this session — do not re-litigate

### Timing: exit and reactivation share Stage 2

The draft's own `fig:Timing` puts **endogenous exit at Stage 2**, right after the aggregate
state is realized, not at Stage 5. `model_equations.md`'s stage table said Stage 5 and was
corrected. The reactivation decision sits at that same node, which is also where `Q_t` is
determined (`e_t = G(Q_t)`), so it conditions on **realized** `Q_t` rather than `E_{t−1}Q_t`.
That is what makes `Λ_{r,t}` procyclical and position destruction countercyclical, and it
confirms `Λ_{r,t}` (not `Λ_{r,t−1}`) and its pairing with `(1−δ_{e,t})`.

**Two agents, two draws, two thresholds.** Exit is the **retailer's** decision on its product
line against `χ_t^c`. Reactivation is the **recruiter's** decision on a vacated position
against `Q_t`, because recruiters hold the vacancies. Same stage, same distribution F, but
distinct draws. Retailer exit still gates reactivation — a position serving a withdrawn product
line dies regardless — which is the `(1−δ_{e,t})` factor.

### α·M is per *vacated* position

`M = ∫_0^{Q/α} χ dF` is a partial (unconditional) expectation, so `α·M = α·Λ_r·E[χ|χ≤Q/α]`
— **it already contains Λ_r, which must not be applied twice.** Closed form and its reading:

```
α·M = ψ_c·Q·(Λ_r − 1 + p_0)
      └────┘  └──────────────┘
    mean cost   probability a vacated position
    among       is a paying reposter
    payers
```

α cancels because it scales cost and threshold in opposite directions. Three distinct objects:
`α·M` (per vacated position — what `X_r` needs), `α·M/Λ_r` (per reactivated position), `ψ_c·Q`
(per position that actually pays). Derivation is in `model_equations.md` §"Deriving M".

### X_r uses the exact integral — no shortcut

The draft's `X^c` previously dropped the `(χ^c/χ_m)^ψ` factor. Applying that shortcut to `X_r`
would be wrong: the relative error is `p_0/(Λ−1+p_0)`, which is **0.27%** for `X^c` (Λ ≈ 0.9986)
but **2.5× at Λ_r = 0.7**, and **qualitatively wrong at the floor** `Λ_r = 1−p_0`, where exact
`X_r = 0` because only the free atom reposts. Good for `X^c` *because* survival is near one;
bad for `X_r` *because* reactivation is not. Both now use exact forms.

### Propositions 1 and 3

**Prop 1 (`prop:independence`) holds with its three conditions unchanged.** `p_0 = 0` collapses
F to a point mass, so `Λ_r = 1`, `Q^rep = Q`, `X_r = 0` — reposting vanishes *with* exit. Had
reactivation drawn from a separate distribution, a fourth condition would have been needed; this
is a second payoff of reusing F. Necessity unaffected: `Λ_r` and `Q^rep` depend only on `Q_t`
(already a labor-block variable) and `X_r` flows labor→formation like `X_t`. **`J_t` was added
to the stated labor block**, matched by `eq:surplus_eq`.

**Prop 2 (`prop:double_limit`) unaffected** by the same argument (it assumes `p_0 = 0`).

**Prop 3 (`prop:ags`) holds, and gained the clone-replacement motivation** it never had — the
motivation was not in any earlier draft version (checked `git log -S`); only
`Notes/AGS_Blanchard_Kahn.md` named the device. Its proof asserted `p_0 = 0` "under clone
replacement", but no-exit follows from **either** `p_0 = 0` **or** `f_e ≥ χ_m` (since
`χ_t^c = f_e` there), and only the first kills reactivation. Now stated explicitly.

### Channel switches — AGS bundles three, so switch one at a time

| To remove | Set | What stays on |
|---|---|---|
| Variety effects | `ζ → 0` (`ε → ∞`) | exit, reactivation |
| Endogenous exit | `ω_δ → 0` | **reactivation** (p_0 untouched; `χ^c = χ_m` gives Λ = 1 with Λ_r < 1) |
| Reactivation | **`α → 0`** | **endogenous exit** (`Q/α → ∞` so Λ_r → 1, α·M → 0; Λ = F(χ^c) untouched) |
| Both exit and reactivation | `p_0 = 0` | — this is the AGS / clone-replacement route |

So **do not read channel contributions off baseline-vs-AGS**: it moves all three. Recorded in
the draft at the AGS paragraph and in R9, which implies a **fourth mechanism comparison**
alongside A–D (an `α → 0` run).

---

## 4. Flagged but not fixed

**Block 2's per-position-draw justification may be unnecessary.** It argues that one draw per
*firm* would accumulate size heterogeneity and break the DS-CES aggregation `ρ = N^{1/(ε−1)}`.
That attributes positions to retailers. Under the segmentation positions belong to
**recruiters**, and retailers rent labor services competitively, so the distribution of vacated
positions across recruiters never reaches retailers — only the aggregate vacancy stock enters
matching, and the LLN delivers `Λ_r` either way. **If so the symmetry concern dissolves** and
per-position draws are simply the natural reading rather than a requirement. The draft's
§Environment now states the aggregation point in the segmentation form. Confirm before relying
on the argument in print.

---

## 5. Environment notes for another machine

- **Python:** no interpreter on PATH under the plain name. Use
  `C:/Users/m.silva/AppData/Local/anaconda3/python.exe`. `pandas`, `numpy`, `statsmodels`
  present; **no** `pypdf`/`PyPDF2`/`fitz`/`pdfplumber`.
- **`pdftotext`** is on PATH (ships with Git for Windows, `mingw64/bin`). Used by
  `Key papers/pdf_to_markdown.py`, which mirrors the PDF library to greppable text with page
  markers under `Key papers/markdown/` (gitignored — regenerate with
  `python pdf_to_markdown.py`). See [`pipeline.md`](pipeline.md) Part 3.
- **LaTeX:** `pdflatex` on PATH. Compile **from inside `Data/Draft/`** — the graphicspath now
  leads with relative entries, which resolve against the working directory. Two passes for
  references; pass `-synctex=1` or `Draft.synctex.gz` gets deleted.
- **Graphics path was fixed this session.** It had accumulated three usernames and the entry
  for the Bartik results folder was commented out while the active one pointed elsewhere,
  breaking four figures. Relative paths now come first; **do not delete them when adding a
  machine.**
- **Heredoc caveat:** the Bash tool collapses `\\` to `\` inside heredocs, which mangles LaTeX
  and regex. For multi-line text with backslashes, write the file with the Write tool and splice
  it in, or build strings with `chr(92)`.
- **File locks:** a viewer holding `Draft.pdf`/`Draft.aux` produces
  `! I can't write on file 'Draft.aux'`. Retry, or close the viewer.

---

## 6. Open decisions (unchanged status unless noted)

| # | Status | Summary |
|---|--------|---------|
| D1 | ✅ Settled | `dest_ann = 0.0320`, δ_e/τ = 0.087, **still not implemented** in `steady_state.jl:613`. Read as a *lower bound*; literature brackets δ at [3.2%, 10%]/yr |
| D2 | 🔴 Open | PATH A vs B — gates Θ_e and the α calibration pass |
| D3 | 🔴 Open | β̂ ↔ β(θ) scaling — gates Block B; overlaps E9 |
| D10 | 🟡 Partial | Block B IRF target; needs E9 (model-side LP) |
| D11 | 🟡 Formulation settled | Reposting margin: formulation settled Sept 28, now written into the draft. The go/no-go for the paper is de facto made on this branch |

---

## 7. Task pointers

Full list in [`pending_tasks.md`](pending_tasks.md). Reposting-specific: **R3** (Prop 5 /
`lem:vpre` re-derivation — start here), **R8** (code implementation, blocked on the two
conventions in §2), **R9** (simulate; now carries the channel-switch table and the fourth
comparison), **R2** (verify α, σ_s, δ_e separately identified — α and δ_e act through the same
channel so cor(u,v) cannot separate them; the LD→v IRF should), **R5** (add the GS /
Shao-Silos / CK-consistency discussion to §5.2), **R6** (whether δ_e moves into Θ_e).

Also owed and independent of reposting: the **D1 cascade** (blocked on D2), **part9** still uses
permanence-adjusted closings for Fact 2, **part10** needs only a re-run.
