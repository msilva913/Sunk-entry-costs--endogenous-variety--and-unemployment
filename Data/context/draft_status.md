# Draft Status — `Draft/Draft.tex`
**Last updated:** October 13, 2026 (P8 theory-section cleanup: `ass:vacancy_fe`, `rem:bgm_nest`, N^max, captions, "vacancy free entry" terminology; numbering below is current). Oct 12: context cleanup; Oct 11–12: R15 option (b) spliced, theory propositions audited and corrected, X^c made exact, `ass:matching` tightened with `eq:Phi` in the main text)

**Build:** compiles with 0 errors and 0 missing graphics. Two undefined references remain,
`app:robustness` and `app:weighting_robustness` (unwritten appendices; P2); no duplicate labels.
Compile from inside `Data/Draft/` (the graphicspath leads with relative entries), **three
`pdflatex` passes** (mathtools `showonlyrefs` needs them), with `-synctex=1`. **Never run
`bibtex`** (it destroys the tracked `Draft.bbl`). Close TeXstudio first: an open file causes
EPERM on save. If an editor auto-compiles on save, the `.aux` can be truncated
mid-build; delete it and rebuild.

**Notation.** χ^c_t (not x^c_t); Λ_t ≡ F(χ^c_t); ς for the survival quantile. Reposting:
Λ_{r,t} = F(Q_t/α), option value Q^rep_t (`eq:Qrep`), shortfall written inline as Q − Q^rep;
the symbols Ψ, D, M are retired (clashes with ψ, D^int, 𝒟, m). 𝓡_t is the expected
discounted shortfall in the wage appendix. Label prefixes `eq:`/`tab:`/`fig:`/`sec:`/`app:`.

**Timing.** The draft is the reference model ([S10](decisions.md)); the code deviates in two
places until R8 ([`model_equations.md`](model_equations.md) "Code status"). Until R8 and the D1
cascade, every §5.3 number comes from the code's timing.

---

## Section status

**Oct 11 efficiency pass, Sections 1–3 (MS-directed).**
- Intro: the variety/segmentation/exit paragraph split into three; estimation described once;
  Fact 1 compressed; the 2001 aside shortened.
- Intro accuracy fixes: the LP summary now matches §4 (unemployment rises, +1.71 pp; vacancy
  rate −0.58 pp, not log points); the δ instrument is built from deaths, not closings;
  roadmap order fixed; `ref:conclusion` used.
- The commented scrap block before §2 was deleted.
- §2: reactivation stated once (MS's paragraph, with the timing anchor fixed); timing restated
  once; the dangling segmentation paragraph rewritten.
- §3: the conditional mean E(χ | χ ≤ χ^c) corrected to p_0ψ_cχ^c(χ^c/χ_m)^ψ/Λ; the X^r
  derivation compressed; the post-definition text de-duplicated; the wage paragraph tightened.
- `eq:exit_threshold_curve` did rely on a shortcut (ψ_c with no p_0, the p_0 = 1 model); fixed Oct 11 with the X^c correction (P5 in `pending_tasks.md`).
- MS kept the DMP paragraph (item 5).
- **Oct 12:** `ass:matching` tightened to (M1) f ≤ min{1, θ} and (M2) ε_f ∈ (0,1) nonincreasing;
  Φ (break-even recruiter compensation) defined in the main text as `eq:Phi`, with an
  interpretation of (C) (the elasticity of Φ in θ does not fall). Appendix J cites `eq:Phi`.


| Section | Status | Outstanding |
|---|---|---|
| 1 Introduction (incl. §1.1 evidence) | ✅ | "Self-correct" prose fixed Oct 11–12 (R11). LP summary matches §4. **Oct 5:** endogenous exit no longer claimed as a recession amplifier (contradicted by §5.3 — Flag 1), reworded to the δ_e decomposition plus the host distribution F for reactivation draws (`ee5a2a5`); new paragraph giving the measurement argument for why the business formation block is needed given GS (`bcee2cb`); Fact 1 gains the entry margin, the C7–C10 correlations and the one-for-one clause (`1e94613`). **Still owed:** the ξ sentence (gated on R14). The planned ξ/CK paragraph on the predetermined share is superseded by `prop:fe`'s user-cost contrast; if D10(a), the LP sentence switches to the h = 0–2 / h = 4–6 features (P7) |
| 2 Environment | ✅ | Convention A `eq:u_lom`; `eq:v_lom` with Λ_{r,t}. Timing prose (Oct 1, late) states the paper's convention: the aggregate state incl. s_t, δ_t is known at the start of t; only incidence is revealed at the end of t. `fig:Timing` unchanged (MS: consistent with the convention); it notes that s hits only incumbent matches |
| 3 Equilibrium | ✅ | Wage dependencies fixed Oct 1: 𝓡_t defined at `eq:repost_shortfall`; `eq:wage_eq` carries −s_t𝓡_t; `eq:surplus_wage` carries −(1−ϕ)s_t𝓡_t. Oct 5: prose notes s_t𝓡_t enters alongside K_t (an extra flow cost of a filled position); recruiter profits Π^int now net of X_t and X^r_t, which fixes the `eq:gdp` income identity (it was off by X even before reposting); `eq:phi_calib` gains −s𝓡. `eq:Lambda_r` defines the reactivation rate. `def:equilibrium` audited for completeness Oct 1: 4 equations for (N, C, χ^c, Y^c), 6 for (θ, K, Q, e, v, u) with `eq:theta_eq` added. The job creation condition `eq:jcc_eq` substitutes the surplus at t+1, so J_t and w_t drop out (`eq:surplus_eq` removed); 𝓡_t is carried with its own equation, so there is no nested expectation. Auxiliary list cites `eq:Lambda_r`, the matching rates, and 𝓡_t. Prop. 1 proof updated to match |
| 3.x Limiting cases | ✅ | **Oct 5:** post-Prop-3 passage rewritten (`2fe27da`). Four channels over AGS, not three — **profit dilution** added (Y^c/N in `eq:euler_bf`, distinct from variety, survives ρ ≡ 1). Second tier added: AGS does not price replacement, and AGS collapses measured and structural destruction (Λ = 1 ⇒ δ_e,t = δ_{t−1}). Ablation paragraph now names Comparison C only, glosses ω_δ at first use, and states that exit and reactivation must be switched separately rather than jointly through p_0. **Oct 11:** restructured (see Propositions). ⚠ **Commits §5.3 to three runs that do not exist**: novar-vs-AGS, α → 0, and a re-specified Comparison B (R8 clause v) |
| 3.x Shock transmission (`sec:mechanism`) | ✅ | **Rebuilt Oct 11** around `lem:theta` + `prop:fe` (free entry, δ–s theorem), `lem:accounting` (h = 1 accounting, any ξ, old Prop. 5 Parts 1 and 3), `rem:finite_xi` (numerical finite-ξ threshold, replaces the unproved "δ̄_e/τ̄ small" sufficient condition), `prop:exit_s` (old Part 2). Numbers: `gs_free_entry_beveridge.txt` (κ = 0), `gs_finite_xi_entry_offset.txt`; the 2.6%/0.7%/81% figures still from `run_prop5_*` (regenerate after R8) |
| 3.x Steady state | ✅ | Oct 5: `eq:e_theta` gains the entry that replaces non-reactivated positions (implicit in θ when α > 0); `eq:curve_jcc` gains −s𝓡(θ); `eq:destruction_curve` corrected for the atom (was valid only for p_0 = 1) |
| 4.1–4.2 Instruments, LP specification | ✅ | — |
| 4.3 Results | ⚠️ | δ–LD paragraph; joint-LP sentence; malformed `\ref` at the end (P1) |
| 5.1 Relation to Coles and Kelishomi | ✅ | δ̄ figures in the D1 cascade |
| 5.2 Calibration and estimation | ⚠️ contested | Documents PATH A, figures use PATH B ([D2](decisions.md)); β̂ ↔ β(θ) unstated ([D3](decisions.md)); α absent; `app:weighting_robustness` unwritten |
| 5.3 Inspecting the mechanism | ✅ | ⏳ Regeneration owed (R8 + D1 + PATH B). Numbers are checked against `mechanism_stats.txt`, which is the authority. The "Choice of shocks" paragraph now says the comparisons assume costless reposting |
| 5.4 Posterior estimates (`sec:posterior`, L3310) | ❌ empty stub | Blocked on the estimation build |
| 6 Conclusion | ❌ not started | Label `ref:conclusion`; all references now use it |

## Appendix status

**`app:loglin` — rewritten Oct 5, 2026 (`9582da1`, `10f3b60`).** ¶9 previously signed the entry response by reading one equation twice: `eq:K_decomp` and `eq:Q_ll_delta` are the same asset-pricing identity (φ_D = D·δ̄_e/(r+δ̄_e) identically), so it cannot sign K̂ and Q̂ independently. New ¶9a derives the dividend from the JCC (`eq:K_jcc_ll`): survival (−), surplus (−), and **congestion relief** (+), the last being the actual entry cushion. ¶9b records that at free entry Q is pinned so the identity *does* sign K (`eq:K_freeentry`, `eq:theta_freeentry`). ¶9c derives why δ̄_e/τ̄ decides the sign of v̂ (`eq:tau_loading`). ê^δ > 0 is now stated as a numerical result from Comparison D rather than derived. `eq:v_ll`'s destruction impulse corrected from a hat to a tilde (level deviation, matching its coefficient v̄^pre). **Owed (R12b):** a figure or `mechanism_stats.jl` number for ê^δ > 0, and the threshold T from the predetermined-vacancies note belongs in ¶7

| Appendix | Status |
|---|---|
| A Relation to literature | ✅ |
| B Data sources | ✅ |
| C Instrument diagnostics | ✅ |
| D Steady state | ✅ — reposting added Oct 1 (shortfall in the job creation condition and wage, Λ_r, Q^rep, X^r, entry `eq:e_ss` with the non-reactivated term); **X^c made exact Oct 11 (P5)**: X^c = Nψ_cωχ^c, ω = 1−(1−p_0)(1−δ)/(1−δ_e); N, L^e, L^c and a sign claim in the derivation subsections also corrected; `eq:profit_share_ss` split from the duplicate label |
| E Derivations — profit share, resource constraint, labor share, **wage determination** (rederived Oct 1: Convention A employment law with δ_{e,t+1}, the reposting shortfall 𝓡_t, the wage-substituted surplus with (1−ϕ)), log-linear propagation (`eq:u_ll` rederived for Convention A) | ✅ |
| F Additional comparisons (incl. Comparison D) | ✅ |
| G Limiting-case proofs (`prop:independence`, `prop:double_limit`, `prop:ags`, `prop:bgm_nest`) | ✅ **corrected Oct 11** (see Propositions) |
| H Proof of `lem:theta`, `prop:fe` (`app:proof_fe`) | ✅ new Oct 11, translated from the note |
| I Proofs of `lem:accounting`, `lem:vpre`, `prop:exit_s` (`app:proof_ds`) | ✅ restructured Oct 11; R12c argument removed |
| J Steady-state equilibria (`app:equilibria_proof`) | ✅ **rewritten Oct 11 (R16)**: log-convexity argument under `ass:matching` + (C) |
| `app:weighting_robustness`, `app:robustness` | ❌ referenced, do not exist |

---

## Propositions (audited Oct 11, 2026)

**Limiting cases restructured Oct 11 (MS):** the subsection now opens with the question (link to CK/GS), then `prop:independence` (1), clone replacement as the thought experiment (how to keep the destruction flow while abstracting from entry and exit) with `prop:ags` (2) and `rem:nesting`, then the double limit (3) as the justification that clone replacement, prima facie a heuristic, is a limit of the baseline. The failed single limit ε → ∞ is motivated by "business formation requires profits" (μ − 1 = ζ under DS-CES, back-reference to the intro's measurement argument); the δ → 0 detour is a footnote; the knife-edge caveat is replaced by the payoff (`lem:theta`/`prop:fe` are statements about this limit; convergence is of equilibrium conditions, not IRFs). Appendix G proofs reordered to match.

| Result | Status |
|---|---|
| `prop:independence` (1) | ✅ Statement now "for all values of the remaining parameters", which makes necessity exact (relaxing a condition gives dependence generically) |
| `prop:double_limit` (3) | ✅ **Claim (i) was wrong, corrected Oct 11.** Old `eq:N_lom_limit` put the limit Y^c → zL inside an O(ε) coefficient, dropping an O(1) term; old `eq:N_ss_limit` (δN = zL lnN/f̄_e) has 0 or 2 roots, not one, and its stability condition N̄ > e^{1−δ} was spurious. Now: N̄ = z̄L̄(1−δ̄)/(f̄_e(r+δ̄)) (= limit of `eq:curve_res_cons`) plus the period-by-period rule `eq:N_rule_limit` from the Euler equation; no autonomous N recursion, interior entry the only regularity condition. Normalization N̄ = 1 ⇔ f̄_e = z̄L̄(1−δ̄)/(r+δ̄), w.l.o.g. because allocations are independent of f̄_e. Check: `double_limit_check.jl`. A paragraph after the statement explains the one-way causation and the AGS link |
| `prop:ags` (2) | ✅ Well-posedness paragraph corrected (the "scalar equation for Q with root β(1−δ̄_e)" claim was false at finite ξ): BK verified numerically at σ = κ = 0 (`gs_finite_xi_entry_offset.jl` asserts 2 stable roots), analytic at ξ → ∞ via `lem:theta`. τ_t defined in the statement |
| `rem:nesting` | ✅ GS isomorphism (matching-function aside removed Oct 11, MS: belongs with probability bounds, not here); CK recovered exactly — CK set their flow vacancy cost c = 0 (p. 124), GS likewise (MS) |
| `rem:bgm_nest` (Remark 2; was `prop:bgm_nest`) | ✅ **Demoted to a remark Oct 13** (MS: it connects the model to BGM, but the subsection was overloaded); the subsection opener now announces it. §5.2's "nests CK as formalized in" now cites `prop:double_limit` + `rem:nesting` (it wrongly cited `prop:bgm_nest`). Proof rewritten: F is a point mass at 0 under p_0 = 0 (was "F(0) = 0"); bracket = N + N^e by N^e = zl^e/f_e and L^c = Y^c/(ρz) (the "gross-output identity" Y^c = zLρ/μ was false off steady state); ν^f = ρf_e/μ as in `eq:free_entry` |
| `rem:curves` (Remark 5; label renamed from `prop:curves` Oct 13) | ✅ Demoted; ε enters all three non-destruction curves through μ (not π_s); destruction curve depends on p_0, ψ, χ_m, δ̄; the identification sentence was removed |
| `prop:equilibria` (6) + `ass:matching` (Assumption 2) + SMCD footnote | ✅ **Oct 13:** the resource-constraint bound is now N^max (it clashed with Prop. 3's N̄), the statement says N_RC < N^max (the old "N_RC → N̄" was false), the footnote names the *low* equilibrium, and both figure captions are trimmed to objects (N6). **R16 done Oct 11.** The repair is conditional on (C), which is checked numerically. Stated under `ass:matching` (CRS; (M1) f ≤ min{1, θ}; (M2) ε_f ∈ (0,1) nonincreasing) and explicit condition (C), ln Φ convex in ln θ, with Φ defined in the main text (`eq:Phi`, Oct 12). Proof: ln N_JC − ln N_RC convex in ln θ ⇒ ≤ 2 zeros; positive at both ends ⇒ 0 or (generically) 2; double-limit uniqueness iff z > Φ(0). The old convexity-in-levels fix fails (Φ is never convex; N_JC not convex near 0). Under extrapolated Cobb–Douglas the low equilibrium needs q ≈ 14–21; under every bounded spec (clipped CD, HRW, 2P-SMCD) it is feasible with q ≈ 1, u ≈ 1/3. The SMCD discussion is a condensed footnote on the sentence introducing `ass:matching` (was `rem:smcd`, moved Oct 11, MS); it cites Silva (2026, SSRN 6956640) and DHRW (doi:10.1257/aer.90.3.482). Program: `equilibria_convexity_check.jl` (loads the SMCD repo's `MatchingSpec` library) |
| `ass:vacancy_fe` (Assumption 1), `lem:theta`, `prop:fe` (4) | ✅ **Oct 13:** new `ass:vacancy_fe` (V1 Q = x_m, e_t > 0; V2 0 < δ̄ ≤ τ̄; V3 (1−δ̄)f̄ + τ̄ < 1, so λ_u ∈ (0,1)) closes the gap that the proof used λ_u ≥ 0 without it in the hypotheses. Titles say "vacancy free entry"; the CK sentence cites p. 125 ("analogous to the free entry case") and p. 120. **Spliced Oct 11** (R15 (b)). **Discussion rewritten Oct 11 (MS):** the conditions are stated as |c_x|(ρ_x − ω_u) > ℓ_x (h = 1) and |c_x|[ρ_x(1−λ_u²)/(1+ρ_xλ_u) − ω_u] > ℓ_x (path), with no Γ_x in the main text (Γ survives only in the proof); the closed forms for c_δ, c_s, D_θ are in `app:proof_fe` (`eq:c_delta`); `tab:fe_irf` contrasts δ (ρ = 0.592) and s (ρ = 0.874) responses; the text leads with the two-channel contrast (tightness 0.79 vs 0.04, unemployment 0.036 vs 0.37 per percent of the shock). Notation translated: a, b → c_δ, c_s; D → D_θ; R_x → Γ_x; Φ → Φ_θ; λ, ω → λ_u, ω_u; ρ → ρ_x; ε_0 → x̃_0. Assumptions added: (1−δ̄)f̄ + τ̄ ≤ 1, interior entry e_t > 0. `cor:allh` is a sentence in the text, proved in Step 2 |
| `lem:accounting` (Lemma 2) | ✅ Old Prop. 5 Parts 1 (s) and 3 (δ) at any entry technology, p_0 = 0; the δ part now also states the unemployment (matching) condition |
| `rem:dmp` (Remark 3, new Oct 13) | ✅ Under `ass:vacancy_fe` with interior entry, the labor market block is isomorphic to a DMP economy with the same timing and margins in which every position bears the flow user cost K_t = x_m(r+δ_t)/(1+r): vacancy cost K_t and match output z_t − K_t. The wage and the job creation condition coincide term by term, v_t = θ_tu_t jumps, and v_pre matters only if e_t ≥ 0 binds. Two limits are stated after it. (1) Under δ shocks the counterpart's cost moves with δ_t, which is the user-cost channel; with δ constant it is DMP with relabeled productivity, so z and s shocks transmit as in DMP. (2) e_t ≥ 0 is the irreversibility of vacancy investment. No number is quoted for when it binds: the 9.3 % in `beveridge_free_entry.txt` is for x_v = 0.5, not the κ = 0 economy (N15). The flow-cost comparator (+0.49) is now described as this remark's counterpart with K fixed. The §2 sentence (~L902) no longer claims DMP is nested |
| `rem:finite_xi` (Remark 4) | ✅ **Reworded Oct 13:** 1 %/yr–τ̄ is the calibratable range; the jump to "every δ̄" lies between ξ = 10 and 100 (MS: no new number); the 85 % (GS economy, ξ = 1) and the next paragraph's 81 % (full model, p_0 = 0, 1/ξ = 1) are labeled. Replaces the old Part 3 sufficient condition (R12c). Numbers from `gs_finite_xi_entry_offset.jl`: δ̄*(ξ) = 3.8%/2.7%/1.6% per year at ξ = 1/2/10 (ρ_δ = 0.592); offset 0.85 at D1, ξ = 1; non-monotone in ξ |
| `prop:exit_s` (5) | ✅ Old Part 2. g(N, L^c) defined once (gross, from `eq:cutoff_eq`, Y^c = ρzL^c), entry labor l^e carried in the joint system; `eq:exit_cond` uses g_L > 0. ⏳ Numerics after R8 |
| `lem:vpre` (Lemma 3) | ✅ Unchanged |

**Limiting cases, Oct 13 (P8):** "three respects" (transmission channels; pricing of replacement; measured vs structural destruction). The comparison-design paragraph is a two-sentence pointer; its unique content (what the Comparison C–AGS gap measures; exit and reactivation switched separately; the α → 0 run) moved to §5.3's comparison list.

**Remaining known issues (not propositions).**
- **Theory gap: no result for partial reposting (R3).**
- **Found Oct 13, not fixed:**
  - §5.3's Comparison A sentence is garbled (`Draft.tex` ~L3142–3144: "The baseline sets ($\delta_e=\bar\delta$ … fixed.").
  - §5.3 describes Comparison B as p_0 = 0 with equal δ̄, while the old limiting-cases text said ω_δ → 0 at fixed steady-state destruction. Settle with R8 clause (v).
  - ✅ **Resolved Oct 13 by `rem:dmp`:** §2 (~L902) claimed infinitely elastic vacancy creation "nests the DMP baseline". It now points to the remark.
- A remark on selection and stability of the low steady state is owed (R18). `app:loglin` ¶5 still quotes δ̄_e/τ̄ ≈ 0.21 (D1 cascade); ¶5's "variance share" reading is loose. ¶9/9a survival re-dated to δ_{e,t+1} (S10); `eq:tau_ll` coefficient corrected to (1−s̄)δ̄_e/τ̄ with s_{t−1}; ¶9b/9c shortened to a pointer to `prop:fe`.

## Companion theory notes (uncited)

`Notes/Baseline_Blanchard_Kahn.md` and `Notes/AGS_Blanchard_Kahn.md`: BK holds under ε > 1 and
endogenous exit makes it easier to satisfy. No proposition needs ε > 2 any more (S9). Cite them
in the draft (P3). `Notes/beveridge_free_entry_GS.tex` is superseded by the draft (`prop:fe`).
