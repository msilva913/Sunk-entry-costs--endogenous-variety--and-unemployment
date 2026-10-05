# Role of each model ingredient
**Written:** October 5, 2026 · **Revised:** October 5, 2026 after checking CK and GS directly
(§2) and the markup–variety identity in `Draft.tex` L560.

Audit of whether each ingredient is indispensable to the paper's objectives. Objectives:
**O1** composition of job loss (δ vs s) drives persistence, because destruction removes vacancy
slots; **O2** clear departure from CK/GS (δ vs s, endogenous entry/exit, variety as a third
state); **O3** fit labor moments (cor(u,v) ≈ −0.80, σ(θ)/σ(LP) ≈ 11.7) via SMM + Bartik IRF
matching.

## 1. The ingredients

| Ingredient | Role | What relaxing it does | Verdict |
|---|---|---|---|
| Markup (μ > 1) | Profits sustain finitely many product lines under sunk entry; motivates the recruiter/retailer split (avoids Stole–Zwiebel) | μ = 1 leaves zero profits, so sunk entry sustains no finite N and a δ shock has no product line to destroy | **Prerequisite, not an ablatable ingredient.** See the identity warning below |
| Love of variety (ζ > 0) | N→ρ→w^int→JCC feedback | Comparison C (forces ρ ≡ 1 via `calibrate_shares_no_variety`, so it does isolate ζ from μ): w^int response to z ≈ ¼ smaller within a year, equal on impact | **More than "moderate"** — it is what makes the δ-shock Beveridge sign survive free entry (§2). Flag 3 |
| Finitely elastic vacancy creation (ξ < ∞) | CK: predetermined vacancies, slots rebuilt at sunk cost. GS: congestion in K delays entry, giving a **hump-shaped vacancy response** (they cite Fujita–Ramey 2007) | Comparison D (**δ shock only**): sign robust even under free entry; ξ sweep moves the δ→u peak only h = 1→2 | **Indispensable — for the s margin, the hump, and the tightness–unemployment relation; not for the δ sign.** Flag 2, §2 |
| Product-line entry (f_e, Euler) | N dynamics: firm stock ≈ 80× more persistent than shocks | No entry / clone replacement (AGS) kills persistent variety loss | Indispensable |
| Separate vacancy and firm entry | u, v, N independent states (vs Schaal one-for-one) | Merging ties N and v one-for-one | Indispensable for O2 |
| Endogenous exit (p_0, ψ, ω_δ) | Intro sells it as recession amplification | Cutoff barely moves; Comparison B: endogenous-exit spec recovers *faster* | **Flag 1** (below) |
| Partial reposting (α) | Beveridge curve: reposting inflow is an identity | Without it cor(u,v) = +0.995 vs −0.80; only alternative (shrink s) contradicts data and CK | Indispensable for O3; **flag 4**: not yet coded/verified |

⚠️ **Markup and variety are one parameter in the baseline.** Under DS-CES,
ζ(N) = μ(N) − 1 = 1/(ε−1) (`Draft.tex` L560): the taste for variety *is* the net markup. Moving
ε cannot ablate one without the other, and ε → ∞ kills both. Only the general CES specification
disentangles them (L564). Comparison C is valid because it forces ρ ≡ 1 by hand rather than
moving ε. **There is no corresponding experiment for the markup**, and there cannot be a useful
one: μ = 1 removes the profits that sustain a finite N, so the object a δ shock destroys ceases
to exist. The markup belongs in a different category from the rest of the table.

## 2. Why ξ matters: the δ and s margins are not symmetric

> ⚠️ **This section's central claim is wrong (found October 5, 2026). Do not build on it.** It asserts that variety is what makes the δ-shock Beveridge sign survive free entry, citing Comparison D. `app:comparison_D` attributes that result to the **low calibrated δ̄_e/τ̄**, not to ρ(N) — and Comparison D holds variety *on* in both arms (it varies only ξ_inv), so it cannot isolate variety at all. `prop:ds_asymmetry` Part 3's sufficient condition is likewise δ̄_e/τ̄ small, with no variety parameter in it. The threshold derivation in [`predetermined_vacancies_and_the_CK_special_case.md`](predetermined_vacancies_and_the_CK_special_case.md) confirms this: T depends on δ̄_e/τ̄ alone, so an **augmented GS model inherits the result**. What is true is that variety *reinforces* the tightness decline through w^int = ρ(N)z/μ; its separate contribution to the sign has never been isolated, which needs Comparison C crossed with Comparison D (ρ ≡ 1 at ξ_inv → 0) and is not run. **Flag 3 and the ξ/variety substitutes conclusion in §5 both rest on this and are suspect.** Rewrite before R14 is designed.

**CK have one separation margin with δ_e = τ**, so every separation destroys the position. In our
taxonomy **CK's separation shock maps to our δ, not our s.** CK's result is that with a less than
infinitely elastic creation process, such a shock makes vacancies fall as unemployment rises, so
no ad hoc zero-separation-shock restriction is needed for Beveridge correlations. They also show
that free entry implies tightness is orthogonal to unemployment conditional on productivity,
which the data reject.

CK need ξ < ∞ because they have **no variety channel**: destroying positions does not lower the
value of creating new ones, so under free entry, entry replaces them. We have one for δ, which is
why Comparison D finds our δ sign robust under free entry and CK's is not.

**But the s margin has no variety channel.** Non-reposted positions are retired without
withdrawing a product line — only δ moves N (`model_equations.md` f[24]). So the s margin sits
exactly where CK's shock sits.

| Margin | What stops entry undoing the destruction |
|---|---|
| δ | **Variety** (N↓ → ρ↓ → Q↓). Works even under free entry — Comparison D |
| s, non-reposted | **Only ξ < ∞.** No variety channel exists here |

Two consequences. ξ and variety are **substitutes** on the δ margin, which is why Comparison D
could retire ξ without changing the sign. And ξ and α are **complements** on the s margin: α
creates the destruction, ξ stops entry from undoing it. Under free entry Q is pinned
(Comparison D: ξ_inv → 0 gives Q ≈ x_m) and entry expands to replace whatever is retired, so
ablating ξ while holding α would destroy the Beveridge gain α exists to deliver — and would
read, wrongly, as "α does nothing."

**Confidence.** The δ row is verified (Comparison D). The s row is inference from CK plus the
structural asymmetry: nothing has varied ξ under an s shock, and α is not coded. One extra cell
in R14 settles it.

## 3. Flag 1 — would estimating `dest_elast_target` make endogenous exit important?

**Probably not by itself.** Evidence from `Programs baseline/run_prop5_parts23_check.jl`
(p_0 = 0.5 PATH B baseline; indicative until R8 fixes the code's exit dating):

1. **Broer et al. (IER 2025) do not argue for a higher elasticity.** Their calibrated
   "separation elasticity" is **1.0**, the elasticity of *job separations* with respect to the
   value of a job (Table A.3), targeted so separations explain 40% of unemployment variance.
   That is the comparable object to our `dest_elast_target` (elasticity of the endogenous exit
   rate to the cutoff, ψς/(1−ς)), **not** our Pareto shape ψ. Ours is already 5 (placeholder),
   five times theirs. ⚠️ Earlier notes ("Broer et al. prefer an implied ψ near 1", decisions D2)
   compared the wrong objects; ψ = 1 here would mean an exit elasticity of several hundred.
   Also, Broer's margin is *all* job separations, which our model treats as exogenous (s).
2. **Small base.** Endogenous exit is 0.0014/month (D1) to 0.0033 (0.0754), vs τ = 0.031.
3. **Anchored cutoff.** χ^c = per-line profit + ν^f; the option value ν^f = f_eρ(N)/μ is
   97–99% of χ^c and moves only with slow-moving N. Peak |χ̂^c| after a 1-SD shock: 7.7e-4 (z),
   6.9e-5 (δ), 2.9e-4 (s) at D1. Peak endogenous-exit response to z: 5.2e-6/month, about 6% of
   a 1-SD exogenous δ impulse (9.1e-5); 9% at dest_ann 0.0754.

**The one data moment that could change this:** the part6b VAR coefficient of exit on lagged
productivity, â_21 = −1.05 (HC3 SE 0.83), i.e. exit elasticity to z ≈ 1. The model's
elasticity of δ_e to z is ≈ 0.2–0.3 (peak |δ̂_e| per 1% z; ⚠️ the program line currently
divides by 100σ_z and prints 0.002–0.003 — fix to divide by σ_z, the IRFs there are plain log
deviations). If estimation targets that moment, `dest_elast_target` would rise roughly 3–5×
(to ~15–25), making endogenous exit a moderate amplifier of z shocks (on the order of a
quarter of a δ impulse) while staying second-order for δ and s shocks. The moment is
imprecise. Linear scaling is a rough guide only.

**Implication.** Either add the exit-on-productivity moment to Block M and let the data speak,
or present endogenous exit as (i) the steady-state split of δ_e and (ii) the host distribution
F for reactivation costs (one parameter, Prop 1 needs no extra condition), not as a recession
amplifier. The intro's second contribution needs rewording in the latter case.

## 4. Other flags

2. **Finitely elastic vacancy creation:** the Comparison D evidence is **δ-only** and does not
   speak to CK's claim, which concerns a position-destroying shock without a variety channel —
   our s margin. Test ξ under an s shock, and the α × ξ interaction, before calling it dispensable.
   Separately, quantify how much δ→v persistence and hump shape it buys.
3. **Variety effects** are the main departure from GS and are load-bearing: they substitute for
   ξ on the δ margin (§2). Report their contribution to Block M moments and to the δ→v sign
   under free entry, not only to w^int. Also resolve the BK inconsistency (notes: ε > 1
   suffices; ξ×ε sweep: failures for ε < 3). ε > 2 is no longer needed by Prop 5.
4. **Partial reposting:** indispensable but unverified (R8, R9), identification vs δ_e (R2);
   turns the δ–s asymmetry into a quantitative result.
5. **Active core may be narrow:** δ→u peak 0.012 pp vs 1.71 data (partly D3 scale), half the
   amplification target, no hump. Documented bite: slot destruction + N persistence, partial
   reposting (pending), variety (moderate on w^int, larger on the δ sign). Endogenous exit and
   κ untested.

## 5. Recommendation

After R8: one-at-a-time ingredient ablation on one calibration — ζ→0, ω_δ→0, α→0, ξ→∞, κ = 0,
σ = 0 — reporting Block M moments and δ IRFs. Task R14 in `../context/pending_tasks.md`.

Two cells R14 must add, from §2:

- **ξ→∞ under an s shock** (not only δ), which is the experiment CK's result actually concerns.
- **ξ→∞ jointly with α > 0**, to test whether free entry undoes the reposting margin. If it
  does, the ingredients are not separable and the one-at-a-time design needs this pair reported
  together.

The markup is not in the ablation list, by §1: it is a prerequisite, and μ = 1 removes the
object under study rather than switching off a channel.
