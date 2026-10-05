# Role of each model ingredient — October 5, 2026

Audit of whether each ingredient is indispensable to the paper's objectives. Objectives:
**O1** composition of job loss (δ vs s) drives persistence, because destruction removes vacancy
slots; **O2** clear departure from CK/GS (δ vs s, endogenous entry/exit, variety as a third
state); **O3** fit labor moments (cor(u,v) ≈ −0.80, σ(θ)/σ(LP) ≈ 11.7) via SMM + Bartik IRF
matching.

| Ingredient | Role | What relaxing it does | Verdict |
|---|---|---|---|
| Markup (μ > 1) | Profits sustain finitely many product lines under sunk entry; motivates the recruiter/retailer split (avoids Stole–Zwiebel) | ε→∞ with f_e fixed sends N→0 (needs `prop:double_limit` scaling) | Indispensable for O2; **not** for the core δ–s mechanism (GS has it without market power) |
| Love of variety (ζ > 0) | N→ρ→w^int→JCC feedback | Comparison C: w^int response to z ≈ ¼ smaller within a year; equal on impact | Needed for O2, **quantitatively moderate** (flag 3) |
| Convex vacancy cost (ξ < ∞) | CK: predetermined vacancies, slots rebuilt at sunk cost | Comparison D: δ-shock negative comovement robust even under free entry; lower convexity raises z amplification; ξ sweep moves the δ→u peak only h=1→2 | **Flag 2**: not needed for the Beveridge signs |
| Product-line entry (f_e, Euler) | N dynamics: firm stock ≈ 80× more persistent than shocks | No entry / clone replacement (AGS) kills persistent variety loss | Indispensable |
| Separate vacancy and firm entry | u, v, N independent states (vs Schaal one-for-one) | Merging ties N and v one-for-one | Indispensable for O2 |
| Endogenous exit (p_0, ψ, ω_δ) | Intro sells it as recession amplification | Cutoff barely moves; Comparison B: endogenous-exit spec recovers *faster* | **Flag 1** (below) |
| Partial reposting (α) | Beveridge curve: reposting inflow is an identity | Without it cor(u,v) = +0.995 vs −0.80; only alternative (shrink s) contradicts data and CK | Indispensable for O3; **flag 4**: not yet coded/verified |

## Flag 1 — would estimating `dest_elast_target` make endogenous exit important?

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

## Other flags

2. **Convex vacancy cost:** role is persistence and the entry cushion, not signs. Check how
   much δ→v persistence it adds (Comparison D numbers).
3. **Variety effects** are the main departure from GS but moderate; report their contribution
   to Block M moments, not only to w^int. Also resolve the BK inconsistency (notes: ε > 1
   suffices; ξ×ε sweep: failures for ε < 3). ε > 2 is no longer needed by Prop 5.
4. **Partial reposting:** indispensable but unverified (R8, R9), identification vs δ_e (R2);
   turns the δ–s asymmetry into a quantitative result.
5. **Active core may be narrow:** δ→u peak 0.012 pp vs 1.71 data (partly D3 scale), half the
   amplification target, no hump. Documented bite: slot destruction + N persistence, partial
   reposting (pending), variety (moderate). Endogenous exit, convex vacancy cost, κ untested.

## Recommendation

After R8: one-at-a-time ingredient ablation on one calibration — ζ→0, ω_δ→0, α→0, ξ→∞,
κ = 0, σ = 0 — reporting Block M moments and δ IRFs. Task R14 in `../context/pending_tasks.md`.
