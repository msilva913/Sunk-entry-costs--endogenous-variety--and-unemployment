# Nesting Map — special cases, limits, and isomorphisms
**Created:** October 13, 2026 · **Branch:** `costly_vacancy_reposting` · **Draft as of:** `5ed1ab1`

This file owns **how the baseline relates to DMP, Coles–Kelishomi (CK), Gabrovski–Silva (GS),
and Bilbiie–Ghironi–Melitz (BGM)**: which restrictions produce which economy, what kind of
equivalence holds, and what is lost. It does **not** own proof status or audit history (that is
[`draft_status.md`](draft_status.md) "Propositions") or numerical results
([`findings.md`](findings.md)). Labels are `Draft.tex` labels; theorem numbers are as of
`5ed1ab1` and drift, so cite by label.

**Terminology (MS, Oct 13).** Bare "free entry" means product lines (`eq:free_entry`). For
vacancies (ξ → ∞, Q = x_m), say "vacancy free entry" or "infinitely elastic vacancy creation".

---

## 1. The chain

```
Baseline (DS-CES or general CES; σ, ζ, p_0, α, ξ, κ)
 │
 ├─ σ = 0, ζ = 0, p_0 = 0 ─────────────► labor block closes on its own        prop:independence
 │
 ├─ p_0 = 0, {L_t, X_t} taken as given ► business block ≅ LR-BGM               rem:bgm_nest
 │
 └─ clone replacement (N_t = 1, μ = 1, p_0 = 0) ► AGS                         prop:ags
      ▲ justified as a limit: DS-CES, p_0 = 0, ε → ∞ with f_e = f̄_e/ε         prop:double_limit
      │
      AGS ── σ = κ = 0 ──► GS + aggregate AR(1) shocks to δ_t, s_t            rem:nesting
              │
              ├─ also s_t = 0 ──► CK (w^int = z)                              rem:nesting
              │
              └─ ξ → ∞, e_t > 0 ──► DMP with flow user cost K_t on every      rem:dmp
                    │                position (vacancy cost K_t, output z − K_t)
                    └─ δ constant ──► standard DMP (relabeled productivity)   rem:dmp
```

The analytical δ–s results (`lem:theta`, `prop:fe`) are proved in the GS node at ξ → ∞
(`ass:vacancy_fe`), which `prop:double_limit` makes a limit of the baseline rather than a
different model.

---

## 2. One row per link

| Link | Restrictions | Equivalence | What survives | What is lost | Label · proof · check |
|---|---|---|---|---|---|
| Baseline → separable labor block | σ = 0 (m = β); ζ = 0 (ρ ≡ 1); p_0 = 0 (Λ = Λ_r = 1, X^r = 0) | **Exact, iff** "for all values of the remaining parameters" | Destruction flow, entry of vacancies, wage bargaining; business block still runs, taking u, X as given | Consumption in the SDF, variety in the JCC, endogenous exit, partial reposting | `prop:independence` · App. G · — |
| Baseline → AGS | Clone replacement: destroyed lines replaced at no cost within the period; N_t = 1, μ = 1 (so ρ = μ = 1), p_0 = 0 | **Construction** (a heuristic on its own) | Destruction as a shock to jobs *and* vacancies; δ_{e,t} = δ_{t−1} | Business formation, variety, exit, reactivation, the pricing of replacement, measured ≠ structural destruction | `prop:ags` · App. G (BK verified numerically at σ = κ = 0) · `gs_finite_xi_entry_offset.jl` asserts 2 stable roots |
| Baseline ⇒ AGS as a limit | DS-CES, p_0 = 0, ε → ∞ **jointly with** f_e = f̄_e/ε | **Limit of equilibrium conditions, for allocations** (N_t still moves via `eq:N_rule_limit` but feeds nothing back) | N̄ = z̄L̄(1−δ̄)/[f̄_e(r+δ̄)] > 0; w^int → z; allocations independent of f̄_e, so N̄ = 1 is w.l.o.g. | Firm value, investment and entry labor vanish at rate 1/ε | `prop:double_limit` · App. G · `double_limit_check.jl` |
| AGS → GS | σ = κ = 0 | **Isomorphism**, extended with AR(1) shocks to δ_t and s_t (GS has constant δ, s) | δ/s split, sunk durable vacancies, predetermined vacancy stock | — | `rem:nesting` · — · — |
| GS → CK | additionally s_t = 0 | **Exact**, with w^int = z | One separation margin that destroys positions | The s margin (costless reposting) | `rem:nesting` · — · — |
| GS at ξ → ∞ → DMP | `ass:vacancy_fe`: Q = x_m, **e_t > 0**, 0 < δ̄ ≤ τ̄, (1−δ̄)f̄ + τ̄ < 1; σ = 0 | **Isomorphism of the labor market block** with a DMP economy (same timing and margins) in which every position bears K_t = x_m(r+δ_t)/(1+r): vacancy cost K_t, match output z_t − K_t. Wage and JCC coincide term by term; v_t = θ_tu_t jumps; v_pre irrelevant | θ, u, v paths | Holds only while e_t ≥ 0 is slack. Resource-cost timing differs (x_m·e_t vs K·v_t), so σ = 0 is needed | `rem:dmp` · algebra in the remark · — |
| → standard DMP | additionally δ constant | **Exact** (DMP with relabeled productivity) | z and s shocks transmit as in DMP | Under δ shocks the counterpart's cost is indexed to δ_t; that indexation **is** the user-cost channel of `prop:fe` | `rem:dmp` · — · flow-cost comparator (K fixed) in `beveridge_free_entry.jl` |
| Baseline business block → LR-BGM | p_0 = 0 (X^c = 0, Λ = 1); LR-BGM = BGM with labor dropped from preferences, {L_t}, {X_t} exogenous | **Isomorphism for any given sequences** {L_t, X_t}, at the level of GDP (`eq:gdp`), not gross output | Entry, variety, the Euler equation for lines | Endogeneity of L_t and X_t (set by the labor block); no intensive labor margin, so BGM itself is not nested | `rem:bgm_nest` · App. G · — |
| DS-CES steady state in the double limit | DS-CES, p_0 → 0, double limit | JC curve vertical; **unique** steady state iff z > Φ(0) (away from the limit: 0 or generically 2, under `ass:matching` + (C)) | — | — | `prop:equilibria` · App. J · `equilibria_convexity_check.jl` |

---

## 3. Caveats that apply across rows

1. **Interior entry.** `rem:dmp`, `lem:theta` and `prop:fe` are local and assume e_t > 0. A large
   δ shock lowers desired vacancies by more than destruction lowers the inherited stock, so
   entry stops, Q_t < x_m, and v_pre becomes a state (irreversibility of vacancy investment).
   No draft number for the binding threshold yet: the 9.3 % of δ̄ in `beveridge_free_entry.txt`
   is for the x_v = 0.5 economy, not κ = 0 (N15). See R17 in [`pending_tasks.md`](pending_tasks.md).
2. **σ = 0** is needed for the separable labor block, for `rem:nesting`'s GS isomorphism, and for
   `rem:dmp`.
3. **The double limit concerns equilibrium conditions, not convergence of impulse responses.**
   The free-entry results describe the limit economy, not the rate at which the baseline
   approaches it.
4. **The CK parameterization ≠ the CK model.** "δ̄ = τ̄" inside our model (all separations destroy
   positions) is CK's separation structure, used at any ξ. CK's own model also has finite ξ
   ("Diamond entry"); their ξ = ∞ is "analogous to the free entry case" (CK p. 125,
   doi:10.1257/mac.20150040). CK's free-entry benchmark corresponds to standard DMP (last two
   rows), not to the durable-vacancy ξ → ∞ case.
5. **The single limit ε → ∞ with f_e fixed is degenerate** (profits vanish, N → 0). Only the joint
   limit works; δ → 0 would restore N but removes the destruction margin.
6. **Exit and reactivation are both removed by p_0 = 0.** To switch them separately, use ω_δ → 0
   (exit) and α → 0 (reactivation); see §5.3's comparison list and R9.

## 4. Sources
- CK: AEJ: Macro 2018, doi:10.1257/mac.20150040 (`../Key papers/markdown/Coles_Kelishomi_2018.md`).
- GS: JEDC 2025, doi:10.1016/j.jedc.2025.105205 (`../Key papers/markdown/Gabrovski_Silva_JEDC.md`).
- BGM: JPE 2012, doi:10.1086/665825 (`../Key papers/markdown/Bilbiie_Ghironi_Melitz_2012.md`).
