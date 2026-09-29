# The recruiter–retailer segmentation: why it is the right structure
**Written:** September 29, 2026 · **Branch:** `costly_vacancy_reposting`
**Status:** modeling-judgment note. Records why the agent segmentation is the right choice for
this paper, its one substantive cost, and how to frame it — especially now that the reposting
margin makes the cost salient.

## 1. What the segmentation is

Retailers are monopolistically competitive firms, one per product line, facing a downward-sloping
demand curve ρ^d(y). They rent labor services in a competitive market at price w_int and choose
employment flexibly (`eq:firm_bellman`, `eq:retailer_optimality`). Recruiters are a separate unit
measure of agents who create vacancies at a sunk cost x ~ G, match with workers, pay the fixed
matching cost κ on filling, and then rent the matched labor to retailers at w_int
(`eq:value_recruiter_unmatched`, `eq:value_recruiter_matched`). Nash bargaining is between the
recruiter and the worker, over w_int taken as given.

## 2. Why it is the right choice — arguably forced

The paper needs two things that are hard to hold together.

1. **Monopolistic competition with diminishing marginal revenue.** Required for the variety
   channel — ρ(N), markups, the N → ρ → w_int → θ externality. Non-negotiable; it is a core
   contribution.
2. **A standard, recognizable DMP labor market.** The clean job creation condition and Nash wage
   that make the Beveridge-curve and amplification results legible.

These conflict. If workers bargained directly with a monopolistically competitive firm facing a
downward-sloping demand curve, the result is Stole–Zwiebel intra-firm bargaining: the firm
strategically over-hires to depress its own marginal revenue and improve its bargaining position,
employment becomes a fixed point, and the wage equation stops looking like DMP. That machinery
would bury the mechanisms the paper is about.

The recruiter segmentation is the **minimal** device that reconciles the two. Inserting a
competitive market for labor services at w_int makes the retailer's marginal revenue product a
scalar the recruiter takes as given, so the recruiter–worker match splits a well-defined surplus
and the model recovers textbook DMP. Every alternative is worse **for this paper**:

- **Stole–Zwiebel directly:** intractable with variety plus endogenous exit.
- **One worker per firm (Shao–Silos):** kills the intensive margin and — fatally — collapses the
  δ-versus-s distinction that is the paper. The draft already criticizes it on empirical grounds.
- **Large-firm search with firm-level employment (Elsby–Michaels, Acemoglu–Hawkins):** more
  realistic, and it would let exit respond to staffing, but it reintroduces firm size as a state
  and fights the DS-CES symmetric aggregation the variety channel relies on. Too heavy, and a
  different paper.

So the segmentation is not a gratuitous trick. It is the standard Christiano-style device
deployed to buy the one thing that otherwise breaks, and referees will recognize it.

## 3. Terminology: "recruiters," not "wholesalers"

We call these agents recruiters because that is precisely their function: they create vacancies,
run the search-and-match, and supply the matched labor. This is more natural than Christiano et
al.'s relabeling of an intermediate good transformed into a final good — the recruiter name says
what the agent does, and it makes the search friction sit where it belongs. Keep it. The segment
is the firm's hiring arm, not a separate species; §5 below aligns the prose accordingly.

## 4. The one substantive cost — firm exit is independent of labor frictions

The reposting cost lands on recruiter profits Π^int = (w_int − w)L, not retail profits
Π^f = Y_c/ε − X_c, and firm exit depends only on Π^f (`eq:firm_bellman`, `eq:cutoff`). So a firm
that finds it expensive to staff its positions does **not** become more likely to exit. The
retailer's value ν_f, its cutoff χ^c, and its dividend d_f carry no vacancy or reposting term;
the dependence is triangular — the recruiter's Q, J, and the reposting option depend on retailer
survival (1−δ)F(χ^c), but not the reverse.

This is a genuine assumption, and the reposting margin makes it sharp. A skeptic can ask: if
reposting is costly enough that firms destroy positions, would that not also push marginal firms
out? The model says no, by construction.

## 5. Why that cost is a feature for this paper

For this paper the independence is arguably a feature, not a bug. The contribution is decomposing
separations into a product-line-destruction channel (δ) and a labor channel (s / reposting). If
firm exit responded to staffing frictions, the δ and s channels would entangle and the clean
decomposition the paper is selling would blur. The segmentation **protects** the separate
identification of the two margins. The limitation lines up with the paper's goals rather than
fighting them.

## 6. What to do — framing, not structure

Do not change the model. Change how it is defended.

1. **State the independence assumption explicitly and own it.** It is currently an implicit
   consequence of the segmentation. Say plainly: firm exit is driven by product-line viability,
   deliberately insulated from labor-market frictions, to keep the δ and s margins separately
   identified; a richer model where exit responds to staffing costs would entangle them and is
   left for future work. This turns a silent assumption into a defended choice.
2. **Align the reposting narrative with the mechanics.** The "missing middle" prose says a
   *surviving firm* cuts a position; in the model the decision sits with the *recruiter* for that
   product line. Economically the recruiter is the firm's hiring arm, but make the identification
   explicit so a reader does not wonder whether "firm" and "recruiter" pull in different
   directions. This matters more now that reposting is a headline mechanism.
3. **Acknowledge the richer alternative.** Large-firm search where exit responds to staffing is
   the natural extension; note it as future work, and note why it is out of scope here (firm size
   as a state fights the symmetric variety aggregation).

## 7. Why this became salient now

Before the reposting margin, the segmentation was a pure tractability device and the
independence-of-exit assumption was invisible. The reposting margin is a labor-side friction that
destroys positions at surviving firms, so it puts pressure exactly on the exit-versus-staffing
link the segmentation severs. The right response is not to re-couple them — that would cost the
decomposition — but to state the separation as a deliberate design choice. See
[`delta_calibration_and_the_reposting_margin.md`](delta_calibration_and_the_reposting_margin.md)
for the reposting margin and [`../context/model_equations.md`](../context/model_equations.md) f[1],
f[3] for the confirmation that reposting stays quarantined in the recruiter block.
