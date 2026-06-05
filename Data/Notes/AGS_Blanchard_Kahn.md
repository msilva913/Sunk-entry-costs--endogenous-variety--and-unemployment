# Blanchard-Kahn Conditions: AGS Model with $\sigma = 0$

**Author:** Mario Silva  
**Date:** June 4, 2026  
**Purpose:** Self-contained derivation of equilibrium uniqueness for the AGS (clone-replacement) model under log-utility / risk-neutrality ($\sigma = 0$). Not for inclusion in paper; for author's reference.

---

## 1. Setup

The AGS model fixes $N_t = 1$, $\rho = \mu = 1$ (no variety effects), and assumes exogenous exit ($p_0 = 0$, $\Lambda_t = 1$, $\delta_{e,t} = \delta_{t-1}$). The restriction $\sigma = 0$ sets the stochastic discount factor to $m_t = \beta$ (constant), eliminating consumption dynamics from the labor market block.

The relevant equilibrium conditions are:

**Job creation condition (JCC):**
$$\kappa + \frac{K_t}{q(\theta_t)} = \beta(1-\delta_t)\,\mathbb{E}_t\!\left[(1-\phi)(z_{t+1}-K_{t+1}-b) - \phi\theta_{t+1}(K_{t+1}+q_{t+1}\kappa) + (1-s_{t+1})\!\left(\kappa+\frac{K_{t+1}}{q_{t+1}}\right)\right]$$

**Flow vacancy value (K definition):**
$$K_t = Q_t - \beta(1-\delta_t)\,\mathbb{E}_t[Q_{t+1}]$$

**Entry:**
$$e_t = \left(\frac{Q_t}{x_m}\right)^\xi$$

**Vacancy LOM** (pre-committed stock + entry):
$$v_t = \underbrace{(1-\delta_{t-1})\bigl[(1-q_{t-1})v_{t-1} + s_{t-1}(1-u_{t-1})\bigr]}_{v_{pre,t}} + e_t$$

**Unemployment LOM** (standard DMP):
$$u_t = (1-f_{t-1})u_{t-1} + \tau_t\bigl[(1-u_{t-1})+f_{t-1}u_{t-1}\bigr]$$

At $\sigma = 0$, the resource constraint $C_t = z_t(1-u_t) - X_t$ determines consumption residually and does not feed back into the JCC or K definition.

---

## 2. State Space

**Predetermined variables** (known at the start of period $t$, before entry decisions):
- $u_t$: unemployment rate (from LOM, depends only on $t-1$ quantities)
- $v_{pre,t}$: pre-committed vacancy stock (determined before entry at $t$)
- Exogenous shocks: $z_t$, $\delta_t$, $s_t$ (each AR(1), realized at start of $t$)

**Jump variable** (one): $Q_t$, the vacancy value. Once $Q_t$ is determined, all remaining endogenous variables follow statically:
- $e_t = (Q_t/x_m)^\xi$
- $v_t = v_{pre,t} + e_t$
- $\theta_t = v_t / u_t$
- $K_t = Q_t - \beta(1-\delta_t)\mathbb{E}_t[Q_{t+1}]$ (requires knowing $\mathbb{E}_t[Q_{t+1}]$)

Since $K_t$ depends on $\mathbb{E}_t[Q_{t+1}]$, the forward-looking structure centers entirely on $Q_t$. The Blanchard-Kahn analysis requires exactly one eigenvalue outside the unit circle for this one jump variable.

---

## 3. Rewriting the JCC in Terms of Posting Cost $\Psi_t$

Define the **total posting cost per expected hire**:
$$\Psi_t \;\equiv\; \kappa + \frac{K_t}{q(\theta_t)}$$

From the K definition: $K_t = q_t(\Psi_t - \kappa)$.

Substitute $K_{t+1} = q_{t+1}(\Psi_{t+1}-\kappa)$ into the three terms on the JCC right-hand side:

**Flow-surplus term:**
$$(1-\phi)(z_{t+1}-K_{t+1}-b) = (1-\phi)(z_{t+1}-b+\kappa q_{t+1}) - (1-\phi)q_{t+1}\Psi_{t+1}$$

**Congestion term:**
$$-\phi\theta_{t+1}(K_{t+1}+q_{t+1}\kappa) = -\phi\theta_{t+1}q_{t+1}\Psi_{t+1} = -\phi f_{t+1}\Psi_{t+1}$$

where $f_{t+1} = \theta_{t+1}q_{t+1}$ is the job-finding rate.

**Continuation term:**
$$(1-s_{t+1})\Psi_{t+1}$$

Collecting, the JCC becomes:
$$\boxed{\Psi_t = \beta(1-\delta_t)\,\mathbb{E}_t\bigl[\pi_{t+1} + \gamma_{t+1}\Psi_{t+1}\bigr]}$$

where:
$$\pi_{t+1} \;\equiv\; (1-\phi)(z_{t+1}-b+\kappa q_{t+1}) \;\geq\; 0 \qquad \text{(flow match surplus)}$$
$$\gamma_{t+1} \;\equiv\; (1-s_{t+1}) - (1-\phi)q_{t+1} - \phi f_{t+1} \qquad \text{(continuation weight)}$$

**Economic reading.** $\Psi_t$ is the expected cost of creating one filled job: vacancy posting $\kappa$ plus the flow value of the vacancy $K_t/q_t$. Equation (JCC-$\Psi$) says this cost equals the discounted expected return, which has two parts: (i) the flow surplus $\pi_{t+1}$ earned when the job is filled, and (ii) the fraction $\gamma_{t+1}$ of the posting cost that "survives" to next period — accounting for the fact that some vacancies are filled (reducing the continuation cost), some matches are destroyed by separation, and hired workers contribute to the firm's going-concern value.

---

## 4. Steady State and the Key Coefficient

At the deterministic steady state (denoting SS values with overbars):
$$\bar\Psi = \beta(1-\bar\delta)\bigl[\bar\pi + \bar\gamma\,\bar\Psi\bigr]$$

Rearranging:
$$\bar\Psi\bigl[1 - \beta(1-\bar\delta)\bar\gamma\bigr] = \beta(1-\bar\delta)\bar\pi \tag{SS-JCC}$$

Define the **persistence coefficient**:
$$a \;\equiv\; \beta(1-\bar\delta)\bar\gamma$$

From (SS-JCC): $\bar\Psi(1-a) = \beta(1-\bar\delta)\bar\pi$.

**Signing $a$.**

*Upper bound ($a < 1$):* In any interior equilibrium, $\bar\Psi > 0$ (positive posting cost) and $\bar\pi > 0$ (positive match surplus — otherwise no job creation). Hence $1-a > 0$, so $a < 1$.

*Lower bound ($a > 0$):* Requires $\bar\gamma > 0$. Since $\bar\gamma = (1-\bar s) - (1-\phi)\bar q - \phi\bar f$, this holds when separation is not too frequent relative to job-creation flows — a mild condition satisfied at all standard calibrations. (At US calibration: $\bar s \approx 0.03$, so $1-\bar s \approx 0.97$, well above $(1-\phi)\bar q + \phi\bar f \approx 0.5$.)

Therefore $0 < a < 1$.

---

## 5. Linearization

Let $\tilde x_t \equiv x_t - \bar x$ denote the deviation of any variable from its steady state. Linearizing (JCC-$\Psi$) around the steady state:

$$\tilde\Psi_t = a\,\mathbb{E}_t[\tilde\Psi_{t+1}] + F_t \tag{Lin-JCC}$$

where the **forcing term** collects all terms that are functions of predetermined variables and shocks only:

$$F_t = \beta(1-\bar\delta)\,\mathbb{E}_t[\tilde\pi_{t+1}] + \beta(1-\bar\delta)\bar\Psi\,\mathbb{E}_t[\tilde\gamma_{t+1}] - \beta\tilde\delta_t\,(\bar\pi + \bar\gamma\bar\Psi)$$

The three components are:
1. $\beta(1-\bar\delta)\mathbb{E}_t[\tilde\pi_{t+1}]$: expected deviation in next-period flow surplus (driven by $z_{t+1}$, $\theta_{t+1}$)
2. $\beta(1-\bar\delta)\bar\Psi\,\mathbb{E}_t[\tilde\gamma_{t+1}]$: expected change in the continuation weight (driven by $s_{t+1}$, $\theta_{t+1}$)
3. $-\beta\tilde\delta_t(\bar\pi+\bar\gamma\bar\Psi)$: current-period effect of the destruction shock on the survival factor $1-\delta_t$

Note that $\tilde\pi_{t+1}$ and $\tilde\gamma_{t+1}$ each contain $\tilde\theta_{t+1}$ (endogenous) alongside shock terms. In a full linearized solution one would substitute $\tilde\theta_{t+1} = \tilde\theta(\tilde\Psi_{t+1}, \tilde u_{t+1}, ...)$ and collect terms — this does not change the eigenvalue $a$ established from the SS condition, because the SS condition already accounts for the equilibrium values $\bar q$, $\bar f$, $\bar\theta$.

---

## 6. Blanchard-Kahn Condition

Write (Lin-JCC) in **backward-map form**:
$$\mathbb{E}_t[\tilde\Psi_{t+1}] = \frac{1}{a}\,\tilde\Psi_t - \frac{1}{a}\,F_t$$

The scalar characteristic root of the backward map is:
$$\lambda = \frac{1}{a} > 1$$

**Blanchard-Kahn count:**

| | |
|---|---|
| Non-predetermined (jump) variables | 1 ($\tilde\Psi_t$, or equivalently $Q_t$) |
| Eigenvalues strictly outside unit circle | 1 ($\lambda = 1/a > 1$) |

The counts match: **the Blanchard-Kahn condition is satisfied**. The model has a unique locally bounded rational-expectations equilibrium in the neighborhood of the steady state.

---

## 7. Unique Bounded Solution

The unique non-explosive solution is obtained by iterating (Lin-JCC) forward and imposing $a^T\mathbb{E}_t[\tilde\Psi_{t+T}] \to 0$ as $T\to\infty$ (which holds since $a < 1$ and $\tilde\Psi$ is bounded in equilibrium):

$$\tilde\Psi_t = \sum_{j=0}^{\infty} a^j\,\mathbb{E}_t[F_{t+j}]$$

Expanding the dominant forcing term:

$$\tilde\Psi_t = \beta(1-\bar\delta)\sum_{j=0}^{\infty} a^j\,\mathbb{E}_t[\tilde\pi_{t+j+1}] + (\text{terms from }\tilde\delta,\,\tilde\gamma)$$

**Economic interpretation.** The total posting cost at date $t$ equals the present discounted value of expected future flow surpluses, discounted at rate $a = \beta(1-\bar\delta)\bar\gamma < 1$. The effective discount rate $a$ is smaller than $\beta(1-\bar\delta)$ because a fraction $1-\bar\gamma$ of the posting-cost "investment" is consumed each period by net job-creation flows — vacancies that are filled reduce the outstanding stock of posting obligations, and separations that occur allow costless reposting but destroy the accumulated match surplus. A larger steady-state separation rate $\bar s$, fill rate $\bar q$, or job-finding rate $\bar f$ all increase the rate at which posting costs are consumed, lowering $\bar\gamma$ and hence $a$, and making the forward sum converge faster.

The solution $Q_t$ (vacancy value) is then uniquely recovered: since $K_t = \Psi_t q_t - \kappa q_t$ and $K_t = Q_t - \beta(1-\delta_t)\mathbb{E}_t[Q_{t+1}]$, the sequence $\{Q_t\}$ is pinned by the unique $\{\Psi_t\}$ together with the initial condition.

---

## 8. Why $\sigma > 0$ Breaks Analytical Tractability

With $\sigma > 0$, the stochastic discount factor is $m_t = \beta(C_t/C_{t-1})^{-\sigma}$, which introduces $C_t$ into the JCC via the discounting of future surpluses. Since $C_t = z_t(1-u_t) - X_t$ and $X_t = \kappa v_t$ depends on $e_t = G(Q_t)$, the scalar reduction of Step 3 no longer holds. The system becomes at minimum two-dimensional in $(\tilde\Psi_t, \tilde C_t)$, and signing both eigenvalues analytically requires tracking the interaction between the labor market block and the household Euler equation. Numerical verification via the solved linear system is the standard approach.

---

## 9. Role of $\kappa$ and Entry Elasticity $\xi$

**Fixed matching cost $\kappa$:** Enters $\bar\pi = (1-\phi)(\bar z - \bar b + \kappa\bar q)$ and $\bar\Psi = \kappa + \bar K/\bar q$. Both shift the level of $a$ through the ratio $\beta(1-\bar\delta)\bar\pi/\bar\Psi$. However, the key condition $a < 1$ follows from $\bar\pi > 0$ regardless of $\kappa \geq 0$. **$\kappa$ plays no role in whether the BK condition holds.**

**Entry elasticity $\xi$:** With $\xi\to\infty$ (free entry), $Q_t = x_m$ is a constant pinned by the entry-cost upper bound, and there is no forward-looking variable in the system — the JCC becomes a static condition pinning $\theta_t$ recursively each period and BK analysis does not apply. For finite $\xi$, $Q_t$ adjusts dynamically and the eigenvalue analysis above applies. The condition $a < 1$ holds for all $\xi \in (0,\infty)$, and as $\xi\to\infty$ the amplitude of $Q$ fluctuations vanishes while the eigenvalue $1/a$ remains unchanged. **$\xi$ governs the amplitude of the response, not the existence or uniqueness of equilibrium.**
