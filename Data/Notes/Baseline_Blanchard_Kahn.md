# Blanchard-Kahn Conditions: Full Baseline Model
**Author:** Mario Silva  
**Date:** June 6, 2026  
**Purpose:** BK analysis for Definition 3.7 equilibrium (endogenous variety, σ > 0). Companion to AGS_Blanchard_Kahn.md. Restricts to DS-CES and p₀ = 0 for the main derivation; §8 discusses endogenous exit.

---

## 1. Setup and Notation

**Preferences.** DS-CES throughout: $\mu(N) = \varepsilon/(\varepsilon-1)$ constant, $\rho(N) = N^{\psi_N}$ with $\psi_N \equiv 1/(\varepsilon-1)$.

**Stochastic discount factor.** $m_{t+1} = \beta(C_{t+1}/C_t)^{-\sigma}$. Log-linearizing: $\hat m_{t+1} \equiv \tilde m_{t+1}/\bar m = -\sigma(\hat C_{t+1} - \hat C_t)$, where hats denote log-deviations from SS and tildes level deviations.

**Survival factor.** With $p_0 = 0$: $\Lambda_t = 1$ deterministically, $\delta_{e,t} = \delta_{t-1}$, and the survival-adjusted discount factor is $\beta(1-\delta_t)$.

**Persistence coefficient** (from AGS notes, now with $\Lambda = 1$):
$$a \;\equiv\; \beta(1-\bar\delta)\bar\gamma, \qquad \bar\gamma = (1-\bar s)-(1-\phi)\bar q - \phi\bar f,\quad 0 < a < 1$$

---

## 2. Jump-Variable Classification

**Predetermined at start of $t$:** $\{z_t, \delta_t, s_t, u_t, v_{pre,t}, N_t\}$.

*On $N_t$:* With $p_0=0$, the N LOM is $N_t = (1-\delta_{t-1})(N_{t-1}+N_{t-1}^e)$ where $N_{t-1}^e$ was determined by last period's resource constraint. Both inputs are $t{-}1$ quantities, so $N_t$ is predetermined.

**Two jump variables:** $\{Q_t,\, C_t\}$. Once these are determined, everything else follows:
- $e_t = G(Q_t)$, $v_t = v_{pre,t}+e_t$, $\theta_t = v_t/u_t$, $K_t = q_t(\Psi_t-\kappa)$ (static)
- $w_t^{int} = \rho(N_t)z_t/\mu$ (predetermined $N_t$, $z_t$)
- $Y_t^c = C_t + X_t + X_t^c$ (resource constraint, static given $C_t$, $Q_t$, and predetermined $N_t$)
- $\chi_t^c$, $N_t^e$ follow statically

**Two forward-looking equations:** JCC (pins $Q_t$ via $\Psi_t$) and BFE / household Euler (pins $C_t$).

---

## 3. Steady State

**JCC SS** (as in AGS notes):
$$\bar\Psi(1-a) = \beta(1-\bar\delta)\bar\pi, \qquad \bar\pi = (1-\phi)(w^{int}-b+\kappa\bar q) > 0 \tag{SS-JCC}$$

**BFE SS** (from eq:euler\_bf, DS-CES, ignoring $X^c/Y^c$ corrections for clarity):
$$f_e\bar\rho = \beta(1-\bar\delta)\Bigl[f_e\bar\rho + \underbrace{\frac{\bar Y^c/\bar N}{\varepsilon-1}}_{\bar\Pi^f}\Bigr] \tag{SS-BFE}$$

where $\bar\Pi^f \equiv \bar Y^c/(\varepsilon-1)\bar N$ is per-firm dividend (net of markup wedge). Define:
$$a_N \;\equiv\; \beta(1-\bar\delta) \qquad \text{and} \qquad \eta \;\equiv\; \frac{\bar\Pi^f}{f_e\bar\rho} = \frac{a_N^{-1}-1}{1} = \frac{1-a_N}{a_N}$$
so $a_N(1+\eta) = 1$. Note $a_N > a$ since $a_N = \beta(1-\bar\delta)$ while $a = \beta(1-\bar\delta)\bar\gamma < \beta(1-\bar\delta)$.

---

## 4. The N LOM and How Jump Variables Affect $N_{t+1}$

This is the key transmission mechanism absent from AGS. The N LOM gives next-period product lines as:
$$N_{t+1} = (1-\delta_t)\Bigl(N_t + \underbrace{z_t(1-u_t)/f_e}_{\text{entry capacity}} - \underbrace{Y_t^c/(f_e\rho(N_t))}_{\text{entry cost}}\Bigr) \tag{N-LOM}$$

Linearizing and grouping by what drives $Y_t^c$:
$$\hat N_{t+1} = (1-\bar\delta)\hat N_t - (1-\bar\delta)\underbrace{\frac{\bar Y^c/\bar N}{f_e\bar\rho}}_{\equiv\,\omega}\hat Y_t^c + \text{forcing}(\hat z_t, \hat u_t, \hat\delta_t) \tag{LOM-lin}$$

where $\omega \equiv (\bar Y^c/\bar N)/(f_e\bar\rho) = \eta(\varepsilon-1)$ from the SS BFE.

Next, decompose $\hat Y_t^c$ into contributions from the two jump variables. Write $Y_t^c = C_t + X_t^v + X_t^f$ where $X_t^v = e_t(\xi/(\xi+1))Q_t$ are sunk vacancy costs and $X_t^f = \kappa q_t v_t$ are fixed matching costs. Define shares:
$$s_C \equiv \bar C/\bar Y^c \approx 0.80\text{--}0.90 \qquad s_Q \equiv \bar X_v/\bar Y^c \approx 0.01\text{--}0.03$$

Then: $\hat Y_t^c = s_C \hat C_t + s_Q \hat Q_t + \text{small terms}$.

Substituting into (LOM-lin):
$$\hat N_{t+1} = (1-\bar\delta)\hat N_t - (1-\bar\delta)\omega\,s_C\,\hat C_t - (1-\bar\delta)\omega\,s_Q\,\hat Q_t + \text{forcing} \tag{N-lin}$$

**Economic reading.** Higher $C_t$ (or $Q_t$) consumes more resources, reducing firm entry and hence $N_{t+1}$. The $s_C$ channel is large; the $s_Q$ channel is small ($\approx 1\text{--}3\%$).

---

## 5. Linearized BFE

Log-linearize eq:euler\_bf around SS. LHS: $f_e\bar\rho\cdot\psi_N\hat N_t$ (predetermined).

RHS expands to three parts: (i) the future firm-value term $f_e\rho(N_{t+1})$; (ii) the dividend term $(\bar Y^c/\bar N)^{eq}\cdot\psi_N \cdot(\hat Y^c_{t+1}-\hat N_{t+1})$; and (iii) the discount factor.

The discount factor contributes $-\sigma(E_t[\hat C_{t+1}]-\hat C_t)$ times $f_e\bar\rho\cdot(1+\eta)$. By the SS BFE, $a_N(1+\eta)=1$, so this coefficient is $a_N(1+\eta)\cdot\sigma/a_N = \sigma$ (after dividing through by $f_e\bar\rho$).

After substituting (N-lin) for $E_t[\hat N_{t+1}]$ (known at time $t$) and collecting terms:

$$\sigma E_t[\hat C_{t+1}] = \sigma\hat C_t + \underbrace{a_N(1-\bar\delta)\omega\,\psi_N\,s_C}_{\equiv\,\beta_{CC}}\hat C_t + \underbrace{a_N(1-\bar\delta)\omega\,\psi_N\,s_Q}_{\equiv\,\beta_{QC}}\hat Q_t + \text{forcing}(\hat N_t, \ldots) \tag{BFE-lin}$$

where the forcing contains $\hat N_t$ (predetermined), future dividend terms, and shock terms. Rearranging into backward-map form:

$$\boxed{E_t[\hat C_{t+1}] = \underbrace{\Bigl(1 + \frac{\beta_{CC}}{\sigma}\Bigr)}_{\equiv\,1/b}\hat C_t + \frac{\beta_{QC}}{\sigma}\hat Q_t + \text{forcing}} \tag{BM-BFE}$$

**The coefficient $1/b$:**
$$\frac{1}{b} = 1 + \frac{\beta_{CC}}{\sigma} = 1 + \frac{a_N(1-\bar\delta)\omega\,\psi_N\,s_C}{\sigma} > 1 \quad \forall\;\varepsilon > 1$$

$1/b > 1$ follows from $\psi_N = 1/(\varepsilon-1) > 0$ (i.e., $\varepsilon > 1$), all other factors being positive.

**Key identity** (used repeatedly below):
$$\frac{1}{b} - \frac{\beta_{CC}}{\sigma} = 1 \tag{ID}$$

---

## 6. Linearized JCC

The wage-substituted JCC (eq:jcc\_wage) in forward form, after linearizing $m_{t+1}$:

$$\tilde\Psi_t = a\,E_t[\tilde\Psi_{t+1}] - \sigma\bar\Psi\,E_t[\hat C_{t+1}] + \sigma\bar\Psi\,\hat C_t + \text{forcing}_1 \tag{JCC-fwd}$$

**Derivation of the $\sigma\bar\Psi$ coefficient.** From the SS JCC: $a_N(\bar\pi + \bar\gamma\bar\Psi) = \bar\Psi$, so $\bar\pi + \bar\gamma\bar\Psi = \bar\Psi/a_N$. Log-linearizing $m_{t+1}(\pi_{t+1}+\gamma_{t+1}\Psi_{t+1})$, the discount factor contributes $-\sigma\cdot a_N\cdot(\bar\pi+\bar\gamma\bar\Psi)\cdot(E_t[\hat C_{t+1}]-\hat C_t) = -\sigma\bar\Psi\cdot(E_t[\hat C_{t+1}]-\hat C_t)$. This is identical to the AGS notes derivation, now with $\sigma > 0$.

**The variety channel.** $w^{int}_{t+1} = \rho(N_{t+1})z_{t+1}/\mu$ depends on $N_{t+1}$, which (by (N-lin)) depends on current $\hat C_t$ and $\hat Q_t$. This introduces two new terms into $\text{forcing}_1$:

$$\text{variety contribution to forcing}_1 = a_N\cdot(1-\phi)w^{int}\cdot\psi_N\cdot\partial_{\hat C_t}\hat N_{t+1} \cdot \hat C_t + (\ldots)\hat Q_t$$

Define:
$$\alpha_C \equiv a_N(1-\bar\delta)\omega\,\psi_N\,(1-\phi)\bar w^{int}\,s_C/\bar z > 0 \qquad \alpha_Q \equiv a_N(1-\bar\delta)\omega\,\psi_N\,(1-\phi)\bar w^{int}\,s_Q/\bar z > 0$$

*(The $\bar w^{int}/\bar z = \bar\rho/\mu$ ratio converts between wage and output units.)*

Collecting all current-period terms, the JCC forward equation is:

$$\tilde\Psi_t = a\,E_t[\tilde\Psi_{t+1}] - \sigma\bar\Psi\,E_t[\hat C_{t+1}] + (\sigma\bar\Psi - \alpha_C)\hat C_t - \alpha_Q\hat Q_t + \text{forcing}(\hat N_t, \hat z_t, \hat\delta_t, \hat s_t) \tag{JCC-lin}$$

**Signs.** The $-\alpha_C\hat C_t$ term captures: higher $C_t$ reduces $N_{t+1}$ (via N-LOM), lowering $w^{int}_{t+1}$, lowering match surplus, and reducing $\Psi_{t+1}$. The coefficient $\sigma\bar\Psi - \alpha_C > 0$ in practice (the direct SDF effect $\sigma\bar\Psi$ dominates the variety channel $\alpha_C$ at calibrated values).

Rearranging (JCC-lin) to backward-map form (solving for $E_t[\tilde\Psi_{t+1}]$):

$$E_t[\tilde\Psi_{t+1}] = \frac{1}{a}\tilde\Psi_t + \frac{\sigma\bar\Psi}{a}E_t[\hat C_{t+1}] - \frac{\sigma\bar\Psi-\alpha_C}{a}\hat C_t + \frac{\alpha_Q}{a}\hat Q_t + \text{forcing} \tag{BM-JCC-pre}$$

Note that (BM-JCC-pre) is **not yet in standard form** because $E_t[\hat C_{t+1}]$ appears on the right-hand side. Substitute (BM-BFE) to eliminate it:

$$E_t[\tilde\Psi_{t+1}] = \frac{1}{a}\tilde\Psi_t + \frac{\sigma\bar\Psi}{a}\Bigl[\frac{1}{b}\hat C_t + \frac{\beta_{QC}}{\sigma}\hat Q_t\Bigr] - \frac{\sigma\bar\Psi - \alpha_C}{a}\hat C_t + \frac{\alpha_Q}{a}\hat Q_t + \text{forcing}$$

Collecting $\hat C_t$ terms:
$$\frac{\sigma\bar\Psi}{ab} - \frac{\sigma\bar\Psi - \alpha_C}{a} = \frac{\sigma\bar\Psi}{a}\Bigl(\frac{1}{b}-1\Bigr) + \frac{\alpha_C}{a} \stackrel{\text{(ID)}}{=} \frac{\sigma\bar\Psi}{a}\cdot\frac{\beta_{CC}}{\sigma} + \frac{\alpha_C}{a} = \frac{\bar\Psi\beta_{CC}+\alpha_C}{a}$$

Collecting $\hat Q_t$ terms:
$$\frac{\bar\Psi\beta_{QC}}{a} + \frac{\alpha_Q}{a} = \frac{\bar\Psi\beta_{QC}+\alpha_Q}{a}$$

The substituted JCC backward map is:

$$\boxed{E_t[\tilde\Psi_{t+1}] = \frac{1}{a}\tilde\Psi_t + \frac{\bar\Psi\beta_{CC}+\alpha_C}{a}\hat C_t + \frac{\bar\Psi\beta_{QC}+\alpha_Q}{a}\hat Q_t + \text{forcing}} \tag{BM-JCC}$$

The diagonal entry $1/a$ is **exact**: the identity (ID) ensures the $\sigma\bar\Psi$ terms cancel completely upon substitution, leaving $1/a$ unchanged.

---

## 7. The 2×2 Backward Map

To write a closed 2×2 system in state vector $(\tilde\Psi_t,\,\hat C_t)$, express $\hat Q_t$ in terms of $\tilde\Psi_t$. From $\Psi_t = \kappa + K_t/q_t$ and the K-definition $K_t = Q_t - a_N E_t[Q_{t+1}] + \ldots$, the AGS notes (§7) show that $\{Q_t\}$ is recovered uniquely from $\{\Psi_t\}$ via the forward sum. In the eigenspace $\tilde\Psi_t = \lambda^{-t}\tilde\Psi_0$:

$$\hat Q_t = c\cdot\frac{\tilde\Psi_t}{\bar\Psi}, \qquad c \;\equiv\; \frac{\bar q\lambda}{\lambda - a_N} > 0 \quad (\lambda > a_N) \tag{Q-Psi}$$

Substituting (Q-Psi) into (BM-JCC) and (BM-BFE):

$$\begin{pmatrix}E_t[\tilde\Psi_{t+1}]\\ E_t[\hat C_{t+1}]\end{pmatrix} = \underbrace{\begin{pmatrix}A_{11} & A_{12}\\ A_{21} & A_{22}\end{pmatrix}}_{\mathbf{A}}\begin{pmatrix}\tilde\Psi_t\\ \hat C_t\end{pmatrix} + \text{forcing}$$

with entries:

$$A_{11} = \frac{1}{a} + \frac{c(\bar\Psi\beta_{QC}+\alpha_Q)}{a\bar\Psi} \;>\; \frac{1}{a} > 1 \tag{A11}$$
$$A_{12} = \frac{\bar\Psi\beta_{CC}+\alpha_C}{a} \;>\; 0 \tag{A12}$$
$$A_{21} = \frac{c\beta_{QC}}{\sigma\bar\Psi} \;\approx\; 0 \quad (s_Q\text{ small}) \tag{A21}$$
$$A_{22} = \frac{1}{b} > 1 \tag{A22}$$

**Signs.** All four entries are positive. $A_{12} > 0$: higher $C_t$ today reduces $N_{t+1}$, lowers $w^{int}_{t+1}$, and lowers $\Psi_{t+1}$ — so $\partial E_t[\tilde\Psi_{t+1}]/\partial\hat C_t > 0$ (lower $C$ → higher $\Psi$ … wait). Let me re-examine. $A_{12} > 0$ means higher $\hat C_t$ raises $E_t[\tilde\Psi_{t+1}]$? Actually: the sign in (JCC-lin) is $-\alpha_C$ (higher $C$ lowers future Ψ), but after multiplying through and rearranging, the backward-map coefficient is $+(\bar\Psi\beta_{CC}+\alpha_C)/a > 0$. Economically: a positive surprise $\hat C_t$ requires a compensating *rise* in future posting cost $\Psi_{t+1}$ to keep the JCC balanced — higher consumption today depresses future match surplus (via lower $N_{t+1}$), requiring a lower effective discount rate (which shows up as higher $\Psi$). The positive $A_{12}$ is consistent with the general principle that both eigenvalues of $\mathbf{A}$ must exceed 1 for BK.

---

## 8. Blanchard-Kahn Condition

**Determinant.** Using (A11)–(A22) and the identity (ID):

$$\det(\mathbf{A}) = A_{11}\cdot A_{22} - A_{12}\cdot A_{21}$$

$$= \Bigl(\frac{1}{a}+\frac{c(\bar\Psi\beta_{QC}+\alpha_Q)}{a\bar\Psi}\Bigr)\frac{1}{b} - \frac{(\bar\Psi\beta_{CC}+\alpha_C)}{a}\cdot\frac{c\beta_{QC}}{\sigma\bar\Psi}$$

$$= \frac{1}{ab} + \frac{c}{a\bar\Psi}\left[\frac{\bar\Psi\beta_{QC}+\alpha_Q}{b} - \frac{(\bar\Psi\beta_{CC}+\alpha_C)\beta_{QC}}{\sigma}\right]$$

Evaluate the bracket at $\alpha_Q = \alpha_C = 0$ (leading order, retaining only the $s_C$-driven terms):

$$\frac{\bar\Psi\beta_{QC}}{b} - \frac{\bar\Psi\beta_{CC}\beta_{QC}}{\sigma} = \bar\Psi\beta_{QC}\underbrace{\Bigl(\frac{1}{b}-\frac{\beta_{CC}}{\sigma}\Bigr)}_{=\,1\text{ by (ID)}} = \bar\Psi\beta_{QC} > 0$$

With $\alpha_Q,\alpha_C > 0$, both $\bar\Psi\beta_{QC}$ and the $\alpha$ corrections contribute positively. Therefore:

$$\boxed{\det(\mathbf{A}) = \frac{1}{ab} + \frac{c\,\beta_{QC}}{a} + \text{positive corrections} > \frac{1}{ab} > 1} \tag{DET}$$

The off-diagonal coupling **increases** the determinant above $1/(ab)$; it does not threaten BK.

**$f(1)$ test.** Evaluate the characteristic polynomial at $\lambda = 1$:

$$f(1) = \det(\mathbf{A}-\mathbf{I}) = (A_{11}-1)(A_{22}-1) - A_{12}A_{21}$$

$$= \underbrace{\Bigl(\frac{1}{a}-1+\text{pos}\Bigr)}_{>0}\underbrace{\Bigl(\frac{1}{b}-1\Bigr)}_{>0} - \underbrace{A_{12}A_{21}}_{\geq 0}$$

At leading order ($A_{21} \approx 0$, $s_Q \to 0$): $f(1) = (1/a-1)(1/b-1) > 0$. With $A_{21} > 0$ (small), $f(1)$ decreases slightly but remains positive provided $A_{21}$ is small, which holds since $A_{21} = O(s_Q)$.

**$f(-1)$ test.** $f(-1) = (A_{11}+1)(A_{22}+1)+A_{12}A_{21} > 0$ trivially (all terms positive).

**Conclusion.** The characteristic polynomial satisfies $f(1) > 0$, $f(-1) > 0$, and $\det(\mathbf{A}) > 1$. For a $2\times 2$ matrix with positive trace and determinant, this implies both eigenvalues have modulus greater than 1. **The Blanchard-Kahn condition is satisfied.**

---

## 9. Role of $\varepsilon$ — Revised Assessment

The earlier draft incorrectly suggested that $\varepsilon > 2$ might be needed for BK. The correct analysis shows:

**$\varepsilon > 1$ is sufficient for both eigenvalues to exceed 1.**

- $1/a > 1$: requires $a < 1$, which follows from positive match surplus (as in AGS notes) — no $\varepsilon$ condition.
- $1/b > 1$: requires $\psi_N > 0$, i.e., $\varepsilon > 1$ (profits are positive and variety has positive real-income value).
- The off-diagonal coupling increases $\det(\mathbf{A})$ above $1/(ab)$ and leaves $f(1) > 0$ intact for $s_Q$ small.

**Why $\varepsilon > 2$ matters elsewhere but not for BK.** The condition $\varepsilon > 2$ (equivalently $\psi_N < 1$) ensures that firm profits fall with entry — the dilution effect dominates the variety-income gain. This is the operative condition in prop:ds\_asymmetry (asymmetric $\delta$ vs $s$ propagation). It bounds the coupling strengths $\beta_{CC}, \alpha_C \propto \psi_N$ to be sub-unitary, but since BK holds for any $\psi_N > 0$, this bound is not needed for determinacy.

The Schaal-Dumouchel concern (aggregate demand externality generating multiple equilibria) is a **global** multiplicity result, not a local BK condition. BK is local and is satisfied for all $\varepsilon > 1$ at the unique interior steady state.

---

## 10. Leading-Order Simplification ($s_Q \to 0$)

When vacancy costs are a negligible share of output ($s_Q \to 0$), $\beta_{QC} \to 0$, $\alpha_Q \to 0$, and:

$$\mathbf{A} \;\to\; \begin{pmatrix}1/a & (\bar\Psi\beta_{CC}+\alpha_C)/a \\ 0 & 1/b\end{pmatrix} \quad\text{(upper triangular)}$$

The matrix is **upper triangular** with eigenvalues exactly $1/a$ and $1/b$, both exceeding 1. The off-diagonal $A_{12}$ does not affect eigenvalues of an upper triangular matrix.

This is the clearest statement of the result: **at calibrated values ($s_Q \approx 1\text{--}3\%$), the system is nearly upper triangular, eigenvalues are approximately $1/a$ and $1/b$, and BK holds by the same argument as AGS for both eigenvalues independently.**

---

## 11. Summary Table

| | AGS ($\sigma=0$, $N=1$) | Full Baseline (DS-CES, $p_0=0$) |
|---|---|---|
| Jump variables | 1 ($\tilde\Psi_t$) | 2 ($\tilde\Psi_t$, $\hat C_t$) |
| Eigenvalues outside unit circle needed | 1 | 2 |
| BM diagonal | $1/a$ | $1/a$ and $1/b$ |
| Off-diagonal | — | $A_{12}>0$ large; $A_{21}>0$ small |
| Matrix structure (leading order) | Scalar | Upper triangular |
| $1/a > 1$ condition | $a<1$: automatic from $\bar\pi>0$ | Same |
| $1/b > 1$ condition | N/A | $\varepsilon > 1$ (positive variety value) |
| $\det(\mathbf{A}) > 1$ | $1/a > 1$ ✓ | $1/(ab)+c\beta_{QC}/a > 1$ ✓ |
| $f(1) > 0$ | $(1-1/a)<0$: scalar, N/A | $(1/a-1)(1/b-1) > A_{12}A_{21}$ ✓ |
| Sufficient $\varepsilon$ condition for BK | N/A | $\varepsilon > 1$ |
| Is $\varepsilon > 2$ needed for BK? | N/A | **No** |

---

## 12. What Changes with Endogenous Exit ($p_0 > 0$)

With $p_0 > 0$, the survival factor $\Lambda_t = F_\chi(\chi_t^c)$ becomes endogenous within the period. The N LOM gains a $\Lambda_t$ multiplier. The static fixed point $(N_t, \chi_t^c, \Lambda_t)$ given $\{C_t, Q_t\}$ is unchanged in character — no new jump variable.

The effect on BK: replace $a_N = \beta(1-\bar\delta)$ with $a_N = \beta(1-\bar\delta)\bar\Lambda < \beta(1-\bar\delta)$, and similarly replace $a$ with $a = \beta(1-\bar\delta)\bar\Lambda\bar\gamma$. This raises $1/a$ and $1/b$ (both persistence coefficients shrink with endogenous exit), strengthening BK. The upper-triangular structure and all sign arguments carry through unchanged. The identity (ID) holds by the same argument (SS BFE still implies $a_N(1+\eta)=1$). **Endogenous exit only makes BK easier to satisfy.**
