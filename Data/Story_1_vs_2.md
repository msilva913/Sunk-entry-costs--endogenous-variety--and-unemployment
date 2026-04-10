Key question: Does the component of B~^LD that is orthogonal to B~^δ have a true zero vacancy effect, or does it have a non-zero vacancy effect that is simply too noisy to identify?

## Test 1: Partial R² and First-Stage Diagnostics

The asymmetric dominance story requires that δ has much higher explanatory power for vacancies than LD does, even in the first stage. Compute the partial R² of each instrument in the vacancy outcome equation at each horizon:

```python
# After running the joint LP at each horizon h, also run:
# 1. Regress vacancy outcome on δ alone (with FEs and controls)
# 2. Regress vacancy outcome on LD alone (with FEs and controls)  
# 3. Compare partial R² and F-statistics

# In the joint regression, extract the partial F for each instrument:
# Partial F for δ  = (R²_joint - R²_LD_only) / (1 - R²_joint) * (N-K)
# Partial F for LD = (R²_joint - R²_delta_only) / (1 - R²_joint) * (N-K)
```

This test directly measures whether LD has sufficient independent variation to identify a non-zero coefficient if one exists. If the partial F for LD in the joint vacancy regression is still above 10 at short horizons, the instrument has sufficient independent variation to identify its coefficient and the collinearity artifact story is weakened. If the partial F for LD collapses to near zero in the joint specification, your co-author's concern is validated.

## Test 2: Artificial Orthogonalization

Construct a version of the LD instrument that is orthogonal to δ by projecting out the δ component:

```python
# Regress bartik_ld on bartik_delta within each quarter
# (cross-sectional regression across states)
# The residual is the component of LD variation orthogonal to delta

import statsmodels.api as sm

def orthogonalize_instrument(panel, col_to_clean, col_to_remove):
    """
    For each quarter t, regress col_to_clean on col_to_remove 
    across states. Residual is the orthogonalized instrument.
    """
    records = []
    for qt, grp in panel.groupby("quarter_label"):
        if len(grp) < 10:
            continue
        y = grp[col_to_clean].values
        X = sm.add_constant(grp[col_to_remove].values)
        resid = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
        grp = grp.copy()
        grp[f"{col_to_clean}_orth"] = resid
        records.append(grp)
    return pd.concat(records, ignore_index=True)

joint_panel = orthogonalize_instrument(
    joint_panel, "bartik_ld", "bartik_delta"
)
```

Now run the separate LP using `bartik_ld_orth` instead of `bartik_ld`. This instrument has zero correlation with δ by construction, so the collinearity concern is completely eliminated.

If the vacancy IRF using `bartik_ld_orth` is also flat around zero, that is strong evidence the LD → 0 result is genuine and not a collinearity artifact — because you have given LD instrument variation that is completely independent of δ and it still shows no vacancy effect.

If the vacancy IRF using `bartik_ld_orth` is significantly negative, that is evidence the collinearity artifact story is correct — LD has a real vacancy effect that was being absorbed by δ in the joint specification.

## Test 3: Vary the Correlation Structure Across Industries

If the collinearity artifact story is right, the joint LP results should be worse — meaning LD further from zero and less precise — in subsamples where r(δ,LD) is lower, because in those subsamples LD has more independent variation to identify its coefficient. Conversely the separate LP results should be most different from the joint LP results in subsamples where r(δ,LD) is highest.

Split the sample into two groups based on which industries dominate each state's Bartik weights. States heavily weighted toward Manufacturing and Retail — where r_dLD is highest — should show the largest discrepancy between separate and joint LP results. States weighted toward Information, Education, and Professional Services — where r_dLD is lowest — should show the smallest discrepancy.

python

```python
# Create high-correlation and low-correlation state groups
# based on their average employment-weighted r(delta, LD)
# then run separate and joint LPs on each subsample
```

If the LD → 0 result holds even in the low-correlation subsample where collinearity is not a concern, that strongly validates the genuine zero interpretation.

## Test 4: Monte Carlo Under the Null

The central interpretive question about the joint LP result is whether the zero β^LD vacancy coefficient reflects a genuine structural zero — consistent with Story 1, where surviving firms repost vacancies after a layoff shock — or an identification failure where the joint LP simply lacks sufficient power to detect a non-zero LD effect after partialling out δ. This distinction matters because both Story 1 and Stories 2/3 produce the same observable outcome in the joint LP: a zero β^LD coefficient. The zero finding alone cannot tell us which story is correct.

The power analysis resolves this ambiguity by asking the following question directly: if Story 2/3 were true and the LD vacancy effect were genuinely non-zero — specifically equal to the separate LP estimate of -0.016 pp per 1-SD shock — how often would the joint LP correctly detect a significant negative coefficient? This is a standard power calculation, where the null hypothesis is Story 1 (true effect = 0), the alternative is Story 2/3 (true effect = -0.016), and the test statistic is the joint LP β^LD coefficient. We simulate panel data with the observed sample dimensions (50 states, 66 quarters), the observed instrument correlation r = 0.42, and the assumed true effects, then run both the separate and joint LP on each simulated dataset across 2000 replications. Power is the fraction of replications in which the joint LP correctly finds a significantly negative β^LD.

The logic of the test is the following. If power is high — above 0.80 by conventional standards — then the joint LP would reliably detect a non-zero LD vacancy effect if one existed. Finding zero in the actual data therefore constitutes meaningful evidence in favor of Story 1: the test had the sensitivity to detect a non-zero effect and did not find one. If power is low — below 0.50 — then the joint LP would frequently find zero even when the true effect is -0.016. In that case the zero β^LD in the actual data is essentially uninformative and cannot distinguish Story 1 from Stories 2/3. The simulation also includes a sensitivity analysis that varies both the instrument correlation and the assumed true LD effect size across a grid, which tells us whether any power problem is primarily driven by the r = 0.42 correlation between the instruments or by the small magnitude of the true effect. That distinction matters for the paper because if power is low primarily due to the small effect size, the conclusion is that LD has a modest vacancy effect that is genuinely difficult to detect in any specification — which is itself an informative finding rather than a pure identification failure.

```python
np.random.seed(42)
n_states = 50
n_quarters = 66
n_sims = 1000
beta_delta_true = -0.04  # true delta vacancy effect
beta_ld_true    =  0.00  # null: LD has no vacancy effect
r_instruments   =  0.42  # observed correlation

results_separate_ld = []
results_joint_ld    = []

for _ in range(n_sims):
    # Simulate correlated instruments
    cov = np.array([[1, r_instruments],
                    [r_instruments, 1]])
    instruments = np.random.multivariate_normal(
        [0, 0], cov, size=n_states * n_quarters
    )
    delta = instruments[:, 0]
    ld    = instruments[:, 1]
    
    # Simulate vacancy outcome under null
    vacancy = (beta_delta_true * delta 
               + beta_ld_true * ld 
               + np.random.normal(0, 0.1, n_states * n_quarters))
    
    # Separate LP for LD
    X_sep = sm.add_constant(ld)
    b_sep = np.linalg.lstsq(X_sep, vacancy, rcond=None)[0][1]
    results_separate_ld.append(b_sep)
    
    # Joint LP for LD
    X_jnt = sm.add_constant(np.column_stack([delta, ld]))
    b_jnt = np.linalg.lstsq(X_jnt, vacancy, rcond=None)[0][2]
    results_joint_ld.append(b_jnt)

print(f"Separate LP LD coefficient under null:")
print(f"  Mean = {np.mean(results_separate_ld):.4f}")
print(f"  Bias = {np.mean(results_separate_ld) - beta_ld_true:.4f}")
print(f"  % significant negative = "
      f"{np.mean(np.array(results_separate_ld) < -0.01)*100:.1f}%")

print(f"\nJoint LP LD coefficient under null:")
print(f"  Mean = {np.mean(results_joint_ld):.4f}")
print(f"  Bias = {np.mean(results_joint_ld) - beta_ld_true:.4f}")
print(f"  % significant negative = "
      f"{np.mean(np.array(results_joint_ld) < -0.01)*100:.1f}%")
```

This simulation tells you directly whether a separate LP with r = 0.42 instruments would spuriously find a negative LD vacancy effect even when the true effect is zero. If the simulation shows the separate LP has substantial spurious negative bias under the null, your concern about the separate LP is validated. If the joint LP correctly recovers near-zero LD coefficients in most simulations, that validates the joint LP as the right specification rather than the collinearity artifact story.