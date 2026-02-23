# Calibration Architecture Improvements

## Summary of Changes

The `calibrate_shares()` function has been refactored to dramatically improve **transparency** of the target-to-parameter mapping and **code maintainability**.

---

## Issues Fixed

### 1. **Opacity in Root-Finding Problems** ✓
**Before:** Loss functions were terse with minimal context about what equations they solve.
```julia
function loss_xc(x)
    π_s = 1 / ε - x
    Yc_YG = (r + δ_e) / (r + δ_e + δ_e * π_s)
    YG_Y = 1 + X_Y + Xc_Y
    return 100 * (Xc_Y / (Yc_YG * YG_Y) - x)
end 
```
(What is this solving for? Why this equation?)

**After:** Each root-finding section now includes:
- **Purpose comment**: What target does this solve?
- **Target-binding relationship**: The equation being solved
- **Variable definitions**: What each quantity represents
- **Bracket justification**: Why bounds are [0.01, 0.5] vs [1e-3, 1.0]

### 2. **Missing Mapping Between Targets and Parameters** ✓
**Before:** No explicit documentation of "Target_X pins Parameter_Y via Equation_Z"

**After:** Each stage includes a clear box:
```
Target: Xc_Y (fixed cost share)
Relationship: Xc_Y = (Xc_Yc) * (Yc_YG) * (YG_Y)
Solving: For Xc_Yc, which pins profit share π_s
```

### 3. **Unclear Ordering and Dependencies** ✓
**Before:** Four root-finding steps appeared in sequence with no indication of whether order matters or how they interact.

**After:** Explicit **5-stage pipeline** with numbered sections:
1. **Direct Conversions**: Annual→monthly, targets→fundamentals (no iteration)
2. **Profit Share**: Root-find (uses: Stage 1 outputs)
3. **Distribution Params**: Root-find (uses: Stage 2 outputs)
4. **Vacancy Value**: Root-find (uses: Stage 2–3 outputs)
5. **Final Calculations**: Extract parameters (uses: all previous stages)

### 4. **Variable Naming and Clarity** ✓
**Before:**
```julia
f = f / (1 - δ_e)  # Redefines input; confusing
q = q / (1 - δ_e)  # Same issue
```

**After:**
```julia
f_corr = f / (1 - δ_e)  # Explicit: "corrected rate"
q_corr = q / (1 - δ_e)  # Clear intent
```

### 5. **Hard-Coded Root-Finding Bounds** ✓
**Before:**
```julia
find_zero(loss_xc, [0.01, 0.5])    # Why these bounds?
find_zero(loss_psi, 0.1)           # Single point; risky
fzero(x -> loss_Q(x)[1], 0.1)      # Arbitrary init
```

**After:**
```julia
find_zero(loss_xc, (0.01, 0.5))   # Bracket: Xc_Yc ∈ (0,1); practical bounds explained
find_zero(loss_psi, (1e-3, 0.95)) # Robust bounds: cons ∈ (0,1)
find_zero(loss_Q_residual, (1e-3, 1.0))  # Conservative bracketing for vacancy value
```

### 6. **Redundant Function Evaluation** ✓
**Before:**
```julia
Q = fzero(x -> loss_Q(x)[1], 0.1)
Q = abs(Q)
out = loss_Q(Q)[2]  # Called AGAIN
```

**After:**
```julia
loss_Q_residual(Q) = loss_Q(Q)[1]
Q_opt = find_zero(loss_Q_residual, (1e-3, 1.0))
Q = abs(Q_opt)
_, Q_results = loss_Q(Q)  # Single subsequent call for extraction
```

### 7. **Intermediate Variables Extraction** ✓
**Before:** Variables computed inside loss functions with unclear role:
```julia
@unpack w_int, κ, z, f_e, K, d_f, ν_f, x_c, X_c, X, C, Y_c, Y = out
# What do these mean? Why are they grouped together?
```

**After:** Explicitly named tuple return with docs on each step:
```julia
# Step A–I: Each intermediate step clearly commented
return loss_val, (; w_int, κ, z, f_e, K, d_f, ν_f, x_c, X_c, X, C, Y_c, Y=Y_new)
```

### 8. **Wage Surplus Decomposition** ✓
**Before:**
```julia
ϕ = (w - b) / (w_int - K + θ * (K + q * κ) - b)
```

**After:**
```julia
w_surplus = w_int - K + θ * (K + q_corr * κ) - b  # Named intermediate
ϕ = (w - b) / w_surplus
```

---

## Information Architecture Improvements

### **Before: Dense Computation**
```
15 targets → direct conversions → 3 nested root-finds → 18 parameters
[Black box]
```

### **After: Transparent Pipeline**
```
15 targets
    ↓
STAGE 1: Direct Conversions (8 derivations)
    → δ_e, δ, r, μ, τ, s, z, ρ
    ↓
STAGE 2: Root-find Xc_Yc → π_s, L_c, L_e
    ↓
STAGE 3: Root-find cons → ψ, f_m
    ↓
STAGE 4: Root-find Q → {K, κ, w_int, z, f_e, ν_f, x_c, X}
    ↓
STAGE 5: Final calcs → ϕ, x_m, dest_el
    ↓
18 output parameters
```

---

## Docstring Enhancement

Added comprehensive docstring with three sections:

1. **Calibration Structure**: Explains 5 stages and which targets pin which parameters
2. **Stage-by-stage mapping**: 
   - Direct conversions list equations
   - Root-finding stages explain target → relationship → pin
3. **Output description**: "Returns a NamedTuple of 18 calibrated parameters for ParaCalib"

---

## Trade-offs and Considerations

| Aspect | Before | After | Impact |
|--------|--------|-------|--------|
| **Code length** | ~50 lines | ~170 lines | More verbose, but clear |
| **Readability** | Terse | Explicit | Easier onboarding |
| **Maintainability** | Hard to modify stages | Easy to audit/modify each stage | Better for future work |
| **Performance** | Slightly faster (fewer comments) | Negligible overhead | Comments don't run |
| **Debuggability** | Black box | Transparent pipeline | Easier to identify which stage fails |

---

## Recommended Next Steps

1. **Validate calibration**: Run `calibrate_shares(TARGETS)` and verify output parameters satisfy the steady-state equation
2. **Add sensitivity analysis**: Plot how small target changes affect output parameters (which are robust? which are sensitive?)
3. **Test edge cases**: What happens if `Xc_Y` is very small or very large? Document numerical stability
4. **Cache common quantities**: If `calibrate_shares` is called repeatedly, pre-compute `μ`, `τ`, etc.
5. **Add convergence diagnostics**: Print solver iteration counts and residuals at each stage for debugging

---

## Files Modified

- `Programs baseline/steady_state_refactored.jl` — `calibrate_shares()` function (~170 lines, fully documented)
