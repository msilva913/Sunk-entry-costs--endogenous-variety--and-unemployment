include("steady_state.jl")
using Printf

############################################################
# --- Calibration Targets and Steady-State Calculation -----
############################################################

# Baseline targets (current calibration)
# Uses PATH B (dest_elast_target=5.0) — Xc_Y is an outcome, not a target
targets = (
    X_Y=0.015,              # recruiting cost share of output
    dest_ann=0.0754,         # annual product destruction rate
    dest_end_frac=0.5,       # endogenous share of destruction rate
    dest_elast_target=5.0,   # destruction elasticity wrt x_c (PATH B)
    p_0=0.5,                 # probability of drawing from continuous part of F
    f=0.41,                  # job-finding rate
    η_L=0.6,                 # elasticity of matching function wrt unemployment
    q=0.8,                   # vacancy filling rate
    sep=0.031,               # aggregate separation rate
    b_ratio=0.9,             # ratio of unemployment benefits to wage
    x_v=0.5,                 # sunk share of total recruiting cost
    ξ_inv=1.0,               # entry elasticity inverse
    r_ann=0.04,              # annual discount rate
    ε=4.3,                   # elasticity of substitution
    σ=1.0,                   # inverse IES
    N=1.0,                   # SS mass of firms (normalization)
    w=1.0                    # SS wage (normalization)
)

# Calibrate and compute steady state
println("=== BASELINE (ξ_inv = 1.0) ===")
cal = calibrate_shares(targets)
steady = steady_state(cal)

############################################################
# ----------- Unpack Calibrated Variables ------------------
############################################################

@unpack θ, δ_e, x_c, ρ, L_c, L_e, w, w_int, L, N, N_e, K, Q, q, f, ν_f, d_f, C, Y_c, Y, X_v, X, X_c,
         labor_share, sunk_vac_cost_share, vacancy_share, x_v, M, e, v, u, π_s,
         profit_share_rec, profit_share_ret, dest_el = steady

@unpack f_e, δ, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, s, f_m, ψ, p_0 = cal

# Composites
ψ_c = ψ/(ψ+1)
ζ = (x_c/f_m)^ψ
τ = 1-(1-δ_e)*(1-s)
μ = ε/(ε-1)


############################################################
# -------------------- Consistency Checks ------------------
############################################################

# Separation elasticity ψ F(x_c)/(1-F(x_c))
dest_elast_1 = dest_elast(cal, x_c)
@show surv_prob = 1-p_0 + p_0*ζ
@assert abs(dest_elast_1 - dest_el) < 1e-12

# Normalizations and target matches
@assert abs(steady.N - targets.N) < 1e-12
@assert abs(steady.f*(1-δ_e) - targets.f) < 1e-12
@assert abs(steady.ρ - 1.0) < 1e-12
@assert abs(targets.X_Y - vacancy_share) < 1e-12

# NOTE: Under PATH B (dest_elast_target), Xc_Y is an outcome, not a target.
# Assertion below is PATH A only. Under PATH B we just report it.
println("  Implied Xc_Y = $(round(X_c/Y*100, digits=2))%  (PATH B: Xc_Y is an outcome, not a target)")

labor_share_alt = 1.0 - profit_share_rec - profit_share_ret
@assert abs(labor_share - labor_share_alt) < 1e-12

# Retail profits
@assert abs(N*d_f - π_s*Y_c) < 1e-12
@assert abs(N*d_f - Y_c/ε + X_c) < 1e-12

# Output identity
@assert abs(Y-C-ν_f*N_e) < 1e-12

# Production and cost identities
@assert abs(Y_c - ρ*z*L_c) < 1e-12
@assert abs(Y_c - C - X - X_c) < 1e-12

@assert abs(ρ*z*L_e/μ-ν_f*N_e) < 1e-12

# Vacancy cost share
@assert abs(κ/(κ+K/q) - (1-x_v)) < 1e-12

# Wage equation and job creation condition
@assert abs(w_int - w - K - (r+τ)/(1-δ_e)*(κ+K/q)) < 1e-12
@assert abs(w - ϕ*(w_int - K + θ*(K + q*κ)) - (1-ϕ)*b) < 1e-12

# Entrant consistency
@assert abs(v - ((1-δ_e)*((1-q)*v + s*(1-u)) + e)) < 1e-12
@assert abs(e - δ_e*(v + 1 - u)) < 1e-12
@assert abs(N_e - L_e*z/f_e) < 1e-12

# Surplus ratio and capital check
surplus_ratio = (r + τ)/(1 - δ_e) * (1/(q*x_v))
@assert abs(K - (1-ϕ)/ϕ*(w-b)/(surplus_ratio + θ*(1/x_v))) < 1e-12

# Consistency x_m with Q
@assert abs(Q - e^ξ_inv*x_m) < 1e-12

# Cutoff
@assert abs(x_c - ρ*f_e/μ*((r+δ_e)/((1-δ_e)*(μ*π_s))*(μ-1) + 1)) < 1e-12

############################################################
# ----------- Display Key Steady-State Ratios --------------
############################################################

@show N_jcc(steady.θ, δ_e,  cal)
@show N_res(steady.θ, δ_e,  cal)
@show labor_share
@show C/Y
@show X/Y
@show sunk_vac_cost_share
@show vacancy_share
@show X/(q*v*w)
@show x_v

# Profits
@show ((w_int-w)*L-X)/Y
@show (N*d_f)/Y

# Value of unemployment benefit
@show b/w_int
@show b/(Y/L)
@show b/z

@show M/(12*Y)

############################################################
# --- Stage-by-Stage Validation ----------------------------
############################################################

println("\n=== STAGE 1: Direct Conversions ===")
δ_e_expected = 1 - (1 - targets.dest_ann)^(1/12)
@assert abs(δ_e - δ_e_expected) < 1e-12 "Stage 1: δ_e ($(δ_e)) does not match conversion ($(δ_e_expected))"
δ_expected = (1 - targets.dest_end_frac) * δ_e_expected
@assert abs(cal.δ - δ_expected) < 1e-12 "Stage 1: δ ($(cal.δ)) does not match ($(δ_expected))"
r_expected = (1 + targets.r_ann)^(1/12) - 1
@assert abs(cal.r - r_expected) < 1e-12 "Stage 1: r ($(cal.r)) does not match ($(r_expected))"
τ_expected = targets.sep
τ_computed = τ
@assert abs(τ_computed - τ_expected) < 1e-12 "Stage 1: τ ($(τ_computed)) does not match target ($(τ_expected))"
println("✓ Stage 1 conversions verified")

println("\n=== STAGE 2: Profit Share (PATH B — dest_elast_target) ===")
# Under PATH B, π_s is derived analytically; ψ × ζ/(1−ζ) = dest_elast_target exactly
@assert abs(dest_el - targets.dest_elast_target) < 1e-10 "Stage 2: dest_el ($(dest_el)) ≠ target ($(targets.dest_elast_target))"
Xc_Yc_recovered = steady.X_c / steady.Y_c
π_s_accounting = 1 / targets.ε - Xc_Yc_recovered
@assert abs(π_s_accounting - π_s) < 1e-12 "Stage 2: π_s accounting ($(π_s_accounting)) ≠ π_s ($(π_s))"
println("✓ Stage 2: dest_elast target matched, π_s consistency verified")

println("\n=== STAGE 3: Distribution Parameters ===")
cons_recovered = ψ_c * cal.p_0
μ_check = cal.ε / (cal.ε - 1)
π_s_equilibrium = (μ_check - 1) / μ_check * (1 - cons_recovered) * (cal.r + δ_e) / (cal.r + δ_e + cons_recovered * (1 - δ_e))
@assert abs(π_s_equilibrium - π_s) < 1e-12 "Stage 3: π_s from equilibrium ($(π_s_equilibrium)) does not match π_s ($(π_s))"
println("✓ Stage 3: Distribution parameters consistent with π_s")

println("\n=== STAGE 4: Vacancy Value (Root-find Q) ===")
@assert abs(steady.X / steady.Y - targets.X_Y) < 1e-12
@assert abs(κ / (κ + steady.K / steady.q) - (1 - targets.x_v)) < 1e-12
println("✓ Stage 4: Recruiting share and vacancy split verified")

println("\n✓ All checks passed for baseline calibration\n")


############################################################
# =========================================================
# FREE-ENTRY COMPARISON: ξ_inv = 0  (ξ → ∞)
# =========================================================
# At free entry Q = x_m (constant), K = (r+δ_e)/(1+r) × Q exactly.
# Identification swap vs. baseline:
#   Dropped:    x_v  (no sunk/flow split at free entry)
#   New target: b/w_int = baseline value (pins w_int → z analytically)
#   Inherited:  κ = baseline value (not recalibrated)
#   Residual:   ϕ recovered from Nash (same as baseline Stage 5)
# Uses dedicated calibrate_shares_free_entry(targets, κ, b/w_int).
############################################################

println("=== FREE ENTRY (ξ_inv = 0) ===")
# κ inherited from baseline; b/w_int targeted at baseline value; x_v dropped; ϕ residual
cal_fe    = calibrate_shares_free_entry(targets, cal.κ, cal.b / steady.w_int)

# ── Diagnostic: evaluate SS residuals at baseline (θ, x_c) under cal_fe ─────
# If (θ_base, x_c_base) is a root of the free-entry system, both residuals
# should be near zero. If not, the free-entry equilibrium is genuinely different.
println("\n── Diagnostic: SS residuals at baseline (θ, x_c) under free-entry parameters ──")
let
    @unpack f_e, δ, s, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, f_m, ψ, p_0 = cal_fe
    μ    = ε / (ε - 1)
    cons = (ψ / (1 + ψ)) * p_0

    θ_test   = steady.θ
    x_c_test = steady.x_c

    δ_e_t = δ_e_fun(x_c_test, cal_fe)
    τ_t   = compute_separation_rate(δ_e_t, s)
    π_s_t = (μ-1)/μ * (1-cons) * (r+δ_e_t) / (r+δ_e_t + cons*(1-δ_e_t))
    f_t   = jf(θ_test, A, η_L)
    q_t   = vf(θ_test, A, η_L)
    K_t   = K_fun(θ_test, δ_e_t, cal_fe)
    L_t   = L_fun(θ_test, δ_e_t, cal_fe)

    lhs_jcc = (κ + K_t/q_t) * (r + τ_t + (1-δ_e_t)*ϕ*q_t*θ_test)
    N_t     = π_s_t * z * L_t * (1-δ_e_t) / (f_e * ((r+δ_e_t)/μ + δ_e_t*π_s_t))
    ρ_t     = N_t^(1/(ε-1))
    w_int_t = ρ_t * z / μ
    rhs_jcc = (1-δ_e_t) * (1-ϕ) * (w_int_t - K_t - b)
    L_c_t   = (r+δ_e_t) * L_t / (r+δ_e_t + δ_e_t*π_s_t*μ)
    Y_c_t   = ρ_t * z * L_c_t
    x_c_new = Y_c_t / (ε*N_t) + ρ_t*f_e/μ

    r1 = rhs_jcc - lhs_jcc
    r2 = log(x_c_new / x_c_test)

    @printf("  θ_test = %.6f,  x_c_test = %.6f\n", θ_test, x_c_test)
    @printf("  K_fe(θ_base) = %.6f  vs  K_base = %.6f   (x_m: %.4f → %.4f)\n",
            K_t, steady.K, cal.x_m, x_m)
    @printf("  N implied    = %.6f  (should be 1.0 if f_e is consistent)\n", N_t)
    @printf("  w_int implied = %.6f  (calibrated w_int = %.6f)\n", w_int_t, steady.w_int)
    @printf("  JCC residual r1 = %.6e  (lhs=%.6f, rhs=%.6f)\n", r1, lhs_jcc, rhs_jcc)
    @printf("  Exit residual r2 = %.6e  (x_c_new=%.6f, x_c=%.6f)\n", r2, x_c_new, x_c_test)
    println("  → If |r1|,|r2| ≈ 0: baseline (θ,x_c) is a free-entry root; solver init issue.")
    println("  → If |r1| large:    K change shifts JCC; free-entry equilibrium genuinely differs.")
end
println()

# ── JCC scan: map residual sign across theta to identify zero bracket ────────
# Uses internal jcc_res from steady_state_free_entry logic, replicated here.
println("\n── JCC scan across theta (cal_fe parameters) ──")
let
    @unpack f_e, δ, s, z, b, ϕ, r, ε, A, η_L, κ, x_m, ψ, p_0 = cal_fe
    μ    = ε / (ε - 1)
    cons = (ψ / (1 + ψ)) * p_0

    function vars_at_scan(θ_v, δ_e_v)
        τ     = compute_separation_rate(δ_e_v, s)
        π_s   = (μ-1)/μ * (1-cons) * (r+δ_e_v) / (r+δ_e_v + cons*(1-δ_e_v))
        K     = x_m * (r+δ_e_v) / (1+r)
        fv, qv = jf(θ_v, A, η_L), vf(θ_v, A, η_L)
        L     = L_fun(θ_v, δ_e_v, cal_fe)
        N     = π_s * z * L * (1-δ_e_v) / (f_e * ((r+δ_e_v)/μ + δ_e_v*π_s))
        ρ     = N^(1/(ε-1))
        L_c   = (r+δ_e_v) * L / (r+δ_e_v + δ_e_v*π_s*μ)
        x_c   = ρ * z * L_c / (ε*N) + ρ*f_e/μ
        w_int = ρ * z / μ
        return (; τ, K, f=fv, q=qv, w_int, x_c)
    end

    solve_de(θ_v) = find_zero(de -> δ_e_fun(vars_at_scan(θ_v, de).x_c, cal_fe) - de, (1e-6, 0.3))

    for θ_t in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.512, 0.6, 0.8, 1.0, 1.5, 1.8, 2.0]
        try
            δ_e_t = solve_de(θ_t)
            v = vars_at_scan(θ_t, δ_e_t)
            lhs = (κ + v.K/v.q) * (v.τ + r + (1-δ_e_t)*ϕ*v.q*θ_t)
            rhs = (1-δ_e_t) * (1-ϕ) * (v.w_int - v.K - b)
            @printf("  theta=%.3f  jcc_res=%+.4f  (lhs=%.4f, rhs=%.4f)\n", θ_t, rhs-lhs, lhs, rhs)
        catch err
            @printf("  theta=%.3f  FAILED: %s\n", θ_t, string(err)[1:min(60,length(string(err)))])
        end
    end
end
println()

steady_fe = steady_state_free_entry(cal_fe)

println("\n✓ Free-entry calibration and steady state computed\n")

# ── Extract free-entry steady state values (avoid name clash with baseline) ──
cal_fe_ϕ     = cal_fe.ϕ
cal_fe_κ     = cal_fe.κ
cal_fe_z     = cal_fe.z
cal_fe_b     = cal_fe.b
cal_fe_f_e   = cal_fe.f_e
cal_fe_x_m   = cal_fe.x_m
cal_fe_ψ     = cal_fe.ψ
cal_fe_f_m   = cal_fe.f_m
cal_fe_dest_el = cal_fe.dest_el

# ── Side-by-side table: calibrated parameters ────────────────────────────────

println("\n")
println("="^70)
println("  SIDE-BY-SIDE: CALIBRATED PARAMETERS")
println("  Baseline (ξ_inv=1)  vs.  Free Entry (ξ_inv→0)")
println("="^70)
@printf("  %-20s  %14s  %14s\n", "Parameter", "Baseline", "Free Entry")
println("-"^70)

# Parameters that differ
@printf("  %-20s  %14.6f  %14.6f\n", "ξ_inv",       cal.ξ_inv,    cal_fe.ξ_inv)
@printf("  %-20s  %14.6f  %14.6f\n", "ϕ",           cal.ϕ,        cal_fe_ϕ)
@printf("  %-20s  %14.6f  %14.6f\n", "κ",           cal.κ,        cal_fe_κ)
@printf("  %-20s  %14.6f  %14.6f\n", "z",           cal.z,        cal_fe_z)
@printf("  %-20s  %14.6f  %14.6f\n", "f_e",         cal.f_e,      cal_fe_f_e)
@printf("  %-20s  %14.6f  %14.6f\n", "x_m",         cal.x_m,      cal_fe_x_m)
@printf("  %-20s  %14.6f  %14.6f\n", "b",           cal.b,        cal_fe_b)
println("  -- (common parameters, unchanged) --")
@printf("  %-20s  %14.6f  %14.6f\n", "δ",           cal.δ,        cal_fe.δ)
@printf("  %-20s  %14.6f  %14.6f\n", "s",           cal.s,        cal_fe.s)
@printf("  %-20s  %14.6f  %14.6f\n", "r",           cal.r,        cal_fe.r)
@printf("  %-20s  %14.6f  %14.6f\n", "ε",           cal.ε,        cal_fe.ε)
@printf("  %-20s  %14.6f  %14.6f\n", "η_L",         cal.η_L,      cal_fe.η_L)
@printf("  %-20s  %14.6f  %14.6f\n", "ψ",           cal.ψ,        cal_fe_ψ)
@printf("  %-20s  %14.6f  %14.6f\n", "f_m",         cal.f_m,      cal_fe_f_m)
@printf("  %-20s  %14.6f  %14.6f\n", "dest_el",     cal.dest_el,  cal_fe_dest_el)
println("="^70)


# ── Side-by-side table: steady-state values ──────────────────────────────────

# Compute derived diagnostics for both
gross_surp_b = steady.w_int - cal.b
gross_surp_f = steady_fe.w_int - cal_fe_b
net_surp_b   = steady.w_int - steady.w - steady.K -
               (cal.r + targets.sep)/(1 - steady.δ_e)*(cal.κ + steady.K/steady.q)
net_surp_f   = steady_fe.w_int - steady_fe.w - steady_fe.K -
               (cal_fe.r + targets.sep)/(1 - steady_fe.δ_e)*(cal_fe_κ + steady_fe.K/steady_fe.q)
tight_rent_b = cal.ϕ * steady.θ * (steady.K / steady.x_v)
tight_rent_f = cal_fe_ϕ * steady_fe.θ * (steady_fe.K / steady_fe.x_v)

println("\n")
println("="^70)
println("  SIDE-BY-SIDE: STEADY-STATE VALUES")
println("  Baseline (ξ_inv=1)  vs.  Free Entry (ξ_inv→0)")
println("="^70)
@printf("  %-24s  %12s  %12s\n", "Variable", "Baseline", "Free Entry")
println("-"^70)

@printf("  %-24s  %12.6f  %12.6f\n", "θ (tightness)",        steady.θ,        steady_fe.θ)
@printf("  %-24s  %12.6f  %12.6f\n", "δ_e (monthly)",        steady.δ_e,      steady_fe.δ_e)
@printf("  %-24s  %12.6f  %12.6f\n", "u",                    steady.u,        steady_fe.u)
@printf("  %-24s  %12.6f  %12.6f\n", "v",                    steady.v,        steady_fe.v)
@printf("  %-24s  %12.6f  %12.6f\n", "N",                    steady.N,        steady_fe.N)
@printf("  %-24s  %12.6f  %12.6f\n", "w_int",                steady.w_int,    steady_fe.w_int)
@printf("  %-24s  %12.6f  %12.6f\n", "K (vacancy value)",    steady.K,        steady_fe.K)
@printf("  %-24s  %12.6f  %12.6f\n", "Q (entry cost)",       steady.Q,        steady_fe.Q)
@printf("  %-24s  %12.6f  %12.6f\n", "ϕ (calibrated)",       cal.ϕ,           cal_fe_ϕ)
println("-"^70)
@printf("  %-24s  %12.6f  %12.6f\n", "b/w_int",              cal.b/steady.w_int,   cal_fe_b/steady_fe.w_int)
@printf("  %-24s  %12.6f  %12.6f\n", "gross surplus (w_int-b)", gross_surp_b,   gross_surp_f)
@printf("  %-24s  %12.6f  %12.6f\n", "net surplus",          net_surp_b,      net_surp_f)
@printf("  %-24s  %12.6f  %12.6f\n", "tightness rent ϕθK/xv", tight_rent_b,   tight_rent_f)
@printf("  %-24s  %12.6f  %12.6f\n", "labor share",          steady.labor_share, steady_fe.labor_share)
@printf("  %-24s  %12.6f  %12.6f\n", "C/Y",                  steady.C/steady.Y,  steady_fe.C/steady_fe.Y)
@printf("  %-24s  %12.6f  %12.6f\n", "Xc_Y (implied)",       steady.X_c/steady.Y, steady_fe.X_c/steady_fe.Y)
println("="^70)


############################################################
# ---- Diagnostic: Free-entry assertions -------------------
############################################################
# With ξ_inv → 0:  Q = e^0 × x_m = x_m,  K = (r + δ_e)/(1+r) × x_m
# Check K is close to (r + δ_e_fe)/(1+r) × Q_fe (the free-entry value)
K_fe_pred = (cal_fe.r + steady_fe.δ_e) / (1 + cal_fe.r) * steady_fe.Q
println("\nFree-entry K check (should hold to machine precision at ξ_inv=0):")
@printf("  K (computed)            = %.8f\n", steady_fe.K)
@printf("  (r+δ_e)/(1+r) × Q      = %.8f\n", K_fe_pred)
@printf("  |difference|            = %.2e\n",  abs(steady_fe.K - K_fe_pred))

# b/w_int verification — should match baseline since ϕ is held fixed
@printf("\nb/w_int check:\n")
@printf("  baseline  b/w_int = %.6f\n", cal.b / steady.w_int)
@printf("  free entry b/w_int = %.6f\n", cal_fe.b / steady_fe.w_int)
println()
