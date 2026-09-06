# run_solution_core.jl
# =============================================================================
# Perturbation solution setup — Gabrovski-Silva endogenous exit model
# =============================================================================
# Based on toolkit by Alvaro Salazar-Perez and Hernán D. Seoane.
# Model: BGM (Bilbiie-Ghironi-Melitz 2012) + Coles-Kelishomi (2018) with
# endogenous exit via heterogeneous continuation costs.
#
# Three aggregate shocks: z (technology), δ (permanent exit), s (separation).
# Shocks are orthogonal by construction (Cholesky, z ordered first in part6b).
# The z→δ_e pathway is captured by equilibrium equations (x_c responds to z),
# not by shock covariance.
#
# Timing (monthly, following CK stage structure):
#   Stage 1: New realizations of (z_t, δ_t, s_t)
#   Stage 2: Bargaining and production — wages, output, profits
#   Stage 3: Vacancy investment — new entrants e_t drawn
#   Stage 4: Matching — m_t matches formed
#   Stage 5: Job separation at rate δ_e_t (endogenous) + s_t (idiosyncratic)
#
# State variables x = [u, N, v_pret, z, δ, s]:
#   u      : unemployment entering period t (pre-matching)
#   N      : mass of incumbent firms entering period t
#   v_pret : surviving vacancies from t-1 (pre-entry); total v_t = v_pret + e_t
#            entrants e_t post vacancies and participate in matching within period t
#            (draft eq:v_lom); v_pret makes vacancies partially predetermined.
#   z, δ, s: exogenous shocks (log-deviations from SS = 1)
#
# Control variables y (31 total) — see declaration below.
#
# Key endogenous exit objects added vs. CK baseline:
#   x_c  : continuation cost cutoff (firms with cost > x_c exit)
#   δ_e  : endogenous destruction rate = 1 - (1-δ*δbar)*F(x_c)
#          where F(x_c) = (1-p_0) + p_0*(x_c/f_m)^ψ
#   π_s  : per-firm profit share; determined by (ε, Xc_Y) calibration targets
# =============================================================================

using MKL
using DataFrames
using Serialization
cd(@__DIR__)

include("solution_functions.jl")
include("steady_state.jl")
include("impulse_response_plots.jl")
include("time_series_fun.jl")

function solution_interface(model, PAR, SS_precomputed::Union{Vector{Float64},Nothing}=nothing;
                            ss_tol::Float64=1e-8)
    # ── World-age fix ────────────────────────────────────────────────────────────
    # eval_* functions are generated inside process_model via eval(Meta.parse(...)).
    # When called from inside another function, Julia's world-age mechanism blocks
    # direct dispatch. Base.invokelatest always uses the latest compiled method.
    eta    = Base.invokelatest(eval_ShockVAR, PAR)
    PAR_SS = Base.invokelatest(eval_PAR_SS,  PAR)

    # ── Numeric SS ───────────────────────────────────────────────────────────────
    # Prefer a pre-computed Float64 SS when the caller supplies one (fast path,
    # avoids SymPy entirely).  Fall back to substituting model.SS symbolically
    # for callers (e.g. run_solution.jl) that don't pass a pre-computed SS.
    #
    # Why the fallback is needed: eval_SS has a bug in its generation loop
    # (`for ip in npar` iterates once instead of 1:npar), so it only substitutes
    # the last parameter; the other 22 stay as live SymPy globals. Passing that
    # Sym vector to eval_deriv causes a Float64 assignment crash.
    if !isnothing(SS_precomputed)
        SS = SS_precomputed
    else
        par_syms = [Sym("PAR[$i]") for i in eachindex(PAR)]
        function sym_to_float(v)
            for (p, val) in zip(model.parameters, PAR)
                v = subs(v, p, val)
            end
            for (ps, val) in zip(par_syms, PAR)
                v = subs(v, ps, val)
            end
            return Float64(v.evalf())
        end
        SS = [sym_to_float(model.SS[j]) for j in 1:length(model.SS)]
    end

    # ── Jacobian and residual (now with numeric SS) ──────────────────────────────
    SS_err = Base.invokelatest(eval_SS_error, PAR_SS, SS)
    deriv  = Base.invokelatest(eval_deriv,    PAR_SS, SS)
    SS_max = maximum(abs.(SS_err))
    println("Max SS residual: $SS_max")

    # ── Guard (added Sept 5, 2026) ────────────────────────────────────────────
    # Perturbation is only valid around an actual steady state. Before this
    # guard existed the residual was printed and ignored, and every mechanism
    # comparison ran with a residual of 12–22 on the free-entry condition f[18]
    # (SS_numeric derived its own f_e inconsistent with PAR's). Never downgrade
    # this to a warning: a first-order solution around a non-steady point still
    # returns a well-formed gx/hx and plausible-looking IRFs, so the failure is
    # silent by construction.
    if SS_max > ss_tol
        worst = argmax(abs.(SS_err))
        error("Steady-state residual $(SS_max) exceeds tolerance $(ss_tol); " *
              "worst equation is f[$(worst)]. The SS passed to the solver does " *
              "not satisfy the model equations — fix the SS before trusting any " *
              "IRF or moment computed from it.")
    end

    sol_mat = solve_model(model, deriv, eta)
    println("Model solved")

    # SS is Vector{Float64}; exp recovers levels from the log-deviation representation
    ss = NamedTuple(zip(model.varnames, exp.(SS[1:model.nvar])))
    return (SS=SS, ss=ss, eta=eta, deriv=deriv, sol_mat=sol_mat)
end

# =============================================================================
# Flags
# =============================================================================
flag_order     = 1
flag_deviation = true
flag_SSsolver  = false

# =============================================================================
# Parameters (symbolic)
# =============================================================================
# Order must match PAR vector in run_solution.jl.
# Removed: F (obsolete scalar CDF — endogenous F(x_c) is computed in equations)
# Added:   ψ, f_m, p_0 (continuation cost distribution parameters)
@syms f_e zbar δbar sbar b ϕ r σ ε A η_L κ ξ_inv x_m ψ f_m p_0 ρ_z σ_z ρ_δ σ_δ ρ_s σ_s

parameters = [f_e; zbar; δbar; sbar; b; ϕ; r; σ; ε; A; η_L; κ; ξ_inv; x_m; ψ; f_m; p_0;
              ρ_z; σ_z; ρ_δ; σ_δ; ρ_s; σ_s]

estimate = []   # filled in when SMM is wired up
position = []
priors   = (;)

# Composite parameters (symbolic)
β   = 1 / (1 + r)
ξ   = 1 / ξ_inv            # entry elasticity
μ   = ε / (ε - 1)          # markup

# =============================================================================
# Variables
# =============================================================================
# State variables (predetermined, known at start of period t)
@syms u N v_pret z δ s
@syms up Np v_pretp zp δp sp

# Control variables (jump variables, determined within period t)
# Added vs. previous version: x_c (exit threshold), δ_e (endogenous exit rate)
@syms θ q L v e K Q ρ N_e ν_f d_f w_int w L_e L_c Y_c C λ Y x_c δ_e labor_prod C_R Y_R Y_cR w_R ls
@syms θp qp Lp vp ep Kp Qp ρp N_ep ν_fp d_fp w_intp wp L_ep L_cp Y_cp Cp λp Yp x_cp δ_ep labor_prod_p C_Rp Y_Rp Y_cRp w_Rp lsp

x  = [u; N; v_pret; z; δ; s]
y  = [θ; q; L; v; e; K; Q; ρ; N_e; ν_f; d_f; w_int; w; L_e; L_c; Y_c; C; λ; Y; x_c; δ_e; labor_prod; C_R; Y_R; Y_cR; w_R; ls]
xp = [up; Np; v_pretp; zp; δp; sp]
yp = [θp; qp; Lp; vp; ep; Kp; Qp; ρp; N_ep; ν_fp; d_fp; w_intp; wp; L_ep; L_cp; Y_cp; Cp; λp; Yp; x_cp; δ_ep; labor_prod_p; C_Rp; Y_Rp; Y_cRp; w_Rp; lsp]

variables = [x; y; xp; yp]
varnames  = vcat(Symbol.(x), Symbol.(y))

nx   = length(x)
ny   = length(y)
nvar = nx + ny
ne   = 3    # three structural shocks: z, δ, s

# =============================================================================
# Shock matrix eta (nx × ne = 6 × 3)
# =============================================================================
# Column j of eta gives the impact of shock j on each state variable.
# Sign convention: a positive z shock raises log(z), so column 1 entry 4 = +σ_z.
# A positive δ shock raises log(δ), so column 2 entry 5 = +σ_δ.
# Shocks are orthogonal (Cholesky z-first in part6b); no off-diagonal terms.
#
# States: row 1=u, row 2=N, row 3=v_pret, row 4=z, row 5=δ, row 6=s
#         z shock   δ shock   s shock
eta = [0    0    0;    # u      (not directly hit)
       0    0    0;    # N      (not directly hit)
       0    0    0;    # v_pret (not directly hit)
       σ_z  0    0;    # z      (+σ_z on z shock)
       0    σ_δ  0;    # δ      (+σ_δ on δ shock)
       0    0    σ_s]  # s      (+σ_s on s shock)

@syms epsilon   # placeholder required by toolkit; ne=3 shocks handled via eta
ex = [epsilon]

# =============================================================================
# Model equations
# =============================================================================
# Notation:
#   Unprimed = period t (current); primed = period t+1 (next period, expected)
#   Shocks z, δ, s enter as multiplicative deviations from SS (= 1 at SS).
#   δ_e = endogenous destruction rate = 1 - (1 - δ*δbar)*F(x_c)
#     where F(x_c) = (1-p_0) + p_0*(x_c/f_m)^ψ   [continuation cost CDF]
#   x_c = continuation cost threshold (firms draw cost each period; exit if > x_c)
#
# Equation count: 6 state + 27 control = 33 total, matching nvar.
# f[1]–f[2]:   endogenous exit (x_c threshold, δ_e rate)
# f[3]–f[5]:   Euler equations (JCC, business formation, K value)
# f[6]–f[8]:   wage block (MRP, Nash wage, vacancy creation)
# f[9]–f[15]:  labor/goods market static conditions
# f[16]–f[24]: resource constraint, identities, LOMs
# f[25]:       labor share
# f[26]–f[30]: data-consistent observables
# f[31]–f[33]: exogenous AR(1) shock processes
#
# Key timing notes:
#   - δ_e,t applies to incumbents at END of period (Stage 5): destroys filled
#     jobs and unfilled vacancies simultaneously.
#   - v_pret is predetermined: surviving vacancies from t-1 before new entry.
#   - e_t new entrants post vacancies in Stage 3 and join the matching pool
#     in Stage 4 of the same period (draft eq:v_lom). Total v_t = v_pret + e_t.
#   - Matching uses total v_t: θ_t = v_t/u_t, f(θ_t) = A*θ_t^(1-η_L). ✓
#   - v_pret_{t+1} = (1-δ_e_t)*[(1-q_t)*v_t + s_t*sbar*(1-u_t)]: surviving
#     unmatched vacancies plus reposted separations, both from total v_t.

function gen_model_equations()
    f = fill(Sym("x"), nvar)

    # ── Endogenous exit composites ─────────────────────────────────────────────
    # Survival probability Λ(x_c) = F(x_c) evaluated at current and next cutoff.
    # Draft eq:Lambda: Λ_t = (1-p_0) + p_0*(x_c_t/f_m)^ψ
    Λ    = (1 - p_0) + p_0 * (x_c  / f_m)^ψ   # survival prob, period t
    Λp   = (1 - p_0) + p_0 * (x_cp / f_m)^ψ   # survival prob, period t+1

    # δ_e_t = 1 - (1-δ_t)*Λ_t  [draft eq:delta_e_lom, δ predetermined]
    # δ_e_{t+1} = 1 - (1-δ_{t+1})*Λ_{t+1}  [next period, not yet realized]
    #   (1-δ_t)*Λ_{t+1} = "fraction of today's firms+jobs surviving to t+1"
    # This uses CURRENT δ (predetermined) and NEXT-PERIOD Λp (depends on x_cp).
    SDF_surv = (1 - δ * δbar) * Λp   # (1-δ_t)*Λ_{t+1}: survival factor for t→t+1

    # ── [1] Exit threshold (eq:cutoff_eq) ─────────────────────────────────────
    # χ_t^c = Y_t^c*(μ-1)/(μ*N_t) + ρ(N_t)*f_e/μ(N_t)
    # = per-firm gross profit R^f + firm value ν_f  (eq:cutoff_free_entry)
    f[1] = x_c - (Y_c * (μ - 1) / (μ * N) + ν_f)

    # ── [2] Endogenous destruction rate (eq:delta_e_lom) ──────────────────────
    # δ_e_t = 1 - (1-δ_t)*Λ_t  where δ_t is predetermined state
    f[2] = δ_e - (1 - (1 - δ * δbar) * Λ)

    # ── [3] Job creation condition (eq:jcc_eq) ────────────────────────────────
    # κ + K/q = E_t[m_{t+1}*(1-δ_t)*Λ_{t+1}*{(1-ϕ)*(w_int'-K'-b)
    #            - ϕ*θ'*(K'+q'*κ) + (1-s'*sbar)*(κ+K'/q')}]
    # Survival factor: (1-δ_t)*Λ_{t+1} = SDF_surv (current δ, next-period Λp).
    f[3] = κ + K/q - β*λp/λ * SDF_surv * (
               (1-ϕ)*(w_intp - Kp - b)
               - ϕ*θp*(Kp + qp*κ)
               + (1 - sp*sbar)*(κ + Kp/qp)
           )

    # ── [4] Business formation Euler (eq:firm_value_char / BGM eq:euler_bf) ───
    # ν_f = β*(λ'/λ)*(1-δ_t)*Λ_{t+1}*(ν_f' + d_f')
    # Interpretation (BGM): marginal cost of creating a firm today (ν_f, via free
    # entry) equals discounted expected payoff tomorrow — dividend d_f' plus
    # continuation value ν_f', survival-weighted by SDF_surv = (1-δ_t)*Λ_{t+1}.
    # Together with free entry f[17] (ν_f = ρ*f_e/μ), this pins the entry margin
    # dynamically: it is not redundant — it is what makes N forward-looking.
    f[4] = ν_f - β*λp/λ * SDF_surv * (ν_fp + d_fp)

    # ── [5] Capital value K (eq:Kdef) ─────────────────────────────────────────
    # K_t = Q_t - (1-δ_t)*E_t[m_{t+1}*Λ_{t+1}*Q_{t+1}]
    # Q is the sunk cost of posting a vacancy; K is its net value after survival.
    f[5] = K - (Q - β*λp/λ * SDF_surv * Qp)

    # ── [6] Marginal revenue product (eq:recruiter_compensation) ──────────────
    # w_int = ρ(N)*z/μ(N) = ρ*z*zbar/μ  (DS-CES: ρ = N^(1/(ε-1)), μ constant)
    f[6] = w_int - ρ * z * zbar / μ

    # ── [7] Nash bargaining wage ───────────────────────────────────────────────
    # w = ϕ*(w_int - K + θ*(K + q*κ)) + (1-ϕ)*b
    f[7] = w - (ϕ*(w_int - K + θ*(K + q*κ)) + (1-ϕ)*b)

    # ── [8] Vacancy creation (eq:entry_eq) ────────────────────────────────────
    # Q_t = x_m * e_t^(ξ_inv): marginal cost of posting the e_t-th vacancy
    f[8] = Q - x_m * e^(ξ_inv)

    # ── Labor market ───────────────────────────────────────────────────────────
    # [9] Market tightness
    f[9] = θ - v/u

    # [10] Vacancy-filling probability (Cobb-Douglas matching)
    f[10] = q - A * θ^(-η_L)

    # [11] Total employment
    f[11] = L - (1 - u)

    # [12] Labor decomposition: production workers + recruiters
    f[12] = L - (L_c + L_e)

    # ── Goods market ───────────────────────────────────────────────────────────
    # [13] Relative price (DS-CES aggregator, N varieties): ρ = N^(1/(ε-1))
    f[13] = ρ - N^(1/(ε-1))

    # [14] Household Euler (λ = marginal utility of consumption)
    f[14] = λ - C^(-σ)

    # [15] Retail production function
    f[15] = Y_c - ρ * z * zbar * L_c

    # ── [16] Resource constraint (eq:rc): Y_c = C + X + X_c ──────────────────
    # X   = e/(1+ξ_inv)*Q + κ*q*v  (sunk entry costs + matching costs [Pissarides 2009])
    #   κ is paid per match (not per vacancy): total matching cost = κ * q(θ)*v
    # X_c = N*p_0*ψ_c*x_c          (aggregate continuation costs)
    # ψ_c = ψ/(ψ+1)
    ψ_c = ψ / (ψ + 1)
    X_c = N * p_0 * ψ_c * x_c
    # Total vacancy creation cost (draft eq. 29, p. 20):
    #   e/(1+ξ_inv)*Q  sunk posting costs — integral of the marginal schedule
    #                  Q = x_m*e^ξ_inv, i.e. ∫₀^e x_m u^ξ_inv du = e*Q/(1+ξ_inv)
    #   κ*q(θ)*v       fixed matching costs (κ per match, q*v matches)
    # FIXED Sept 5, 2026: was `e * ξ_inv/(1+ξ_inv) * Q`, a symbol-for-symbol
    # substitution of ξ→ξ_inv into the original ξ/(ξ+1) instead of an algebraic
    # inversion. Since ξ = 1/ξ_inv, ξ/(ξ+1) = 1/(1+ξ_inv). The two coincide at
    # ξ_inv = 1, so only Comparison D's ξ_inv = 0.1 arm was affected (10× error).
    # steady_state.jl:362,742 and SS_numeric (line 439) already had it right.
    X   = e / (1 + ξ_inv) * Q + κ * q * v
    f[16] = Y_c - (C + X + X_c)

    # [17] New entrants: N_e = z*zbar*L_e/f_e
    f[17] = N_e - z * zbar * L_e / f_e

    # ── [18] Free entry condition (eq:free_entry) ──────────────────────────────
    # ν_f = ρ(N)*f_e/μ(N) = ρ*f_e/μ  (DS-CES)
    f[18] = ν_f - ρ * f_e / μ

    # ── [19] GDP (eq:gdp): Y = C + ν_f*N_e ───────────────────────────────────
    f[19] = Y - (C + ν_f * N_e)

    # ── [20] Income identity: Y = w_int*L + N*d_f ─────────────────────────────
    # Retail profits net of continuation costs: d_f = Y_c/(ε*N) - X_c/N
    # Under DS-CES: Y_c/ε = R^f*N + markupˉ¹*Y_c ... equivalently:
    # d_f*N = Y_c*(μ-1)/μ - X_c = (Y_c - w_int*L_c) - X_c = π^f (net profits)
    # This is an identity given other equations; use as check/residual for d_f.
    f[20] = N * d_f - (Y_c * (μ - 1) / μ - X_c)

    # ── Laws of motion ─────────────────────────────────────────────────────────
    # [21] Total vacancies (eq:v_lom split): v = v_pret + e
    f[21] = v - (v_pret + e)

    # ── [22] Pre-entry vacancies tomorrow (predetermined part of v_lom) ────────
    # v_pret_{t+1} = (1-δ_e_t)*[(1-q_t)*v_t + s_t*sbar*(1-u_t)]
    # (draft eq:v_lom): surviving unmatched vacancies + reposted separations.
    # Uses total v_t (including e_t): entrants participate in matching this period.
    f[22] = v_pretp - (1 - δ_e)*((1 - q)*v + s*sbar*(1 - u))

    # ── [23] Unemployment LOM (eq:u_lom) ──────────────────────────────────────
    # u_{t+1} = [1 - (1-δ_{e,t})f(θ_t)]*u_t + τ_t*(1 - u_t)
    # where f(θ_t) = A*θ_t^(1-η_L) is the job-finding rate using total v_t
    # (entrants e_t participate in matching within period t), and
    # τ_t = δ_e + s*sbar*(1-δ_e)  [total separation rate, eq:tau_lom]
    #
    # CONVENTION A (pre-matching), restored September 5, 2026.
    # Matches formed at t reach t+1 only if the firm survives — hence the
    # (1-δ_e) on the outflow from unemployment — while τ_t applies to the
    # INHERITED stock (1-u_t) only. A match therefore cannot dissolve before it
    # has produced. This is the base the vacancy LOM f[22] reposts on, so the
    # pair conserves job positions exactly: v_pre' + L' = (1-δ_e)(v + L).
    #
    # This equation had been written with a post-matching base,
    #     up - ((1 - f_t)*u + τ_t*((1 - u) + f_t*u)),
    # which applies τ to (1-u)+f·u while f[22] reposts only on (1-u). That
    # mixture leaks positions at rate (1-δ_e)·s·q·v: workers matched at t who
    # separate at t leave employment without their job returning as a vacancy.
    # It is also inconsistent with the steady-state formulas eq:u_ss and
    # eq:e_ss, and with Gabrovski-Silva, which this model is meant to nest.
    # See Notes/LOM_timing_consistency.md.
    τ_t = δ_e + s*sbar*(1 - δ_e)
    f_t  = A * θ^(1 - η_L)         # job-finding rate: f(θ_t), θ_t uses total v_t
    f[23] = up - ((1 - (1 - δ_e)*f_t)*u + τ_t*(1 - u))

    # ── [24] Firm LOM (eq:N_lom_eq simplified) ────────────────────────────────
    # N_{t+1} = (1-δ_{e,t})*(N_t + N_e_t)
    # Both incumbents and entrants face current-period destruction δ_e. ✓
    f[24] = Np - (1 - δ_e)*(N + N_e)

    # [25] Labor share
    f[25] = ls - w*L/Y

    # ── Data-consistent observables ────────────────────────────────────────────
    # [26] Labor productivity (real output per worker)
    f[26] = labor_prod - Y/(ρ*L)

    # [27]–[30] Real (price-deflated) aggregates
    f[27] = C_R    - C/ρ
    f[28] = Y_R    - Y/ρ
    f[29] = Y_cR   - Y_c/ρ
    f[30] = w_R    - w/ρ

    # ── Exogenous AR(1) shock processes (eq:ar1_z – eq:ar1_s) ─────────────────
    # Log-linear around SS = 1 for each shock (log-deviation = 0 at SS).
    # Shock innovations are orthogonal by Cholesky construction (part6b).
    f[31] = log(zp) - ρ_z * log(z)   # technology
    f[32] = log(δp) - ρ_δ * log(δ)   # structural exit (⊥ z)
    f[33] = log(sp) - ρ_s * log(s)   # worker separation

    return f
end

# =============================================================================
# Symbolic steady state
# =============================================================================
# Provides the analytic SS to the perturbation toolkit, avoiding numerical
# root-finding. Uses calibrate_shares → steady_state pipeline from steady_state.jl.
#
# Key changes from previous version:
#   - Removed labor_share target (it is an outcome, not an input in new calibration)
#   - Uses dest_end_frac to split δ_e into exogenous/endogenous components
#   - x_c and δ_e now appear explicitly in the SS vector
#   - δ_e replaces all occurrences of δbar in LOM expressions
#   - p_0, ψ, f_m replace scalar F

# Numeric version of SS computation (used when PAR is available)
function SS_numeric(PAR::Vector{Float64}, targets)
    # =========================================================================
    # REWRITTEN September 5, 2026.
    #
    # This function used to derive the steady state analytically, from scratch,
    # in parallel with calibrate_shares/steady_state in steady_state.jl. Its
    # entry cost was
    #     f_e = zbar*(μ-1)*L/N*(1-δ_e)/(δ_e*μ + r)
    # which is the correct closed form ONLY when X_c = 0 (no stochastic
    # continuation costs, i.e. p_0 = 0). It uses the gross Lerner share
    # (μ-1)/μ = 1/ε where the model with endogenous exit requires the profit
    # share NET of continuation costs, π_s = 1/ε - X_c/Y_c. The correct form is
    # steady_state.jl:745,
    #     f_e = π_s*z*L_c*(1-δ_e)*μ/(N*(r+δ_e)).
    #
    # Consequence: f_e was overstated by the factor ((μ-1)/μ)/π_s ≈ 2.66, hence
    # ν_f = ρ*f_e/μ likewise, so the free-entry condition f[18] evaluated against
    # PAR's f_e left a residual of 12-22 depending on calibration. Every
    # mechanism comparison (A-D) was linearized around that non-steady point.
    # The error cascaded: inflated ν_f -> inflated x_c (f[1]) -> inflated
    # X_c = N*p_0*ψ_c*x_c -> π_s = 1/ε - X_c/Y_c turned NEGATIVE.
    #
    # The fix is structural, not arithmetic: having a second independent SS
    # derivation was the root cause. This is now a thin mapper over
    # steady_state(cal), the construction that steady_state_checks.jl certifies
    # against f[16], f[18], f[19], f[20], the JCC, the Nash wage and both LOMs
    # at 1e-12. DO NOT reintroduce a parallel derivation here.
    # =========================================================================

    # Rebuild the calibration from targets. calibrate_shares is deterministic,
    # so this reproduces exactly the `cal` the caller used to build PAR.
    cal = calibrate_shares(targets)

    # Guard: PAR must be the vector implied by this calibration. If a caller
    # ever builds PAR differently, fail here rather than silently linearizing
    # around a point inconsistent with the parameters handed to the solver —
    # that inconsistency is precisely the bug this rewrite removes.
    par_expected = [cal.f_e, cal.z, cal.δ, cal.s, cal.b, cal.ϕ, cal.r, cal.σ,
                    cal.ε, cal.A, cal.η_L, cal.κ, cal.ξ_inv, cal.x_m,
                    cal.ψ, cal.f_m, cal.p_0]
    par_names = ["f_e","zbar","δbar","sbar","b","ϕ","r","σ","ε","A","η_L",
                 "κ","ξ_inv","x_m","ψ","f_m","p_0"]
    for i in eachindex(par_expected)
        if !isapprox(PAR[i], par_expected[i]; rtol=1e-10, atol=1e-14)
            error("SS_numeric: PAR[$i] ($(par_names[i])) = $(PAR[i]) does not match " *
                  "calibrate_shares(targets) = $(par_expected[i]). PAR and the steady " *
                  "state would describe different economies.")
        end
    end

    st = steady_state(cal)
    @unpack u, N, v_pret, θ, q, L, v, e, K, Q, ρ, N_e, ν_f, d_f, w_int, w,
            L_e, L_c, Y_c, C, Y, x_c, δ_e, labor_share = st

    # Shocks are multiplicative deviations with SS = 1 (log-deviation = 0)
    z_s, δ_s, s_s = 1.0, 1.0, 1.0

    # Marginal utility (f[14]) and the data-consistent observables f[26]-f[30].
    # Recomputed here rather than taken from `st` so that this block is a literal
    # transcription of the model equations.
    λ_s          = C^(-cal.σ)
    labor_prod_s = Y / (ρ * L)      # f[26]
    C_R_s        = C / ρ            # f[27]
    Y_R_s        = Y / ρ            # f[28]
    Y_cR_s       = Y_c / ρ          # f[29]
    w_R_s        = w / ρ            # f[30]
    ls_s         = labor_share      # f[25] = w*L/Y

    # ── Assemble SS vector ────────────────────────────────────────────────────
    # Order MUST match varnames = vcat(Symbol.(x), Symbol.(y)):
    #   x = [u, N, v_pret, z, δ, s]
    #   y = [θ, q, L, v, e, K, Q, ρ, N_e, ν_f, d_f, w_int, w, L_e, L_c, Y_c,
    #        C, λ, Y, x_c, δ_e, labor_prod, C_R, Y_R, Y_cR, w_R, ls]
    SS_block = [log(val) for val in [
        u, N, v_pret, z_s, δ_s, s_s,
        θ, q, L, v, e, K, Q,
        ρ, N_e, ν_f, d_f,
        w_int, w, L_e, L_c, Y_c,
        C, λ_s, Y,
        x_c, δ_e,
        labor_prod_s, C_R_s, Y_R_s, Y_cR_s, w_R_s, ls_s
    ]]

    SS = vcat(SS_block, SS_block)   # current and future (same at SS)
    return SS
end

# ⚠️ CARRIES THE SAME DEFECT AS THE OLD SS_numeric (see the note there).
# Its f_e (below) uses the gross Lerner share (μ-1) in place of the net profit
# share π_s·μ, i.e. it solves the p_0 = 0 / X_c = 0 model. The defect is LATENT,
# not active: process_model uses this only to generate `eval_SS(PAR)`, which is
# bypassed whenever the caller passes SS_precomputed — and every current runner
# does. It cannot simply call steady_state(cal), which is numeric, so it is left
# in place. Do NOT use its output as a linearization point. If the symbolic SS
# path is ever revived, fix this first.
function SS_symbolics(parameters::Vector{Sym{PyObject}}, targets)

    f_e, zbar, δbar, sbar, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, ψ, f_m, p_0,
        ρ_z, σ_z, ρ_δ, σ_δ, ρ_s, σ_s = parameters

    @unpack N, w, f, q, X_Y, Xc_Y, dest_ann, dest_end_frac = targets

    # ── Monthly rates ─────────────────────────────────────────────────────────
    N_s    = N
    w_s    = w

    # Endogenous destruction rate (monthly) at SS: δ_e = 1-(1-dest_ann)^(1/12)
    δ_e_s   = 1 - (1 - dest_ann)^(1/12)

    # Exogenous destruction rate at SS: δ = (1-dest_end_frac)*δ_e
    δ_s_ss  = (1 - dest_end_frac) * δ_e_s

    # Effective job finding/filling (conditional on firm survival)
    f_corr  = f / (1 - δ_e_s)
    q_corr  = q / (1 - δ_e_s)

    # Aggregate separation rate and SS s
    sep_s   = targets.sep
    s_ss    = (sep_s - δ_e_s) / (1 - δ_e_s)

    # Matching and labor market
    θ_s     = f_corr / q_corr
    u_s     = sep_s / (sep_s + (1 - δ_e_s) * f_corr)
    L_s     = 1 - u_s
    v_s     = θ_s * u_s
    e_s     = δ_e_s * (v_s + 1 - u_s)
    v_prets = v_s - e_s

    # Firms and relative price ρ = N^(1/(ε-1))  [consistent with steady_state.jl]
    μ_s      = ε / (ε - 1)
    ρ_val    = N_s^(1/(ε-1))   # local name avoids collision with ρ_s (s-shock persistence)

    # ── x_c and δ_e at SS ─────────────────────────────────────────────────────
    # x_c_ss = Y_c*(μ-1)/(μ*N) + ν_f  (eq:cutoff_eq = f[1] at SS)
    # This is computed after ν_f_s and Y_c_s are known below.
    # δ_e_sym = δ_e_s (numeric calibration target)

    # ── Value functions ────────────────────────────────────────────────────────
    β_s      = 1 / (1 + r)
    τ_s      = sep_s
    surplus_ratio = (r + τ_s) / (1 - δ_e_s) * (1 / (q_corr * targets.x_v))
    K_s      = (1 - ϕ) / ϕ * (w_s - b) / (surplus_ratio + θ_s / targets.x_v)
    κ_s      = (1 - targets.x_v) / targets.x_v * K_s / q_corr
    Q_s      = K_s * (1 + r) / (r + δ_e_s)

    # w_int from JCC (interior labor productivity)
    w_int_s  = w_s + K_s + (r + τ_s) / (1 - δ_e_s) * (κ_s + K_s / q_corr)

    # zbar pinned by w_int = ρ*z*zbar/μ at SS (z=1)
    zbar_s   = μ_s * w_int_s / ρ_val

    # Entry cost and recruiters
    f_e_s    = zbar_s * (μ_s - 1) * L_s / N_s * (1 - δ_e_s) / (δ_e_s * μ_s + r)
    ν_f_s    = ρ_val * f_e_s / μ_s
    d_f_s    = (r + δ_e_s) / (1 - δ_e_s) * ν_f_s
    N_e_s    = δ_e_s / (1 - δ_e_s) * N_s
    L_e_s    = N_e_s * f_e_s / zbar_s
    L_c_s    = L_s - L_e_s
    Y_c_s    = ρ_val * zbar_s * L_c_s

    # x_c at SS: consistent with f[1] = x_c - (Y_c*(μ-1)/(μ*N) + ν_f)
    x_c_s    = Y_c_s * (μ_s - 1) / (μ_s * N_s) + ν_f_s
    δ_e_sym  = δ_e_s   # numeric calibration target

    # Continuation cost aggregate at SS (consistent with f[16])
    ψ_c_s    = ψ / (ψ + 1)
    X_c_s    = N_s * p_0 * ψ_c_s * x_c_s

    # x_m: pinned by entry condition f[6]: Q = x_m * e^(ξ_inv) at SS
    x_m_s    = Q_s / e_s^(ξ_inv)

    # Consumption and output (resource constraint f[16]: Y_c = C + X + X_c)
    # X = e/(1+ξ_inv)*Q + κ*q*v  (sunk entry costs + matching costs: κ per match, q*v matches)
    X_s      = e_s / (1 + ξ_inv) * Q_s + κ_s * q_corr * v_s
    C_s      = Y_c_s - X_s - X_c_s
    λ_s      = C_s^(-σ)
    Y_s      = C_s + ν_f_s * N_e_s   # GDP: eq:gdp consistent with f[19]

    # ── Shock SS values = 1 (log-linear AR(1) around log=0) ──────────────────
    z_s  = 1.0
    δ_s  = 1.0   # shock deviation (not the rate δbar; rate is δ_e_s)
    s_s  = 1.0

    # ── Data-consistent observables ───────────────────────────────────────────
    labor_prod_s = Y_s / (ρ_val * L_s)
    C_R_s        = C_s / ρ_val
    Y_R_s        = Y_s / ρ_val
    Y_cR_s       = Y_c_s / ρ_val
    w_R_s        = w_s / ρ_val
    ls_s         = w_s * L_s / Y_s

    # ── Assemble SS vector ─────────────────────────────────────────────────────
    # Order must match [x; y]: [u,N,v_pret,z,δ,s, θ,q,L,v,e,K,Q,ρ,N_e,ν_f,d_f,w_int,w,L_e,L_c,Y_c,C,λ,Y,x_c,δ_e,labor_prod,C_R,Y_R,Y_cR,w_R,ls]
    SS_block = [log(x) for x in [
        u_s, N_s, v_prets, z_s, δ_s, s_s,          # states
        θ_s, q_corr, L_s, v_s, e_s, K_s, Q_s,      # controls 1–7
        ρ_val, N_e_s, ν_f_s, d_f_s,                 # controls 8–11
        w_int_s, w_s, L_e_s, L_c_s, Y_c_s,          # controls 12–16
        C_s, λ_s, Y_s,                               # controls 17–19
        x_c_s, δ_e_sym,                              # controls 20–21 (new)
        labor_prod_s, C_R_s, Y_R_s, Y_cR_s, w_R_s, ls_s  # controls 22–27
    ]]

    SS = vcat(SS_block, SS_block)   # current and future (same at SS)
    return SS
end

# =============================================================================
# Targets (pre-estimation; estimated parameters at placeholder/prior values)
# =============================================================================
# Note: labor_share removed — it is an outcome in the new calibration.
# p_0 and dest_end_frac are Θ_e (estimated); default values used here.

targets = TARGETS   # load from steady_state.jl constant
