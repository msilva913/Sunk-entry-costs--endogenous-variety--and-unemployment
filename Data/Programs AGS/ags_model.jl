# ags_model.jl
# =============================================================================
# The AGS model with costly reactivation (Part 1 of D15; equations in Notes/ags_equilibrium.md;
# calibration in ags_steady_state.jl): symbolic model, first-order solution, and simulation
# helpers. Toolkit by Salazar-Perez and Seoane (Programs baseline/solution_functions.jl).
# Monthly. Included by ags_run_solution.jl, the sweeps and ags_svar_sim.jl, as the baseline's
# run_solution_core.jl is included by its runners.
#
# Timing within t (the baseline's, run_solution_core.jl):
#   1. z_t, δ_t, s_t are realized.
#   2. Positions vacated at the end of t−1 are reactivated (rate Λ_r,t = F_r(Q_t)) or retired;
#      new positions e_t are created.
#   3. Production; u_t workers and v_t = v_pre,t + e_t vacancies match.
#   4. Separation s_t and destruction δ_t take effect.
# States: Ẽ_t = employment entering t; Ṽ^u_t = unmatched vacancies carried over;
#   Ṽ^s_t = positions vacated by separation at the end of t−1 (reactivated at t or retired);
#   the shocks z, δ, s (multiplicative deviations around 1; δ·δ̄, s·s̄ are the rates).
# v_pre,t = Ṽ^u_t + Λ_r,t Ṽ^s_t depends on Q_t, so v_pre is a control, not a state.
# =============================================================================

cd(@__DIR__)
const BASELINE_DIR = joinpath(@__DIR__, "..", "Programs baseline")
include(joinpath(BASELINE_DIR, "solution_functions.jl"))
include(joinpath(BASELINE_DIR, "time_series_fun.jl"))
include("ags_steady_state.jl")
using DataFrames, Statistics, Printf, Plots

# Same values as SHOCKS_BASE in run_solution_core.jl (S14, filter-consistent; part6b).
const SHK = (ρz = 0.9628, σz = 0.0062, ρδ = 0.5319, σδ = 0.0641, ρs = 0.8048, σs = 0.0437)
const T_IR, T_SM, LAM = 60, 60_000, 1_600.0

# =============================================================================
# Model
# =============================================================================
flag_order, flag_deviation, flag_SSsolver = 1, true, false

@syms δbar sbar r A η_L b ϕ κ x_m ξ_inv zbar a p_r ψ_r σ ρ_z σ_z ρ_δ σ_δ ρ_s σ_s
parameters = [δbar; sbar; r; A; η_L; b; ϕ; κ; x_m; ξ_inv; zbar; a; p_r; ψ_r; σ;
              ρ_z; σ_z; ρ_δ; σ_δ; ρ_s; σ_s]
estimate, position, priors = [], [], (;)

"""PAR in the order of `parameters`. p_r = 0 is costless reactivation (Λ_r ≡ 1, Q^rep = Q)."""
build_PAR(cal, shk) = [cal.δ; cal.s; cal.r; cal.A; cal.η_L; cal.b; cal.ϕ; cal.κ; cal.x_m;
                       cal.ξ_inv; cal.zbar; cal.a; cal.p_r; cal.ψ_r; cal.σ;
                       shk.ρz; shk.σz; shk.ρδ; shk.σδ; shk.ρs; shk.σs]

@syms Et Vu Vs z δ s Etp Vup Vsp zp δp sp
@syms u θ q v v_pre e Q K Λ_r Q_rep H w C V Rf
@syms up θp qp vp v_prep ep Qp Kp Λ_rp Q_repp Hp wp Cp Vp Rfp
x  = [Et; Vu; Vs; z; δ; s]
y  = [u; θ; q; v; v_pre; e; Q; K; Λ_r; Q_rep; H; w; C; V; Rf]
xp = [Etp; Vup; Vsp; zp; δp; sp]
yp = [up; θp; qp; vp; v_prep; ep; Qp; Kp; Λ_rp; Q_repp; Hp; wp; Cp; Vp; Rfp]
variables = [x; y; xp; yp]
varnames  = vcat(Symbol.(x), Symbol.(y))
nx, ny = length(x), length(y)
nvar, ne = nx + ny, 3
eta = [0 0 0; 0 0 0; 0 0 0; σ_z 0 0; 0 σ_δ 0; 0 0 σ_s]
@syms epsilon
ex = [epsilon]

function ags_equations()
    f    = fill(Sym("x"), nvar)
    β    = 1 / (1 + r)
    m    = β * (Cp / C)^(-σ)                # household discount factor m_{t+1}; β when σ = 0
    S    = 1 - δ * δbar                     # survival of destruction from t to t+1
    R_sh, R_shp = Q - K - H, Qp - Kp - Hp   # 𝓡 = E[m(1−δ)(Q' − Q^rep')], by f[7] and f[11]
    ψ_cr = ψ_r / (ψ_r + 1)
    # Investment: creating positions and hiring (X), reactivating vacated positions (X^r)
    X    = e * Q / (1 + ξ_inv) + κ * q * v
    X_r  = Vs * (Q * Λ_r - Q_rep)
    Xp   = ep * Qp / (1 + ξ_inv) + κ * qp * vp
    X_rp = Vsp * (Qp * Λ_rp - Q_repp)
    Dp   = (zbar * zp - wp) * (1 - up) - Xp - X_rp     # dividends at t+1
    # Post-reactivation stocks and matching
    f[1]  = u - (1 - Et)                                        # unemployment entering matching
    f[2]  = θ - v / u                                           # tightness
    f[3]  = q - A * θ^(-η_L)                                    # vacancy-filling rate
    f[4]  = v - (v_pre + e)                                     # E1: vacancy pool
    f[5]  = v_pre - (Vu + Λ_r * Vs)                             # inherited vacancies
    # Positions: creation, user cost, reactivation
    f[6]  = Q - x_m * e^ξ_inv                                   # creation cost Q(e)
    f[7]  = K - (Q - m * S * Qp)                                # user cost
    f[8]  = κ + K / q - m * S * (                               # E2: job creation
                (1 - ϕ) * (zbar * zp - Kp - b) - ϕ * θp * (Kp + qp * κ)
                + (1 - sp * sbar) * (κ + Kp / qp) - (1 - ϕ) * sp * sbar * R_shp)
    f[9]  = Λ_r - ((1 - p_r) + p_r * (Q / a)^ψ_r)               # Λ_r = F_r(Q); 1 if p_r = 0
    f[10] = Q_rep - Q * (Λ_r - ψ_cr * (Λ_r - 1 + p_r))          # E max(Q − c, 0); Q if p_r = 0
    f[11] = H - m * S * Q_repp                                  # H = E[m(1−δ) Q^rep']
    f[12] = w - (ϕ * (zbar * z - K + θ * (K + q * κ) - s * sbar * R_sh) + (1 - ϕ) * b)  # wage
    # Laws of motion (E3, E4)
    f[13] = Etp - S * ((1 - s * sbar) * (1 - u) + A * θ^(1 - η_L) * u)
    f[14] = Vup - S * (1 - q) * v
    f[15] = Vsp - S * s * sbar * (1 - u)
    # Shocks
    f[16] = log(zp) - ρ_z * log(z)
    f[17] = log(δp) - ρ_δ * log(δ)
    f[18] = log(sp) - ρ_s * log(s)
    # Household (Notes/ags_equilibrium.md §2): resources, the share Euler, the bond Euler
    f[19] = C - (zbar * z * (1 - u) - X - X_r)                  # resource constraint
    f[20] = V - m * (Vp + Dp)                                   # shares: V = E[m(V' + D')]
    f[21] = 1 - Rf * m                                          # riskless bond: 1 = R^f E[m]
    return f
end

model = (parameters = parameters, estimate = estimate, estimation = position,
         npar = length(parameters), ns = length(estimate), priors = priors,
         x = x, y = y, xp = xp, yp = yp, variables = variables, varnames = varnames,
         nx = nx, ny = ny, nvar = nvar, e = ex, eta = eta, ne = ne,
         f = ags_equations(), nf = nvar, SS = Sym.(zeros(2nvar)), PAR_SS = parameters[:],
         flag_order = flag_order, flag_deviation = flag_deviation, flag_SSsolver = flag_SSsolver)
process_model(model)

"""Log steady state in the order of `varnames` (current and next period)."""
function SS_ags(cal)
    st = steady_state_ags(cal)
    Vs = (1 - cal.δ) * cal.s * (1 - st.u)
    Vu = (1 - cal.δ) * (1 - st.q) * st.v
    blk = log.([1 - st.u, Vu, Vs, 1.0, 1.0, 1.0,
                st.u, st.θ, st.q, st.v, st.v_pre, st.e, st.Q, st.K, st.Λ_r, st.Qrep, st.H, st.w,
                st.C, st.V, st.Rf])
    return vcat(blk, blk), st
end

"""Solve at targets `tg`; errors on a steady-state residual or a failed Blanchard–Kahn count."""
function solve_ags(tg; shk = SHK)
    cal = calibrate_ags(tg)
    PAR = build_PAR(cal, shk)
    SS, st = SS_ags(cal)
    PAR_SS = Base.invokelatest(eval_PAR_SS, PAR)
    err = maximum(abs, Base.invokelatest(eval_SS_error, PAR_SS, SS))
    if !(err < 1e-10)
        error("AGS steady-state residual $err")
    end
    deriv = Base.invokelatest(eval_deriv, PAR_SS, SS)
    sol   = solve_model(model, deriv, Base.invokelatest(eval_ShockVAR, PAR))
    if !(sol.qzflag == 1)
        error("Blanchard–Kahn failed (qzflag = $(sol.qzflag))")
    end
    return (cal = cal, st = st, SS = SS, sol_mat = sol)
end

shock_col(i, sd) = (c = zeros(nx, 1); c[3 + i, 1] = sd; c)     # i = 1 z, 2 δ, 3 s
const SDS = (SHK.σz, SHK.σδ, SHK.σs)

irf(sol, i) = 100 .* DataFrame(simulate_model((; model..., ne = 1), sol.sol_mat, T_IR,
                                              shock_col(i, SDS[i]), sol.SS, true, true), varnames)

hpq(sim, v) = hp_filter(Float64.(monthly_to_quarterly(sim[!, [v]])[!, v]), LAM)
function moments(sol, cols)
    eta_m = hcat([shock_col(i, SDS[i]) for i in cols]...)
    sim = DataFrame(simulate_model((; model..., ne = length(cols)), sol.sol_mat, T_SM, eta_m,
                                   sol.SS, false, true), varnames)
    hu, hv = hpq(sim, :u), hpq(sim, :v)
    return (c = cor(hu, hv), su = 100std(hu), sv = 100std(hv))
end

# =============================================================================
# Model counterparts of the data series (quarterly)
# =============================================================================

"""
    model_quarterly(sim, cal) → DataFrame

Quarterly averages of the model counterparts of the data series, in log deviations, from a
monthly simulation `sim` (log deviations, columns `varnames`). Not filtered.
- u, v, θ, z (labor productivity: output per worker is z̄·z), s (the separation rate);
- f = (1 − δ_t) A θ_t^(1−η), the gross job-finding rate (the data's jf);
- δ: measured exit δ_e,t = δ_{t−1} (BED Deaths dates destruction one month later).
"""
function model_quarterly(sim, cal)
    m = DataFrame()
    m.u = sim.u
    m.v = sim.v
    m.θ = sim.θ
    m.z = sim.z
    m.s = sim.s
    m.f = log.(1 .- cal.δ .* exp.(sim.δ)) .- log(1 - cal.δ) .+ (1 - cal.η_L) .* sim.θ
    m.δ = vcat(0.0, sim.δ[1:end-1])
    return monthly_to_quarterly(m)
end

"""HP-1600 cycles of every column of a quarterly frame."""
hp_frame(q) = DataFrame([c => hp_filter(Float64.(q[!, c]), LAM) for c in names(q)])

"""Quarterly HP cycles of the model counterparts (model_quarterly, then hp_frame)."""
model_cycles(sim, cal) = hp_frame(model_quarterly(sim, cal))
