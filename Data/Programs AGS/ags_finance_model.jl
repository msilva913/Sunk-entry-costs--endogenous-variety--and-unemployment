# ags_finance_model.jl
# =============================================================================
# AGS with the financial friction of Notes/ags_equilibrium.md §10.3 (branch d15_ags_finance).
# Standalone: the toolkit's derivative functions are global, so this model cannot share a Julia
# process with ags_model.jl. Monthly; timing, labor block and household as in ags_model.jl.
#
# The friction (§10.3). Lines are owned by the household. Within month t a line
#   1. holds internal funds n_t (aggregate W_t, a state);
#   2. buys new positions at P_t = x_m e_t^ξ_inv (a market price) and borrows the rest within
#      the month: B_t = P_t e_t − W_t, leverage b_t = B_t/W_t;
#   3. produces; its cash flow Π_t = (z̄z_t − w_t)(1 − u_t) − X^r_t − κ q_t v_t pays hiring and
#      reactivation (only creation is financed);
#   4. is destroyed with probability δ_t δ̄: it defaults and its cash is lost (a deadweight
#      liquidation loss, §10.3); lenders price this, R_t = Ψ(b_t)/(1 − δ_t δ̄);
#   5. a survivor repays and pays out ϑ_t W_t at cost φ(ϑ_t)W_t, φ(ϑ) = ϑ + (κ_d/2)(ϑ − ϑ̄)².
# Premium Ψ(b) = 1 + prem, prem = ψ̄(b/b̄)^χ (deadweight agency cost). Marginal premia:
#   borrowing at given funds  Ψ^m = 1 + prem(1 + χ),
#   funds at given outlays    Ψ^n = 1 + prem(1 + χ(1 + b)).
# Optimality: payout λ_t φ'(ϑ_t) = 1; retention λ_t = E_t[m_{t+1} ω_{t+1}]; value of funds
#   ω_t = λ_t[Ψ^n(b_t) + (1 − δ_tδ̄)(κ_d/2)(ϑ_t² − ϑ̄²)]; creation Q_t = Ψ^m(b_t) x_m e_t^ξ_inv;
#   firms discount with m̃_{t+1} = m_{t+1} λ_{t+1}/λ_t (f[7], f[8], f[11]). The wage equation is
#   unchanged (bargaining in the line's cash units; an open choice, §10.3).
# Funds: W_{t+1} = (1 − δ_tδ̄)[Π_t − φ(ϑ_t)W_t] − Ψ(b_t)B_t (replacement lines start with no funds).
#
# Calibration: AGS's (calibrate_ags), with x_m divided by Ψ̄^m so that the labor-market steady
# state (θ, u, v, Q, K, w, ...) is identical to AGS. Steady state: λ̄ = 1, Ψ^n(b̄) = 1/β, so
# ψ̄ = r/(1 + χ(1 + b̄)); ϑ̄ from the stationarity of W.
# =============================================================================

cd(@__DIR__)
const BASELINE_DIR = joinpath(@__DIR__, "..", "Programs baseline")
include(joinpath(BASELINE_DIR, "solution_functions.jl"))
include(joinpath(BASELINE_DIR, "time_series_fun.jl"))
include("ags_steady_state.jl")
using DataFrames, Statistics, Printf

const SHOCK_CSV = joinpath(@__DIR__, "..", "Bartek analysis", "data", "results", "shock_calibration.csv")
const SHK = shock_processes(SHOCK_CSV; prefix = "bprime")
const LAM = 1_600.0

"""Friction targets (§10.3): leverage b̄ (one third of creation outlays borrowed, as in
Petrosky-Nadeau 2014), slope χ of the premium in leverage, payout cost κ_d (per unit of the
monthly payout ratio; 0 = funds move freely)."""
const FIN_TARGETS = (b_bar = 0.5, χ = 1.0, κ_d = 1.0)

flag_order, flag_deviation, flag_SSsolver = 1, true, false

@syms δbar sbar r A η_L b ϕ κ x_m ξ_inv zbar a p_r ψ_r σ ψbar χ bbar κ_d ϑbar ρ_z σ_z ρ_δ σ_δ ρ_s σ_s
parameters = [δbar; sbar; r; A; η_L; b; ϕ; κ; x_m; ξ_inv; zbar; a; p_r; ψ_r; σ;
              ψbar; χ; bbar; κ_d; ϑbar; ρ_z; σ_z; ρ_δ; σ_δ; ρ_s; σ_s]
estimate, position, priors = [], [], (;)

build_PAR(cal, shk) = [cal.δ; cal.s; cal.r; cal.A; cal.η_L; cal.b; cal.ϕ; cal.κ; cal.x_m;
                       cal.ξ_inv; cal.zbar; cal.a; cal.p_r; cal.ψ_r; cal.σ;
                       cal.ψbar; cal.χ; cal.b_bar; cal.κ_d; cal.ϑbar;
                       shk.ρz; shk.σz; shk.ρδ; shk.σδ; shk.ρs; shk.σs]

@syms Et Vu Vs W z δ s Etp Vup Vsp Wp zp δp sp
@syms u θ q v v_pre e Q K Λ_r Q_rep H w C V Rf B ϑ λ ω prem
@syms up θp qp vp v_prep ep Qp Kp Λ_rp Q_repp Hp wp Cp Vp Rfp Bp ϑp λp ωp premp
x  = [Et; Vu; Vs; W; z; δ; s]
y  = [u; θ; q; v; v_pre; e; Q; K; Λ_r; Q_rep; H; w; C; V; Rf; B; ϑ; λ; ω; prem]
xp = [Etp; Vup; Vsp; Wp; zp; δp; sp]
yp = [up; θp; qp; vp; v_prep; ep; Qp; Kp; Λ_rp; Q_repp; Hp; wp; Cp; Vp; Rfp; Bp; ϑp; λp; ωp; premp]
variables = [x; y; xp; yp]
varnames  = vcat(Symbol.(x), Symbol.(y))
nx, ny = length(x), length(y)
nvar, ne = nx + ny, 3
eta = [0 0 0; 0 0 0; 0 0 0; 0 0 0; σ_z 0 0; 0 σ_δ 0; 0 0 σ_s]
@syms epsilon
ex = [epsilon]

function fin_equations()
    f    = fill(Sym("x"), nvar)
    β    = 1 / (1 + r)
    m    = β * (Cp / C)^(-σ)                # household discount factor
    mt   = m * λp / λ                       # the lines' discount factor m̃
    S    = 1 - δ * δbar
    R_shp = Qp - Kp - Hp
    ψ_cr = ψ_r / (ψ_r + 1)
    X_r  = Vs * (Q * Λ_r - Q_rep)
    Xc   = x_m * e^(1 + ξ_inv) / (1 + ξ_inv) + κ * q * v     # resource cost of creation and hiring
    Π    = (zbar * z - w) * (1 - u) - X_r - κ * q * v        # lines' operating cash flow
    φ    = ϑ + κ_d / 2 * (ϑ - ϑbar)^2                        # payout cost per unit of funds
    lev  = B / W
    f[1]  = u - (1 - Et)
    f[2]  = θ - v / u
    f[3]  = q - A * θ^(-η_L)
    f[4]  = v - (v_pre + e)
    f[5]  = v_pre - (Vu + Λ_r * Vs)
    f[6]  = Q - (1 + prem * (1 + χ)) * x_m * e^ξ_inv           # creation: Q = Ψ^m P
    f[7]  = K - (Q - mt * S * Qp)
    f[8]  = κ + K / q - mt * S * (
                (1 - ϕ) * (zbar * zp - Kp - b) - ϕ * θp * (Kp + qp * κ)
                + (1 - sp * sbar) * (κ + Kp / qp) - (1 - ϕ) * sp * sbar * R_shp)
    f[9]  = Λ_r - ((1 - p_r) + p_r * (Q / a)^ψ_r)
    f[10] = Q_rep - Q * (Λ_r - ψ_cr * (Λ_r - 1 + p_r))
    f[11] = H - mt * S * Q_repp
    f[12] = w - (ϕ * (zbar * z - K + θ * (K + q * κ) - s * sbar * (Q - K - H)) + (1 - ϕ) * b)
    f[13] = Etp - S * ((1 - s * sbar) * (1 - u) + A * θ^(1 - η_L) * u)
    f[14] = Vup - S * (1 - q) * v
    f[15] = Vsp - S * s * sbar * (1 - u)
    f[16] = Wp - (S * (Π - φ * W) - (1 + prem) * B)             # internal funds
    f[17] = log(zp) - ρ_z * log(z)
    f[18] = log(δp) - ρ_δ * log(δ)
    f[19] = log(sp) - ρ_s * log(s)
    # Resources: output less creation, hiring and reactivation, the agency cost, the payout
    # cost, and destroyed lines' cash (liquidation loss)
    f[20] = C - (zbar * z * (1 - u) - Xc - X_r - prem * B
                 - S * κ_d / 2 * (ϑ - ϑbar)^2 * W - (1 - S) * Π)
    f[21] = V - m * (Vp + Cp - wp * (1 - up))                  # claim to all non-labor income
    f[22] = 1 - Rf * m
    f[23] = B - (x_m * e^(1 + ξ_inv) - W)                      # borrowing: outlays P e less funds
    f[24] = prem - ψbar * (lev / bbar)^χ                        # premium Ψ − 1
    f[25] = λ * (1 + κ_d * (ϑ - ϑbar)) - 1                      # payout: λ φ'(ϑ) = 1
    f[26] = λ - m * ωp                                          # retention
    f[27] = ω - λ * (1 + prem * (1 + χ * (1 + lev)) + S * κ_d / 2 * (ϑ^2 - ϑbar^2))
    return f
end

model = (parameters = parameters, estimate = estimate, estimation = position,
         npar = length(parameters), ns = length(estimate), priors = priors,
         x = x, y = y, xp = xp, yp = yp, variables = variables, varnames = varnames,
         nx = nx, ny = ny, nvar = nvar, e = ex, eta = eta, ne = ne,
         f = fin_equations(), nf = nvar, SS = Sym.(zeros(2nvar)), PAR_SS = parameters[:],
         flag_order = flag_order, flag_deviation = flag_deviation, flag_SSsolver = flag_SSsolver)
process_model(model)

"""AGS calibration plus the friction: ψ̄ from Ψ^n(b̄) = 1/β, x_m rescaled by 1/Ψ̄^m, and the
steady state (labor block from steady_state_ags at AGS's x_m, then the financial block)."""
function calibrate_fin(tg = AGS_TARGETS, ft = FIN_TARGETS)
    cal0 = calibrate_ags(tg)
    st0  = steady_state_ags(cal0)
    @unpack b_bar, χ, κ_d = ft
    @unpack δ, r, ξ_inv, κ, zbar = cal0
    ψbar = r / (1 + χ * (1 + b_bar))
    Ψm   = 1 + ψbar * (1 + χ)
    x_m  = cal0.x_m / Ψm
    S    = 1 - δ
    O    = x_m * st0.e^(1 + ξ_inv)                     # creation outlays P e
    W    = O / (1 + b_bar)
    B    = O - W
    Π    = (zbar - st0.w) * (1 - st0.u) - st0.X_r - κ * st0.q * st0.v
    ϑbar = (S * Π - (1 + ψbar) * B - W) / (S * W)      # stationarity of W, φ(ϑ̄) = ϑ̄
    if ϑbar <= 0
        error("calibrate_fin: steady-state payout ratio ϑ̄ = $ϑbar ≤ 0.")
    end
    Xc   = x_m * st0.e^(1 + ξ_inv) / (1 + ξ_inv) + κ * st0.q * st0.v
    C    = zbar * (1 - st0.u) - Xc - st0.X_r - ψbar * B - δ * Π
    V    = (C - st0.w * (1 - st0.u)) / r
    cal  = merge(cal0, (; x_m, ψbar, χ, b_bar, κ_d, ϑbar))
    fin  = (; W, B, Π, ϑ = ϑbar, λ = 1.0, ω = 1 + r, prem = ψbar, C, V, Xc, O)
    return cal, st0, fin
end

"""Log steady state in the order of `varnames` (current and next period)."""
function SS_fin(cal, st, fin)
    Vs = (1 - cal.δ) * cal.s * (1 - st.u)
    Vu = (1 - cal.δ) * (1 - st.q) * st.v
    blk = log.([1 - st.u, Vu, Vs, fin.W, 1.0, 1.0, 1.0,
                st.u, st.θ, st.q, st.v, st.v_pre, st.e, st.Q, st.K, st.Λ_r, st.Qrep, st.H, st.w,
                fin.C, fin.V, 1 + cal.r, fin.B, fin.ϑ, fin.λ, fin.ω, fin.prem])
    return vcat(blk, blk)
end

"""Solve at targets; errors on a steady-state residual or a failed Blanchard–Kahn count."""
function solve_fin(tg = AGS_TARGETS, ft = FIN_TARGETS; shk = SHK)
    cal, st, fin = calibrate_fin(tg, ft)
    PAR = build_PAR(cal, shk)
    SS  = SS_fin(cal, st, fin)
    PAR_SS = Base.invokelatest(eval_PAR_SS, PAR)
    err = maximum(abs, Base.invokelatest(eval_SS_error, PAR_SS, SS))
    if !(err < 1e-10)
        error("Finance steady-state residual $err")
    end
    deriv = Base.invokelatest(eval_deriv, PAR_SS, SS)
    sol   = solve_model(model, deriv, Base.invokelatest(eval_ShockVAR, PAR))
    if !(sol.qzflag == 1)
        error("Blanchard–Kahn failed (qzflag = $(sol.qzflag))")
    end
    return (cal = cal, st = st, fin = fin, SS = SS, sol_mat = sol)
end

shock_col(i, sd) = (c = zeros(nx, 1); c[4 + i, 1] = sd; c)     # i = 1 z, 2 δ, 3 s

"""Monthly impulse response (log deviations) to a one-SD shock i, T months after impact."""
irf_fin(sol, i, T; shk = SHK) = DataFrame(simulate_model((; model..., ne = 1), sol.sol_mat, T + 1,
        shock_col(i, (shk.σz, shk.σδ, shk.σs)[i]), sol.SS, true, true), varnames)[2:end, :]

"""Quarterly averages of a monthly vector (first nq quarters)."""
quarterly(m, nq) = [mean(m[3k-2:3k]) for k in 1:nq]

"""Model counterparts of the data series, quarterly averages (as model_quarterly in
ags_model.jl), plus the spread: prem in basis points, annualized (12 × monthly)."""
function model_quarterly(sim, cal)
    m = DataFrame()
    m.u = sim.u
    m.v = sim.v
    m.θ = sim.θ
    m.z = sim.z
    m.s = sim.s
    m.f = log.(1 .- cal.δ .* exp.(sim.δ)) .- log(1 - cal.δ) .+ (1 - cal.η_L) .* sim.θ
    m.δ = vcat(0.0, sim.δ[1:end-1])
    m.ebp = 12e4 .* cal.ψbar .* (exp.(sim.prem) .- 1)          # deviation in bp, annualized
    return monthly_to_quarterly(m)
end
hp_frame(q) = DataFrame([c => hp_filter(Float64.(q[!, c]), LAM) for c in names(q)])
model_cycles(sim, cal) = hp_frame(model_quarterly(sim, cal))
