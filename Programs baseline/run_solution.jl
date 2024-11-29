# Based on original code by Alvaro Salazar-Perez and Hernán D. Seoane
# Modified by Mario Silva

using MKL
using DataFrames
using NaNStatistics
using Serialization
cd(@__DIR__)
#v1.7- 
#BLAS.vendor() 
#:mkl

include("solution_functions.jl")
include("steady_state.jl")
include("impulse_response_plots.jl")
include("time_series_fun.jl")

function solution_interface(model, PAR)
    eta     =   eval_ShockVAR(PAR)
    PAR_SS  =   eval_PAR_SS(PAR)
    SS      =   eval_SS(PAR_SS)
    SS_err  =   eval_SS_error(PAR_SS, SS)
    deriv   =   eval_deriv(PAR_SS, SS)
    SS_max = maximum(abs.(SS_err))
    argmax(abs.(SS_err))
    println("Residuals: $SS_max")

    ss = NamedTuple(zip(model.varnames, exp.(SS[1:model.nvar])))
    # dev = zeros(model.nvar)
    # common_keys = intersect(keys(ss), keys(steady))
    # for (i, field) in enumerate(common_keys)
    #     dev[i] = ss[field] -steady[field]
    # end
    # print(maximum(abs.(dev)))
    # @btime sol_mat = solve_model(model, deriv, eta)
    sol_mat = solve_model(model, deriv, eta)
    println("Model solved")
    out = (SS=SS, ss=ss, eta=eta, deriv=deriv, sol_mat=sol_mat)
    return out
end

## Model
# Adjustments
    flag_order      = 1
    flag_deviation  = true
    flag_SSsolver   = false

# Parameters
    @syms f_e zbar δbar sbar b ϕ ρ σ ε A η_L F κ ξ_inv ρ_z σ_z ρ_δ σ_δ ρ_s σ_s
    parameters      = [f_e; zbar; δbar; sbar; b; ϕ; ρ; σ; ε; A; η_L; F; κ; ξ_inv; ρ_z; σ_z; ρ_δ; σ_δ; ρ_s; σ_s ]
    estimate        = []
    position        = []
    priors          = (;)

    # Transformations
    β = 1/(1+ρ)
    ξ = 1/ξ_inv
    τbar = 1 - (1-δbar)*(1-sbar)
    μ = ε/(ε-1)

# Variables
@syms  z δ s u θ q K L u v v_pret e K Q  N p N_e ν_f d_f w_int w L_e L_c Y_c C λ Y labor_prod C_R Y_R Y_cR w_R ls
@syms zp  δp sp θp qp Kp Lp up vp v_pretp ep Kp Qp Np pp N_ep ν_fp d_fp w_intp wp L_ep L_cp Y_cp Cp λp Yp labor_prod_p C_Rp Y_Rp Y_cRp w_Rp lsp


x               = [u; N; v_pret; z; δ; s] # predetermined
y               = [θ; q; L; v; e; K; Q; p; N_e; ν_f; d_f; w_int; w; L_e; L_c; Y_c; C; λ; Y; labor_prod; C_R; Y_R; Y_cR; w_R; ls]
xp              = [up; Np; v_pretp; zp; δp; sp]
yp              = [θp; qp; Lp; vp; ep; Kp; Qp; pp; N_ep; ν_fp; d_fp; w_intp; wp; L_ep; L_cp; Y_cp; Cp; λp; Yp; labor_prod_p; C_Rp; Y_Rp; Y_cRp; w_Rp; lsp]
variables       = [x; y; xp; yp]
varnames = vcat(Symbol.(x), Symbol.(y))

# Shock
@syms epsilon
ex               = [epsilon]
eta = Array([0.0; 0.0; 0.0; -σ_z; σ_δ; σ_s]) # size of state space

    
nx = length(x) 
ny = length(y)
nvar = nx + ny 
ne = length(ex) 

# Array of symbolics storing model equtions 
f = fill(Sym("x"), nvar)
# Equilibrium conditions
    # Job creation condition -> θ
    f[1]  =  κ + K/q - β*λp/λ*(1-δbar*δ)*(w_intp-wp-Kp+(1-sbar*sp)*(κ+Kp/qp))
    # Marginal revenue product (welfare-based labor prod. measure) -> w_int
    f[2] = w_int - p*z*zbar/μ
    # Wage equation -> w
    f[3] = w - (ϕ*(w_int-K+θ/(1-δbar*δ)*(K+q*κ)) +(1-ϕ)*b)
    # Value of a vacancy -> Q
    f[4] = Q - (e/F)^(ξ_inv)
    # Expected discounted difference in vacancy value -> K
    f[5] = K - (Q-β*λp/λ*(1-δbar*δ)*Qp)
    # Market tightness -> v
    f[6] = θ - v/u 
    # Vacancy filling probability -> q
    f[7] = q - A*θ^(-η_L) 
    # Aggregate labor -> L
    f[8] = L - (1-u)
    # Composition of labor -> L_c
    f[9] = L - (L_c+L_e) 
    # Relative price -> p
    f[10] = p - N^(1/(ε-1))
    # Lagrangian multiplier -> λ
    f[11] = λ - C^(-σ)
    # Firm value -> ν_f
    f[12] = ν_f - β*(1-δbar*δ)*λp/λ*(ν_fp+d_fp)
    # Retail output: resources -> Y_c
    f[13] = Y_c - p*z*zbar*L_c
    # Retail output: expenditure -> C
    f[14] = Y_c - (C+F/(1+ξ_inv)*(e/F)^(1+ξ_inv)+κ*v*q)
    # Business entrants -> N_e
    f[15] = N_e - L_e*z*zbar/f_e 
    # Firm value relative to price 
    f[16] = ν_f - p*f_e/μ 
    # Output = expenditure
    f[17] = Y - (Y_c+ν_f*N_e)
    # Output = income 
    f[18] = Y - (w_int*L+N*d_f)
    # LOM of vacancies 
    f[19] = v - (v_pret + e)
    # Predetermined vacancies 
    f[20] = v_pretp - (1-δbar*δ)*((1-q)*v+sbar*s*(1-u))
    # LOM of unemployment
    f[21] = up - ((1-(1-δbar*δ)*(θ*q))*u + (1-(1-δ*δbar)*(1-sbar*s))*(1-u))
    # LOM of firms 
    f[22] = Np - (1-δbar*δ)*(N+N_e)
    # Profits 
    #f[23] = d_f - Y_c/(N*ε)
    # Labor share of income
    f[23] = ls - w*L/Y

    # Data-consistent variables (_R)
    f[24] = labor_prod - Y/(p*L) # labor productivity
    f[25] = C_R - C/p 
    f[26] = Y_R - Y/p 
    f[27] = Y_cR - Y_c/p
    f[28] = w_R - w/p

    # Exogenous processes
    f[29]  =   log(zp) - ρ_z * log(z)
    f[30] =    log(δp) -  ρ_δ * log(δ)
    f[31] = log(sp) - ρ_s*log(s)

# Steady state     
# Values 

# Initial parameters: targets and normalizations/ leave parameters as symbolic to be populated with calibration
N_s = 1.0
w_s = 1.0
z_s = 1.0
δ_s = 1.0
s_s = 1.0

fbar = 0.41
qbar = 0.8
x_v = 0.1
labor_share = 0.66

ls_s = labor_share
p_s = N_s^(1/(ε-1))
N_es = δbar/(1-δbar)*N_s
q_s = qbar/(1-δbar)
θ_s = fbar/qbar
u_s = τbar/(τbar+(1-δbar)*(θ_s*q_s))
L_s = 1 - u_s 
v_s = θ_s*u_s
e_s = δbar*(v_s+1-u_s)
v_prets = v_s - e_s 
recruiter_share= (δbar+(ρ+δbar)*(ε-1))/(δbar+(ρ+δbar)*ε)
w_wint = labor_share/recruiter_share
w_ints = w_s/(w_wint)

# Rescale zbar to be consistent with wage=1
zbar = (μ/p_s)*w_ints
surplus_ratio = (ρ+τbar)/(1-δbar)*(1/(q_s*x_v))
K_s = (w_ints-w_s)/(1+surplus_ratio) 
κ   = (1-x_v)/x_v*K_s/q_s
f_e = (μ-1)*zbar*(L_s/N_s)*(1-δbar)/(δbar*μ+ρ)
ν_fs =p_s*f_e/μ
d_fs = (ρ+δbar)/(1-δbar)*ν_fs
L_es = N_es*f_e/zbar
L_cs = L_s - L_es 
Y_cs = p_s*zbar*L_cs
Q_s = K_s*(1+ρ)/(ρ+δbar)
C_s = Y_cs -  F/(1+ξ_inv)*(e_s/F)^(1+ξ_inv) -  κ*v_s*q_s 
λ_s = C_s^(-σ)
Y_s = Y_cs + ν_fs*N_es

# data consistent
labor_prod_s = Y_s/(p_s*L_s)
C_Rs = C_s/p_s
Y_Rs = Y_s/p_s
Y_cRs = Y_cs/p_s
w_Rs = w_s/p_s
#x               = [u; N; v_pret; z] # predetermined
#y               = [θ; q; L; v; e; K; Q; p; N_e; ν_f; d_f; w_int; w; L_e; L_c; Y_c; C; λ; Y; labor_prod]
# Vector
SS_block  = [log(x) for x in [u_s, N_s, v_prets, z_s, δ_s, s_s, θ_s, q_s, L_s, v_s, e_s, K_s, Q_s, p_s, N_es, ν_fs, d_fs, w_ints, w_s, L_es, L_cs, Y_cs, C_s, λ_s, Y_s, 
        labor_prod_s, C_Rs, Y_Rs, Y_cRs, w_Rs, ls_s]]
# vertical concatenate: represent both current and future variables
SS = vcat(SS_block, SS_block)

#PAR_SS = [ALPHA; BETA; DELTA; RHO; SIGMA; MUU; AA]
PAR_SS = parameters[:]


# Procesing the model (no adjustment needed)                
model = (parameters = parameters, estimate = estimate, estimation = position,
        npar = length(parameters), ns = length(estimate), 
        priors = priors,
        x = x, y = y, xp = xp, yp = yp, variables = variables,
        varnames=varnames, #store symbols of variable names
        nx = nx, ny = ny, nvar = nvar,
        e = ex, eta = eta,
        ne = ne,
        f = f,
        nf = nvar,
        SS = SS, PAR_SS = PAR_SS,
        flag_order = flag_order, flag_deviation = flag_deviation, flag_SSsolver = flag_SSsolver)
process_model(model)


## Solution
# Parametrization: need to handle dependent parameters 
#parameters      = [f_e; δ; s; zbar; b; ϕ; ρ; σ; ε; A; η_L; F; κ; ξ_inv; ρ_z; μ_z]
targets = (labor_share=labor_share, 
           dest_ann=0.1, 
           r_ann=0.04, 
           f =fbar, 
           η_L=0.6, 
           q=qbar, 
           sep=0.031, 
           b_ratio=0.71, 
           x_v=x_v, 
           ξ_inv=1/0.265, 
           ε=4.3, σ=1.0, 
           N=N_s, w=w_s)
cal = calibrate_labor_share(targets)


# Shock values (from Coles and Kelishomi), monthly frequency
ρ_z = 0.965
σ_z = 0.007
ρ_δ = 0.875
σ_δ = 0.042
ρ_s = 0.875
σ_s = 0.0042

@unpack  f_e, δ, s, z, b, ϕ, ρ, σ, ε, A, η_L, F, κ, ξ_inv = cal

zbar = z 
δbar = δ
sbar = s
PAR     =   [f_e; zbar; δbar; sbar; b; ϕ; ρ; σ; ε; A; η_L; F; κ; ξ_inv; ρ_z; σ_z; ρ_δ; σ_δ; ρ_s; σ_s ]

sol = solution_interface(model, PAR)
@unpack ss, SS, sol_mat, eta = sol
# Export: model, targets, PAR, sol
model_output = (model, targets, PAR, sol)
serialize("model_output.jls", model_output)


## Impulse responses ##
eta_z = zero(eta) # Tech shock
eta_δ = zero(eta) # Product destruction shock
eta_s = zero(eta) # Idiosyncratic separation shock
eta_τ = zero(eta) # Common separation shock

eta_z[4] = eta[4]
eta_δ[5] = eta[5]
eta_s[6] = eta[6]

eta_τ[5] = eta[5]
eta_τ[6] = eta[6]
flag_IR = true
flag_logdev = true
T_IR = 120 # 10 years


# Technology shock
irf_z= simulate_model(model, sol_mat, T_IR, eta_z, SS, flag_IR, flag_logdev)
irf_z = 100 .*DataFrame(irf_z, varnames)
gen_irf(irf_z)
Plots.savefig("z_shock.pdf")
savefig("z_shock.png")

# Destruction rate shock: consistent with Beveridge curve
irf_δ= simulate_model(model, sol_mat, T_IR, eta_δ, SS, flag_IR, flag_logdev) 
irf_δ = 100 .*DataFrame(irf_δ, varnames)
gen_irf(irf_δ)
Plots.savefig("dest_shock.pdf")

# Idiosyncratic job separation shock
irf_s= simulate_model(model, sol_mat, T_IR, eta_s, SS, flag_IR, flag_logdev) 
irf_s = 100 .*DataFrame(irf_s, varnames)
gen_irf(irf_s)
Plots.savefig("s_shock.pdf")

#Commmon separation shock
irf_τ= simulate_model(model, sol_mat, T_IR, eta_τ, SS, flag_IR, flag_logdev) 
irf_τ = 100 .*DataFrame(irf_τ, varnames)
gen_irf(irf_τ)
Plots.savefig("common_shock.pdf")

### Impulse response comparison ###
# ε
targets2 = (targets..., ε=100.0 )
cal2 = calibrate_labor_share(targets2)
@unpack  f_e, δ, s, z, b, ϕ, ρ, σ, ε, A, η_L, F, κ, ξ_inv = cal2
PAR2     =   [f_e; z; δ; s; b; ϕ; ρ; σ; ε; A; η_L; F; κ; ξ_inv; ρ_z; σ_z; ρ_δ;
 σ_δ; ρ_s; σ_s ]
sol2 = solution_interface(model, PAR2)
sol_mat2 = sol2.sol_mat
SS2 = sol2.SS

sim_IR2 = simulate_model(model, sol_mat2, T_IR, eta_z, SS2, flag_IR, flag_logdev)
irf_z2 = 100 .*DataFrame(sim_IR2, varnames)
gen_irf_comp(irf_z, irf_z2, ["Baseline", " ε=100"])
Plots.savefig("irf_comp_epsi.pdf")

# High δ calibration 
targets3 = (targets..., dest_ann=0.2)
cal3 = calibrate_labor_share(targets3)
@unpack  f_e, δ, s, z, b, ϕ, ρ, σ, ε, A, η_L, F, κ, ξ_inv = cal3
PAR3     =   [f_e; s; zbar; δbar; sbar; b; ϕ; ρ; σ; ε; A; η_L; F; κ; ξ_inv; ρ_z; σ_z; ρ_δ; σ_δ; ρ_s; σ_s ]
sol3 = solution_interface(model, PAR3)
sol_mat3 = sol3.sol_mat
SS3 = sol3.SS

sim_IR3 = simulate_model(model, sol_mat3, T_IR, eta_z, SS2, flag_IR, flag_logdev)
irf_z3 = 100 .*DataFrame(sim_IR3, varnames)
gen_irf_comp(irf_z, irf_z3, ["Baseline", " 20% annual destruction rate"])
Plots.savefig("irf_comp_delta.pdf")

# High sunk vacancy posting costs
targets4 = (targets..., x_v=0.4)
cal4 = calibrate_labor_share(targets4)
@unpack  f_e, δ, s, z, b, ϕ, ρ, σ, ε, A, η_L, F, κ, ξ_inv = cal4
PAR4     =   [f_e; s; zbar; δbar; sbar; b; ϕ; ρ; σ; ε; A; η_L; F; κ; ξ_inv; ρ_z; σ_z; ρ_δ; σ_δ; ρ_s; σ_s ]
sol4 = solution_interface(model, PAR4)
sol_mat4 = sol4.sol_mat
SS4 = sol4.SS

sim_IR4 = simulate_model(model, sol_mat4, T_IR, eta_z, SS2, flag_IR, flag_logdev)
irf_z4 = 100 .*DataFrame(sim_IR4, varnames)
gen_irf_comp(irf_z, irf_z3, ["Baseline", " 20% annual destruction rate"])
Plots.savefig("irf_comp_delta.pdf")