using Plots
using Parameters, CSV, StatsBase, Statistics, Random
using NLsolve
using Roots, Optim, LeastSquaresOptim
using PrettyPrinting
using Distributions
cd(@__DIR__)
include("steady_state.jl")


targets = (labor_share=0.66, dest_ann=0.06, r_ann=0.04, f =0.41, η_L=0.6, q=0.8, sep=0.031, b_ratio=0.71, x_v=0.2, ξ_inv=1, ε=4, σ=1.5, N=1, w=1)
cal = calibrate_labor_share(targets)
steady = steady_state(cal)

T = 100
θ_grid = range(1e-10, stop=1.5, length=T)

N_jcc_grid = N_jcc.(θ_grid, Ref(cal))
N_res_grid = N_res.(θ_grid, Ref(cal))
   
# Initial plot of $\theta$ and $N$
# Initial plot
plt = plot(θ_grid, N_jcc_grid, label="N: job creation condition", linewidth=2)  # Adjust the alpha value here
plot!(θ_grid, N_res_grid, label="N: aggregate resource constraint", linewidth=2)  # Adjust the alpha value here
xlabel!(L"θ")
ylabel!(L"N")
vline!([steady.θ], line=:dash, color=:black, label="Steady "* L"θ", linewidth=0.5)  # Adjust the alpha value here
hline!([steady.N], line=:dash, color=:black, label="Steady "*"N", linewidth=0.5)  # Adjust the alpha value here
plot!(legend=true)
savefig("equilibrium_curves.pdf")
display(plt)

# Check for second Equilibrium
#steady = steady_state(cal, init=0.01)

#N_res is increasing with markup (decreasing with ε)-> incentive to place more resources in business formation
para = (cal..., ε=cal.ε*0.9)
N_res_grid2 = N_res.(θ_grid, Ref(para))

fig, ax = plt.subplots()
ax.plot(θ_grid, N_res_grid, label="N: aggregate resource constraint")
ax.plot(θ_grid, N_res_grid2, label="N: aggregate resource constraint, higher markup μ")
ax.set_xlabel("θ")
ax.set_ylabel("N")
ax.legend()
display(fig)

#N_res decreasing with δ (reduces labor resources, increases effective discounting, decreases ratio of N to N_e)

para = (cal..., δ=cal.δ*1.1)
N_res_grid2 = N_res.(θ_grid, Ref(para))

fig, ax = plt.subplots()
ax.plot(θ_grid, N_res_grid, label="N: aggregate resource constraint")
ax.plot(θ_grid, N_res_grid2, label="N: aggregate resource constraint, higher destruction rate δ")
ax.set_xlabel("θ")
ax.set_ylabel("N")
ax.legend()
display(fig)


# Baseline exercise: increase in f_e

para = (cal..., f_e=cal.f_e*1.2)
steady2 = steady_state(para)
N_jcc_grid2 = N_jcc.(θ_grid, Ref(para))
N_res_grid2 = N_res.(θ_grid, Ref(para))

# New plot
fig, ax = plt.subplots()
ax.plot(θ_grid, N_jcc_grid, label="N: job creation condition")
ax.plot(θ_grid, N_res_grid, label="N: aggregate resource constraint")
ax.plot(θ_grid, N_res_grid2, label="N: aggregate resource constraint, higher sunk entry cost f_e")
ax.axvline(steady.θ, linestyle="--", color="black", linewidth=0.5)
ax.axhline(steady.N, linestyle="--", color="black", linewidth=0.5)
ax.axvline(steady2.θ, linestyle="--", color="red", linewidth=0.5)
ax.axhline(steady2.N, linestyle="--", color="red", linewidth=0.5)
ax.legend()
display(fig)

# Increase in F, amount of firms which can create vacancies

para = (cal..., F=cal.F*1.2)
steady2 = steady_state(para)
N_jcc_grid2 = N_jcc.(θ_grid, Ref(para))
N_res_grid2 = N_res.(θ_grid, Ref(para)) # only shifts N_jcc



# New plot
fig, ax = plt.subplots()
ax.plot(θ_grid, N_jcc_grid, label="N: job creation condition")
ax.plot(θ_grid, N_jcc_grid2, label="N: job creation condition, higher F")
ax.plot(θ_grid, N_res_grid, label="N: aggregate resource constraint")
ax.axvline(steady.θ, linestyle="--", color="black", linewidth=0.5)
ax.axhline(steady.N, linestyle="--", color="black", linewidth=0.5)
ax.axvline(steady2.θ, linestyle="--", color="red", linewidth=0.5)
ax.axhline(steady2.N, linestyle="--", color="red", linewidth=0.5)
ax.legend()
display(fig)



# function comparative_statics(Para, sym::Symbol; T=10, scal=1.5)

#     # Create copy of mutable struct instance
#     para = deepcopy(Para)
#     par = getfield(para, sym)
#     par_shock = range(1.0, stop=scal, length=T)
#     par_shock = par_shock.*par

#     θ_seq = zeros(T)
#     u_seq = zeros(T)
#     v_seq = zeros(T)
#     e_seq = zeros(T)
#     N_seq = zeros(T)
#     K_seq = zeros(T)

#     for i in 1:T
#         setfield!(para, sym, par_shock[i])
#         res, sol = steady_state(para)
#         θ_seq[i] = res.θ
#         u_seq[i] = res.u
#         v_seq[i] = res.v 
#         e_seq[i] = res.e
#         N_seq[i] = res.N 
#         K_seq[i] = res.K
#     end

#     out = [100*log.(x_seq./x_seq[1]) for x_seq in [par_shock, θ_seq, u_seq, v_seq, e_seq, N_seq, K_seq]]
#     par_shock_per, θ_per, u_per, v_per, e_per, N_per, K_per = out

#     res = (par_shock_per=par_shock_per, θ_per=θ_per, u_per=u_per, v_per=v_per, e_per=e_per, N_per=N_per, K_per=K_per)
#     return res
# end


# res = comparative_statics(para, :f_e) # Increase in firm entry costs
# res = comparative_statics(para, :F, scal=0.5)
# res = comparative_statics(para, :δ, scal=2.0)
# res = comparative_statics(para, :ε, scal=2.0)
# res = comparative_statics(para, :z, scal=1.5)


# fig, ax = plt.subplots()
# ax.plot(res.par_shock_per, res.u_per, label="u")
# ax.plot(res.par_shock_per, res.θ_per, label="θ")
# ax.plot(res.par_shock_per, res.v_per, label="v")
# ax.plot(res.par_shock_per, res.N_per, label="N")
# ax.plot(res.par_shock_per, res.e_per, label="e")
# ax.plot(res.par_shock_per, res.K_per, label="K")
# ax.plot(res.par_shock_per, res.K_per, label="C")
# ax.legend()
# display(fig)

# # Decrese in number of recruiters who can post vacancies

function ϵ_e_fun(θ, para)
    @unpack τ, δ, A, η_L = para
    f = jf(θ, A, η_L)
    #ϵ_e_θ = (θ*τ + (1-δ)*(1-η_L)*f)/(θ*τ+(1-δ)*f) - ((1-δ)*(1-η_L)*f)/(τ+(1-δ)*f)
    u = 1 - L_fun(θ, para)
    ϵ_e_θ = u*(τ+(1-δ)*f)/(θ*τ+(1-δ)*f)*(θ+(1-η_L)*(1-θ)*(1-u))
    return ϵ_e_θ
end

# function ϵ_e_fun_approx(θ, para)
#     @unpack τ, δ, A, η_L = para
#     f = jf(θ, A, η_L)
#     #ϵ_e_θ = τ*(θ+(1-η_L)*(1-θ))/(τ*(1+θ)+(1-δ)*f)
#     u = 1 - L_fun(θ, para)
#     ϵ_e_θ = u*(2*τ+(1-δ)f)*(θ+(1-η_L)*(1-θ))/(τ*(1+θ)+(1-δ)*f)
#     return ϵ_e_θ
# end

fig, ax = plt.subplots()
ax.plot(θ_grid, ϵ_e_fun.(θ_grid, Ref(para)), label="ϵ_e_θ")
#ax.plot(θ_grid, ϵ_e_fun_approx.(θ_grid, Ref(para)), label="approximation")
ax.plot(θ_grid,  1 .-L_fun.(θ_grid, Ref(para)), label="u")
ax.set_xlabel("θ")
ax.set_ylabel("e")
ax.legend()
display(fig)


###############
