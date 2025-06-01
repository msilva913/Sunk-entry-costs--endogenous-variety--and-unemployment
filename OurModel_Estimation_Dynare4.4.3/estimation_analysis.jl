
using CSV, ArgParse
using Latexify
#using LinearAlgebra, Roots, Optim, NLsolve
using Printf
using DataFrames
using PyPlot
using MAT
using TexTables
using Dates, Distributions
#using StatsModels

homedir()
#include("time_series_fun.jl")
cd(@__DIR__)
columns(M) = (view(M, :, i) for i in 1:size(M, 2))

posterior = matopen("posterior_density.mat")
struc = read(posterior, "posterior_density")

print(keys(struc))
n = length(keys(struc))

# Posterior mean
posterior_mean = matopen("posterior_mean.mat")
posterior_mean_struc = read(posterior_mean, "posterior_mean")

posterior_mean_struc_shocks = read(posterior_mean_shocks, "posterior_mean_shocks")

# Posterior standard deviation
# posterior_std = matopen("posterior_std.mat")
# posterior_std_shocks = matopen("posterior_std_shocks.mat")
# posterior_std_struc = read(posterior_std, "posterior_std")
# posterior_std_struc_shocks = read(posterior_std_shocks, "posterior_std_shocks")

# Extract keys


function beta_map(m::Float64, s::Float64)
    if m <= 0 || m >= 1 || s < 0 || s >= 1
        error("Mean must be in (0, 1) and StdDev must be non-negative and less than 1.")
    end
    α = (m * (1 - m) / s^2 - 1) * m
    β = (m * (1 - m) / s^2 - 1) * (1 - m)
    return α, β
end

function gamma_map(m::Float64, s::Float64)
    if m <= 0 || s <= 0
        throw(ArgumentError("Mean and standard deviation must be positive."))
    end

    beta = s^2 / m
    alpha = m^2 / s^2
    return alpha, beta
end

# b_ratio,    beta_pdf,       0.71, 0.2;
# x_v,        beta_pdf,       0.5, 0.25;
# xi_inv,     gamma_pdf,      1.0, 2;
# delta,      beta_pdf,       0.0083, 0.005; 0.0083 
# epsi,       gamma_pdf,      4.2, 1.5;
# rho_z,      beta_pdf,       0.8, 0.2;
# rho_delta,  beta_pdf,       0.93, 0.1;
# rho_s,      beta_pdf,       0.8, 0.2;
# sigma_z, 0.01, 0.000001, 0.2,    inv_gamma_pdf,  0.01, 1.0;
# sigma_delta, 0.01, 0.00001, 0.2,     inv_gamma_pdf,  0.02, 0.01; 
# sigma_s, 0.01, 0.00001, 0.2,            inv_gamma_pdf,  0.01, 1.0;

# xi_inv prior 
α, β = gamma_map(1.0, 2.0)
gamma_dist = Gamma(α, β)
ξ_inv_x = 0:0.1:6.0
ξ_prior_pdf = pdf(gamma_dist, ξ_inv_x)

# epsi prior
α, β = gamma_map(4.2, 1.5)
gamma_dist = Gamma(α, β)
epsi_x = 2:0.1:6.0
epsi_prior_pdf = pdf(gamma_dist, epsi_x)

#δ prior
α, β = beta_map(0.0083, 0.005)
beta_dist = Beta(α, β)
delta_x = 0.006:0.00025:0.025
delta_prior_pdf = pdf(beta_dist, delta_x)


# Table: prior mean, prior std, posterior mean, posterior std
ξ_inv_vals, ξ_inv_density = columns(struc["xi_inv"])
epsi_vals, epsi_density = columns(struc["epsi"])
delta_vals, delta_density = columns(struc["delta"])


fig = plt.figure(figsize=(14, 4))
# First subplot for xi_inv
ax1 = fig.add_subplot(1, 3, 1)
ax1.plot(ξ_inv_vals, ξ_inv_density, linewidth=1.5, color="orange", label="Posterior", zorder=2)
ax1.fill_between(ξ_inv_vals, ξ_inv_density, color="gold", alpha=0.3, zorder=1)
ax1.plot(ξ_inv_x, ξ_prior_pdf, linewidth=1, linestyle="--", color="blue", label="Prior", zorder=2)
ax1.fill_between(ξ_inv_x, ξ_prior_pdf, color="cyan", alpha=0.3, zorder=1)
ax1.set_xlabel("ξ_inv", fontsize=14)
ax1.set_ylabel("Density", fontsize=12)
ax1.legend(loc="upper right", fontsize=10)
ax1.grid(linestyle="--", alpha=0.7)
ax1.set_xlim(minimum(ξ_inv_x), maximum(ξ_inv_x))
ax1.set_ylim(0, max(maximum(ξ_inv_density), maximum(ξ_inv_density)) * 1.1)

# Second subplot for ε
ax2 = fig.add_subplot(1, 3, 2)
ax2.plot(epsi_vals, epsi_density, linewidth=1.5, color="orange", label="Posterior", zorder=2)
ax2.fill_between(epsi_vals, epsi_density, color="gold", alpha=0.3, zorder=1)
ax2.plot(epsi_x, epsi_prior_pdf, linewidth=1, linestyle="--", color="blue", label="Prior", zorder=2)
ax2.fill_between(epsi_x, epsi_prior_pdf, color="cyan", alpha=0.3, zorder=1)
ax2.set_xlabel("ε", fontsize=14)
ax2.set_ylabel("Density", fontsize=12)
ax2.legend(loc="upper right", fontsize=10)
ax2.grid(linestyle="--", alpha=0.7)
ax2.set_xlim(minimum(epsi_x), maximum(epsi_x))
ax2.set_ylim(0, max(maximum(epsi_density), maximum(epsi_density)) * 1.1)

# Third subplot for δ
ax3 = fig.add_subplot(1, 3, 3)
ax3.plot(delta_vals, delta_density, linewidth=1.5, color="orange", label="Posterior", zorder=2)
ax3.fill_between(delta_vals, delta_density, color="gold", alpha=0.3, zorder=1)
ax3.plot(delta_x, delta_prior_pdf, linewidth=1, linestyle="--", color="blue", label="Prior", zorder=2)
ax3.fill_between(delta_x, delta_prior_pdf, color="cyan", alpha=0.3, zorder=1)
ax3.set_xlabel("δ", fontsize=14)
ax3.set_ylabel("Density", fontsize=12)
ax3.legend(loc="upper right", fontsize=10)
ax3.grid(linestyle="--", alpha=0.7)
ax3.set_xlim(minimum(delta_x), maximum(delta_x))
ax3.set_ylim(0, max(maximum(delta_density), maximum(delta_density)) * 1.1)

plt.tight_layout()
display(fig)
plt.savefig("posterior_prior_plots.pdf")


# phi prior
α, β = beta_map(0.32, 0.2)
beta_dist = Beta(α, β)
ϕ_x = 0:0.01:1
ϕ_prior_pdf = pdf(beta_dist, ϕ_x)

# eta prior
α, β = gamma_map(0.2, 0.15)
gamma_dist = Gamma(α, β)
η_x = 0:0.01:1
gamma_prior_pdf = pdf(beta_dist, η_x)

# νR_prior
α, β = beta_map(0.2, 0.1)
beta_dist = Beta(α, β)
νR_x = 0:0.01:1
νR_prior_pdf = pdf(beta_dist, νR_x )

# Distribution: structural parameters
key_map = ["σ_a", "ζ", "η", "ρ_ZI", "ρ_N", "ρ_D", "θ", "Ψ_K", "ρ_C", "ρ_g"]



fig = plt.figure(figsize=(14, 4))
# First subplot for ϕ
ax2 = fig.add_subplot(1, 3, 1)
ax2.plot(ϕ_vals, ϕ_density, linewidth=1.5, color="orange", label="Posterior", zorder=2)
ax2.fill_between(ϕ_vals, ϕ_density, color="gold", alpha=0.3, zorder=1)
ax2.plot(ϕ_x, ϕ_prior_pdf, linewidth=1, linestyle="--", color="blue", label="Prior", zorder=2)
ax2.fill_between(ϕ_x, ϕ_prior_pdf, color="cyan", alpha=0.3, zorder=1)
ax2.set_xlabel("ϕ", fontsize=14)
ax2.set_ylabel("Density", fontsize=12)
ax2.legend(loc="upper right", fontsize=10)
ax2.grid(linestyle="--", alpha=0.7)
ax2.set_xlim(minimum(ϕ_x), maximum(ϕ_x))
ax2.set_ylim(0, max(maximum(ϕ_density), maximum(ϕ_prior_pdf)) * 1.1)

# Second subplot for η
ax1 = fig.add_subplot(1, 3, 2)
ax1.plot(η_vals, η_density, linewidth=1.5, color="orange", label="Posterior", zorder=2)
ax1.fill_between(η_vals, η_density, color="gold", alpha=0.3, zorder=1)
ax1.plot(η_x, gamma_prior_pdf, linewidth=1, linestyle="--", color="blue", label="Prior", zorder=2)
ax1.fill_between(η_x, gamma_prior_pdf, color="cyan", alpha=0.3, zorder=1)
ax1.set_xlabel("η", fontsize=12)
ax1.set_ylabel("Density", fontsize=14)
ax1.legend(loc="upper right", fontsize=10)
ax1.grid(linestyle="--", alpha=0.7)
ax1.set_xlim(minimum(η_x), maximum(η_x))
ax1.set_ylim(0, max(maximum(η_density), maximum(gamma_prior_pdf)) * 1.1)

ax1 = fig.add_subplot(1, 3, 3)
ax1.plot(νR_vals, νR_density, linewidth=1.5, color="orange", label="Posterior", zorder=2)
ax1.fill_between(νR_vals, νR_density, color="gold", alpha=0.3, zorder=1)
ax1.plot(νR_x, νR_prior_pdf, linewidth=1, linestyle="--", color="blue", label="Prior", zorder=2)
ax1.fill_between(νR_x, νR_prior_pdf, color="cyan", alpha=0.3, zorder=1)
ax1.set_xlabel("νR", fontsize=14)
ax1.set_ylabel("Density", fontsize=12)
ax1.legend(loc="upper right", fontsize=10)
ax1.grid(linestyle="--", alpha=0.7)
ax1.set_xlim(minimum(νR_x), maximum(νR_x))
ax1.set_ylim(0, max(maximum(νR_density), maximum(νR_prior_pdf)) * 1.1)

# Adjust layout and display
plt.tight_layout()
display(fig)
plt.savefig("posterior_prior_plots.pdf")


function cumulate(x)
    irf_levels = 100*cumsum(x)
    irf_levels =  cat(zeros(1), irf_levels, dims=1)
end

extract_series(str, dic::Dict) = vec(dic[str])

function irf_fun(vars, irf_dic; shock="e_D", length=20)
    periods = 1:length
    irf_array = []
    for (i, key) in enumerate(vars)
        str = key*"_"*shock
        irf = irf_dic[str]
        push!(irf_array, irf)
    end
    return irf_array
end

function irf_fun_plot(irf_array, vars_list, vars_list_label; shock, savefig=true)
    fig = plt.figure(figsize=(16, 10))
    periods = 1:length(irf_array[1][:])
    periods = 1:length(irf_array[1])
    for (i, key) in enumerate(vars_list)
        irf = transpose(irf_array[i])
        ax = fig.add_subplot(3, 3, i)
        ax.plot(periods, 100*irf, linewidth=3, color="teal")
        ax.tick_params(labelsize=12)
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
        ax.set_title(vars_list_label[i])
        plt.tight_layout()
        fig.suptitle("A 1 standard-deviation shock to "*shock, fontsize=14)
    end 
    display(fig)
    if savefig
        figname = "irf_"*shock*".pdf"
    end
    plt.savefig(figname)
end

function irf_fun_plot_grouped(irf_dic; shock, savefig=true)
    fig = plt.figure(figsize=(16, 8))
    periods = 1:length(irf_array[1][:])
    periods = 1:length(irf_array[1])
    # lists of grouped variables
    list_1 = [:C, :I]
    list_2 = [:N_C, :N_I]
    list_3 = [:util_ND, :util_D, :util]
    lists = (list_1, list_2, list_3, [:SR], [:D, :h, :util], [:p_I])
    linestyles = ["solid", "dashed", "dotted"]
    j = 1
    for (n, list) in enumerate(lists)
        ax = fig.add_subplot(2, 3, n)
        ax.tick_params(labelsize=12)
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
        for (i, key) in enumerate(list)
            irf = transpose(irf_dic[key])
            ax.plot(periods, 100*irf, linewidth=3, label=key, linestyle = linestyles[i], alpha=0.7)
            ax.legend()
            #ax.set_title(vars_list_label[i])
            j += 1
        end 
    end
    fig.suptitle("A 1 standard-deviation shock to "*shock, fontsize=14)
    plt.tight_layout()
    display(fig)
    if savefig
        figname = "irf_"*shock*".pdf"
    end
    plt.savefig(figname)
end

# Standard impulse responses: basic model
#vars_list = ["log_C", "log_I", "log_Y_N", "log_p_I", "log_NC", "log_NI"]
#vars_list_label=[:C, :I, :Y_N, :p_I, :N_C, :N_I]
vars_list = ["C_obs", "I_obs", "NC_obs", "NI_obs", "util_ND_obs", "util_D_obs", "util_obs", "SR_obs", "D_obs", "h_obs", "p_I_obs"]
vars_list_label = [:C, :I, :N_C, :N_I, :util_ND, :util_D, :util, :SR, :D, :h, :p_I]

x = matopen("irf.mat")
vars = read(x, "irf")
irf_dic  = vars



# Shopping preference shock
irf_array = irf_fun(vars_list, irf_dic, shock="e_D", length=20 )
irf_dic_spec = Dict(zip(vars_list_label, irf_array))
irf_fun_plot_grouped(irf_dic_spec; shock="e_D")

# Stationary technology shock
irf_array = irf_fun(vars_list, irf_dic, shock="e_Z", length=20 )
irf_dic_spec = Dict(zip(vars_list_label, irf_array))
irf_fun_plot_grouped(irf_dic_spec; shock="e_Z")

# Discount factor shock
irf_array = irf_fun(vars_list, irf_dic, shock="e_b", length=20 )
irf_dic_spec = Dict(zip(vars_list_label, irf_array))
irf_fun_plot_grouped(irf_dic_spec; shock="e_b")

# Permanent technology shock
irf_array = irf_fun(vars_list, irf_dic, shock="e_g", length=20 )
irf_dic_spec = Dict(zip(vars_list_label, irf_array))
irf_fun_plot_grouped(irf_dic_spec; shock="e_g")



# In levels
irf_array = irf_fun(vars_list, irf_dic, shock="e_b", length=20, cumulate_resp=true )
irf_fun_plot(irf_array, vars_list, vars_list_label; shock="e_b", savefig=false)

# I-specific shopping shock
irf_array = irf_fun(vars_list, irf_dic, shock="e_DI", length=20 )
irf_fun_plot(irf_array, vars_list, vars_list_label; shock="e_DI")
# Shopping preference shock


# Wage markup shock

irf_array = irf_fun(vars_list, irf_dic, shock="e_muC", length=20 )
irf_fun_plot(irf_array, vars_list, vars_list_label; shock="e_muC")

irf_array = irf_fun(vars_list, irf_dic, shock="e_muI", length=20 )
irf_fun_plot(irf_array, vars_list, vars_list_label; shock="e_muI")

## Smoothed variables
sv = matopen("sv.mat")
sv = read(sv, "sv")
pprint(keys(sv))


# 1. Subplot of shopping effort and utilization series
# 2. Subplot decomposing utilization
# 3. Subplot expressing utilization and Solow residual in growth rates

# Create datetime object 
# Create a range of dates starting from 1964 Q1
start_date = Date(1964, 1, 1)
end_date = start_date + Dates.Quarter(222)  # 223 elements long
date_range = collect(start_date:Dates.Quarter(1):end_date)

D, D_obs, util, util_obs, C_obs, tech_obs, SR_obs, util_ND, util_D, util_sc = 
    [100 * vec(sv[symbol]) for symbol in ["D", "D_obs", "util", "util_obs", "C_obs", "tech_obs", 
    "SR_obs", "util_ND", "util_D", "util_sc"]]

SR_cum, tech_cum, C_cum = [100 * cumsum(vec(sv[symbol])) for symbol in ["SR_obs", "tech_obs", "C_obs"]]

#
var(SR_obs)
var(tech_obs)/var(SR_obs)
df = DataFrame()
df[!, :C_obs] = C_obs 
df[!, :tech_obs] = tech_obs
df[!, :C_obs_lag1] = lag(df.C_obs, 1)
df[!, :tech_obs_lag1] = lag(df.tech_obs, 1)
df[!, :tech_obs_lag2] = lag(df.tech_obs, 2)
df[!, :tech_obs_lag3] = lag(df.tech_obs, 3)
df[!, :tech_obs_lag4] = lag(df.tech_obs, 4)
model = lm(@formula(tech_obs ~ tech_obs_lag1 + tech_obs_lag2+tech_obs_lag3+tech_obs_lag4+ C_obs_lag1), df)

# Annual percentage changes in shopping effort 
function quarter_perc_ann(x)
    n = length(x)
    last_full_set = n ÷ 4
    remaining_quarters = n % 4
    ann_changes_full_sets = [prod(1 .+ x[(4*(i-1)+1):(4*i)]) - 1 for i in 1:last_full_set]
    
    last_year_changes = prod(1 .+ x[(4*last_full_set+1):end]) - 1
    ann_changes = vcat(ann_changes_full_sets, last_year_changes)
    return ann_changes
end

function to_annual(quarterly_date_range)
    annual_dates = [quarterly_date_range[i] for i in 1:4:length(quarterly_date_range)]
    return Dates.Date.(annual_dates)
end

D_obs_ann = quarter_perc_ann(D_obs)
C_obs_ann = quarter_perc_ann(C_obs)
tech_obs_ann = quarter_perc_ann(tech_obs)
SR_obs_ann = quarter_perc_ann(SR_obs)

# Filtering data
dates_ann = to_annual(date_range)
dates_red = filter(date -> year(date) in 2003:2019, dates_ann)
indices = findall(date -> year(date) in 2003:2019, dates_red)


# Annualized shopping effort percentage changes
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(dates_red, D_obs_ann[indices])
display(fig)

fig = plt.figure(figsize=(16, 4))
# Shopping effort and utilization
ax1 = fig.add_subplot(1, 3, 1)
ax1.plot(date_range, D.*util[1]./D[1], linewidth=1.5, color="orange", label="Shopping effort", zorder=2, alpha=0.7)
ax1.plot(date_range, util, linewidth=1.5, color="blue", label="Utilization", alpha=0.7)
ax1.set_xlabel("Time", fontsize=12)
ax1.set_ylabel("Units", fontsize=12)
ax1.legend(loc="upper right", fontsize=10)
ax1.grid(linestyle="--", alpha=0.7)


# Consumption, TFP, and Technology
ax2 = fig.add_subplot(1, 3, 2)
ax2.plot(dates_ann, C_obs_ann, linewidth=1.5, color="orange", label="Consumption", zorder=2, alpha=0.7)
ax2.plot(dates_ann, SR_obs_ann, linewidth=1.5, color="blue", label="Solow residual", alpha=0.7)
ax2.plot(dates_ann, tech_obs_ann, linewidth=1.5, color="red", label="Technology (purified Solow residual)", alpha=0.7)
ax2.set_xlabel("Time", fontsize=12)
#ax2.set_ylabel("Units", fontsize=12)
#ax2.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
ax2.legend(loc="upper right", fontsize=10)
ax2.grid(linestyle="--", alpha=0.7)

# Utilization series
ax3 = fig.add_subplot(1, 3, 3)
ax3.plot(date_range, util_ND, linewidth=1.5, color="orange", label="Utilization: non-durables", zorder=2, alpha=0.7)
ax3.plot(date_range, util_D, linewidth=1.5, color="blue", label="Utilizaton: durables", alpha=0.7)
#ax3.plot(date_range, util_sc, linewidth=1.5, color="red", label="Utilizaton: services", alpha=0.7)
ax3.set_xlabel("Time", fontsize=12)
ax3.set_ylabel("Units", fontsize=12)
#ax2.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
ax3.legend(loc="upper right", fontsize=10)
ax3.grid(linestyle="--", alpha=0.7)
display(fig)














