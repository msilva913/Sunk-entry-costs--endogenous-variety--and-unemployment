
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

function inverse_gamma_map(m::Float64, s::Float64)
    if m <= 0 || s <= 0
        throw(ArgumentError("Mean and standard deviation must be positive."))
    end
# Calculate alpha using the mean formula
alpha = 2 + (m^2 / s^2)

# Calculate beta using the formula derived from the mean
beta = m * (alpha - 1)
return alpha, beta
end 



# sigma,      gamma_pdf,      1.5, 0.5;
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

# sigma prior 
α, β = gamma_map(1.5, 0.5)
gamma_dist = Gamma(α, β)
sigma_x = 0.5:0.1:2.0
sigma_prior_pdf = pdf(gamma_dist, sigma_x)

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

#b prior 
α, β = beta_map(0.71, 0.2)
beta_dist = Beta(α, β)
b_x = 0.6:0.005:0.96
b_prior_pdf = pdf(beta_dist, b_x)

# Shocks 


# Table: prior mean, prior std, posterior mean, posterior std
sigma_vals, sigma_density = columns(struc["sigma"])
ξ_inv_vals, ξ_inv_density = columns(struc["xi_inv"])
epsi_vals, epsi_density = columns(struc["epsi"])
delta_vals, delta_density = columns(struc["delta"])
b_vals, b_density = columns(struc["b_ratio"])

# Density of shock parameters
rho_z_vals, rho_z_density = columns(struc["rho_z"])
rho_delta_vals, rho_delta_density = columns(struc["rho_delta"])
rho_s_vals, rho_s_density = columns(struc["rho_s"])
sigma_z_vals, sigma_z_density = columns(struc["sigma_z"])
sigma_delta_vals, sigma_delta_density = columns(struc["sigma_delta"])
sigma_s_vals, sigma_s_density = columns(struc["sigma_delta"])

# Vector of dictionaries
params = [
    Dict("vals" => sigma_vals, "density" => sigma_density, "x" => sigma_x, "prior_pdf"=> sigma_prior_pdf, "xlabel"=> "σ"),
    Dict("vals" => ξ_inv_vals, "density" => ξ_inv_density, "x" => ξ_inv_x, "prior_pdf" => ξ_prior_pdf, "xlabel" => "ξ_inv"),
    Dict("vals" => epsi_vals, "density" => epsi_density, "x" => epsi_x, "prior_pdf" => epsi_prior_pdf, "xlabel" => "ε"),
    Dict("vals" => delta_vals, "density" => delta_density, "x" => delta_x, "prior_pdf" => delta_prior_pdf, "xlabel" => "δ"),
    Dict("vals" => b_vals, "density" => b_density, "x" => b_x, "prior_pdf" => b_prior_pdf, "xlabel" => "b")
]

fig, axs = subplots(2, 3, figsize=(18, 6))
axs = axs[:]

for i in 1:length(params)
    param = params[i]
    ax = axs[i]
    ax.plot(param["vals"], param["density"], linewidth=1.5, color="orange", label="Posterior", zorder=2)
    ax.fill_between(param["vals"], param["density"], color="gold", alpha=0.3, zorder=1)
    ax.plot(param["x"], param["prior_pdf"], linewidth=1, linestyle="--", color="blue", label="Prior", zorder=2)
    ax.fill_between(param["x"], param["prior_pdf"], color="cyan", alpha=0.3, zorder=1)
    ax.set_xlabel(param["xlabel"], fontsize=14)
    ax.set_ylabel("Density", fontsize=12)
    ax.legend(loc="upper right", fontsize=10)
    ax.grid(linestyle="--", alpha=0.7)
    ax.set_xlim(minimum(param["x"]), maximum(param["x"]))
    ax.set_ylim(0, maximum(param["density"]) * 1.1)
end
fig.delaxes(axs[end])
tight_layout()
display(fig)
plt.savefig("posterior_prior_plots.pdf")


# Shock processes
#rho_z
α, β = beta_map(0.8, 0.1)
beta_dist = Beta(α, β)
rho_z_x = 0.7:0.001:0.999
rho_z_prior_pdf = pdf(beta_dist, rho_z_x)

#rho_delta
α, β = beta_map(0.93, 0.05)
beta_dist = Beta(α, β)
rho_delta_x = 0.8:0.01:0.99
rho_delta_prior_pdf = pdf(beta_dist, rho_delta_x)

# rho_s
α, β = beta_map(0.8, 0.1)
beta_dist = Beta(α, β)
rho_s_x = 0.7:0.01:0.99
rho_s_prior_pdf = pdf(beta_dist, rho_s_x)

# sigma_z prior 
α, β = inverse_gamma_map(0.01, 1.0)
inverse_gamma_dist = InverseGamma(α, β)
sigma_z_x = 0.001:0.001:0.06
sigma_z_prior_pdf = pdf(inverse_gamma_dist, sigma_z_x)

# sigma_delta prior 
α, β = inverse_gamma_map(0.02, 0.01)
inverse_gamma_dist = InverseGamma(α, β)
sigma_delta_x = 0.001:0.005:0.06
sigma_delta_prior_pdf = pdf(inverse_gamma_dist, sigma_delta_x)

α, β = inverse_gamma_map(0.01, 1.0)
inverse_gamma_dist = InverseGamma(α, β)
sigma_s_x = 0.001:0.001:0.06
sigma_s_prior_pdf = pdf(inverse_gamma_dist, sigma_s_x)


params = [
    Dict("vals" => rho_z_vals, "density" => rho_z_density, "x" => rho_z_x, "prior_pdf" => rho_z_prior_pdf, "xlabel" => " ρ_z"),
    Dict("vals" => rho_delta_vals, "density" => rho_delta_density, "x" => rho_delta_x, "prior_pdf" => rho_delta_prior_pdf, "xlabel" => "ρ_δ"),
    Dict("vals" => rho_s_vals, "density" => rho_s_density, "x" => rho_s_x, "prior_pdf" => rho_s_prior_pdf, "xlabel" => "ρ_s"),
    Dict("vals" => sigma_z_vals, "density" => sigma_z_density, "x" => sigma_z_x, "prior_pdf" => sigma_z_prior_pdf, "xlabel" => " σ_z"),
    Dict("vals" => sigma_delta_vals, "density" => sigma_delta_density, "x" => sigma_delta_x, "prior_pdf" => sigma_delta_prior_pdf, "xlabel" => "σ_δ"),
    Dict("vals" => sigma_s_vals, "density" => sigma_s_density, "x" => sigma_s_x, "prior_pdf" => sigma_s_prior_pdf, "xlabel" => "σ_s")
]

fig, axs = subplots(2, 3, figsize=(14, 4))
axs = [axs[i, j] for j in 1:size(axs, 2), i in 1:size(axs, 1)]
for (i, ax) in enumerate(axs)
    param = params[i]
    ax.plot(param["vals"], param["density"], linewidth=1.5, color="orange", label="Posterior", zorder=2)
    ax.fill_between(param["vals"], param["density"], color="gold", alpha=0.3, zorder=1)
    ax.plot(param["x"], param["prior_pdf"], linewidth=1, linestyle="--", color="blue", label="Prior", zorder=2)
    ax.fill_between(param["x"], param["prior_pdf"], color="cyan", alpha=0.3, zorder=1)
    ax.set_xlabel(param["xlabel"], fontsize=14)
    ax.set_ylabel("Density", fontsize=12)
    ax.legend(loc="upper right", fontsize=10)
    ax.grid(linestyle="--", alpha=0.7)
    ax.set_xlim(minimum(param["x"]), maximum(param["x"]))
    ax.set_ylim(0, maximum(param["density"]) * 1.1)
end

tight_layout()
display(fig)
savefig("posterior_priors_shocks.pdf")






