
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
#cd(@__DIR__)
cd("C:/Users/msilva913/Documents/GitHub/Sunk_entry_costs_endogenous_variety_unemployment")
columns(M) = (view(M, :, i) for i in 1:size(M, 2))

posterior = matopen("posterior_density.mat")
struc = read(posterior, "posterior_density")

print(keys(struc))
n = length(keys(struc))

# Priors 
priors = matopen("priors.mat")
priors = read(priors, "priors")
priors_mean = priors["mean"]

# Posterior mean
posterior_mean = matopen("posterior_mean.mat")
posterior_mean_struc = read(posterior_mean, "posterior_mean")

# Posterior mode 
posterior_mode = matopen("posterior_mode.mat")
posterior_mode_struc = read(posterior_mode, "posterior_mode")


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


"""
sigma,      gamma_pdf,      1.5, 0.5;
b_ratio,    beta_pdf,       0.7, 0.1;
xi_inv,     gamma_pdf,      2.0, 2.0;
x_v,        beta_pdf,       0.5, 0.2;
epsi,       gamma_pdf,      4.2, 1.5;
rho_z,      beta_pdf,       0.8, 0.10;
rho_s,      beta_pdf,       0.6, 0.20;
rho_delta,  beta_pdf,       0.9, 0.05;
sigma_z,    inv_gamma_pdf,  0.01, 0.01;
sigma_s,    inv_gamma_pdf,  0.06, 0.06;
sigma_delta,inv_gamma_pdf,  0.02, 0.02;
"""

# sigma prior 
α, β = gamma_map(1.5, 0.5)
gamma_dist = Gamma(α, β)
sigma_x = 0.5:0.1:2.0
sigma_prior_pdf = pdf(gamma_dist, sigma_x)

# xi_inv prior 
α, β = gamma_map(2.0, 2.0)
gamma_dist = Gamma(α, β)
ξ_inv_x = 0:0.05:1.0
ξ_prior_pdf = pdf(gamma_dist, ξ_inv_x)

# epsi prior
α, β = gamma_map(4.2, 1.5)
gamma_dist = Gamma(α, β)
epsi_x = 2:0.1:10
epsi_prior_pdf = pdf(gamma_dist, epsi_x)

#δ prior

# α, β = beta_map(0.0083, 0.005)
# beta_dist = Beta(α, β)
# delta_x = 0.006:0.00025:0.025
# delta_prior_pdf = pdf(beta_dist, delta_x)

#b prior 
α, β = beta_map(0.7, 0.2)
beta_dist = Beta(α, β)
b_x = 0.6:0.005:0.96
b_prior_pdf = pdf(beta_dist, b_x)

# Shocks 


# Table: prior mean, prior std, posterior mean, posterior std
sigma_vals, sigma_density = columns(struc["sigma"])
ξ_inv_vals, ξ_inv_density = columns(struc["xi_inv"])
epsi_vals, epsi_density = columns(struc["epsi"])
#delta_vals, delta_density = columns(struc["delta"])
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
    Dict("vals" => b_vals, "density" => b_density, "x" => b_x, "prior_pdf" => b_prior_pdf, "xlabel" => "b")
    #Dict("vals" => delta_vals, "density" => delta_density, "x" => delta_x, "prior_pdf" => delta_prior_pdf, "xlabel" => "δ"),
   
]

fig, axs = subplots(2, 2, figsize=(10, 6))
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
#fig.delaxes(axs[end])
tight_layout()
display(fig)
plt.savefig("posterior_prior_plots_structural.pdf")


# Shock processes
#rho_z
α, β = beta_map(0.8, 0.1)
beta_dist = Beta(α, β)
rho_z_x = 0.7:0.001:0.999
rho_z_prior_pdf = pdf(beta_dist, rho_z_x)

#rho_delta
α, β = beta_map(0.9, 0.05)
beta_dist = Beta(α, β)
rho_delta_x = 0.8:0.01:0.99
rho_delta_prior_pdf = pdf(beta_dist, rho_delta_x)

# rho_s
α, β = beta_map(0.6, 0.2)
beta_dist = Beta(α, β)
rho_s_x = 0.7:0.01:0.99
rho_s_prior_pdf = pdf(beta_dist, rho_s_x)

# sigma_z prior 
α, β = inverse_gamma_map(0.01, 0.01)
inverse_gamma_dist = InverseGamma(α, β)
sigma_z_x = 0.001:0.001:0.03
sigma_z_prior_pdf = pdf(inverse_gamma_dist, sigma_z_x)

# sigma_delta prior 
α, β = inverse_gamma_map(0.02, 0.02)
inverse_gamma_dist = InverseGamma(α, β)
sigma_delta_x = 0.001:0.005:0.03
sigma_delta_prior_pdf = pdf(inverse_gamma_dist, sigma_delta_x)

# sigma_s prior
α, β = inverse_gamma_map(0.06, 0.06)
inverse_gamma_dist = InverseGamma(α, β)
sigma_s_x = 0.001:0.001:0.1
sigma_s_prior_pdf = pdf(inverse_gamma_dist, sigma_s_x)


params = [
    Dict("vals" => rho_z_vals, "density" => rho_z_density, "x" => rho_z_x, "prior_pdf" => rho_z_prior_pdf, "xlabel" => " ρ_z"),
    Dict("vals" => rho_delta_vals, "density" => rho_delta_density, "x" => rho_delta_x, "prior_pdf" => rho_delta_prior_pdf, "xlabel" => "ρ_δ"),
    Dict("vals" => rho_s_vals, "density" => rho_s_density, "x" => rho_s_x, "prior_pdf" => rho_s_prior_pdf, "xlabel" => "ρ_s"),
    Dict("vals" => sigma_z_vals, "density" => sigma_z_density, "x" => sigma_z_x, "prior_pdf" => sigma_z_prior_pdf, "xlabel" => " σ_z"),
    Dict("vals" => sigma_delta_vals, "density" => sigma_delta_density, "x" => sigma_delta_x, "prior_pdf" => sigma_delta_prior_pdf, "xlabel" => "σ_δ"),
    Dict("vals" => sigma_s_vals, "density" => sigma_s_density, "x" => sigma_s_x, "prior_pdf" => sigma_s_prior_pdf, "xlabel" => "σ_s")
]

fig, axs = subplots(2, 3, figsize=(12, 5))
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
savefig("posterior_prior_plots_shocks.pdf")

####### Posterior-prior tables ########
posterior_HPDlow = matopen("posterior_HPDlow.mat")
posterior_HPDlow = read(posterior_HPDlow, "posterior_HPDlow")

posterior_HPDhigh = matopen("posterior_HPDhigh.mat")
posterior_HPDhigh = read(posterior_HPDhigh, "posterior_HPDhigh")

posterior_std = matopen("posterior_std.mat")
posterior_std = read(posterior_std, "posterior_std")


prior_dict = Dict(
    "sigma"       => ("Gamma", 1.5, 0.5),
    "b_ratio"     => ("Beta", 0.7, 0.1),
    "xi_inv"      => ("Gamma", 2.0, 2.0),
    "x_v"         => ("Beta", 0.5, 0.2),
    "epsi"        => ("Gamma", 4.2, 1.5),
    "rho_z"       => ("Beta", 0.8, 0.10),
    "rho_s"       => ("Beta", 0.6, 0.20),
    "rho_delta"   => ("Beta", 0.9, 0.05),
    "sigma_z"     => ("Inverse Gamma", 0.01, 0.01),
    "sigma_s"     => ("Inverse Gamma", 0.06, 0.06),
    "sigma_delta" => ("Inverse Gamma", 0.02, 0.02)
)


print(keys(posterior_mean_struc))
keys_r = ["sigma", "b_ratio", "xi_inv", "x_v", "epsi", "rho_z", "rho_s", "rho_delta", "sigma_z", "sigma_s", "sigma_delta"]

posterior_table = DataFrame(
    Parameter = keys_r,
    PriorDistribution = [prior_dict[param][1] for param in keys_r],
    PriorMean =  [prior_dict[param][2] for param in keys_r],
    PriorStd = [prior_dict[param][3] for param in keys_r],
    PosteriorMean = [posterior_mean_struc[param] for param in keys_r],
    PosteriorStd = [posterior_std[param] for param in keys_r],
    PosteriorMode = [posterior_mode_struc[param] for param in keys_r],
    LowerHPD = [posterior_HPDlow[param] for param in keys_r],
    HigherHPD =[posterior_HPDhigh[param] for param in keys_r]
)

for col in names(posterior_table)
    if eltype(posterior_table[!, col]) <: AbstractFloat  # Assuming the first column is non-numeric (like 'Parameter')
    posterior_table[!, col] = round.(posterior_table[!, col], sigdigits=3)
    end
end

print(posterior_table)


function dataframe_to_latex_multicolumn(df::DataFrame)
    header = [
        "\\begin{table}[ht]",
        "\\centering",
        "\\begin{tabular}{c|ccc|ccccc}",
        "\\toprule",
        "& \\multicolumn{3}{c|}{Priors} & \\multicolumn{5}{c|}{Posteriors} \\\\ \\midrule",
        "Parameter & Distribution & Mean & Std & Mean & Mode & Std & Lower HPD & Higher HPD \\\\"
    ]

    rows = [
        string(
            row.Parameter, " & ", row.PriorDistribution, " & ",
            round(row.PriorMean, sigdigits=2), " & ", round(row.PriorStd, sigdigits=2), " & ",
            round(row.PosteriorMean, sigdigits=2), " & ", round(row.PosteriorMode, sigdigits=2), " & ",
            round(row.PosteriorStd, sigdigits=2), " & ", round(row.LowerHPD, sigdigits=2), " & ", round(row.HigherHPD, sigdigits=2),
            " \\\\"
        ) for row in eachrow(df)
    ]

    footer = [
        "\\bottomrule",
        "\\end{tabular}",
        "\\caption{Posterior estimates for model parameters with prior information.}",
        "\\label{tab:posterior_estimates}",
        "\\end{table}"
    ]

    return join([header; rows; footer], "\n")
end

latex_table = dataframe_to_latex_multicolumn(posterior_table)
println(latex_table)