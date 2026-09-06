
using LaTeXStrings, KernelDensity
using Parameters, CSV, StatsBase, Statistics, Random
using DataFrames
using ShiftedArrays
using MappedArrays
using TexTables
using TypedTables
using GLM
using LinearAlgebra
using SparseArrays

columns(M) = (view(M, :, i) for i in 1:size(M, 2))

function moments(dat, var_RSD, var_corr; lags =2, verbose=true)
    #sd = DataFrames.colwise(std, dat)
    # Drop missing values
    dat = dropmissing(dat)
    sd = map(std, eachcol(dat))
    RSD = sd./std(dat[!, var_RSD])

    corrs = zeros(ncol(dat), length(var_corr))
    lab = []
    for (i, var) in enumerate(var_corr)
        corrs[:, i] = [cor(dat[:, i], dat[:, var]) for i in 1:ncol(dat)]
        lab
    end
    cor_lab = ["Cor(x, $(var_corr[n]))" for n = 1:length(var_corr)]
    cor_lab = permutedims(cor_lab)
    
  
    ac = zeros(ncol(dat), lags)

    for k in 2:(lags+1)
        ac[:, k-1] .= [autocor(Float64.(dat[:, i]))[k] for i in 1:size(dat, 2)]
    end

    mom = [names(dat) sd RSD corrs ac]
    #mom = DataFrames.DataFrame(mom, :auto)
    mom = DataFrames.DataFrame(mom, :auto)
    name = ["Variable" "SD" "RSD" cor_lab "Cor(x, x_{-1})" "Cor(x, x_{-2})"]
    DataFrames.rename!(mom, vec(name))
    #Table(mom)
    if verbose
        mom[!, 2:end] = round.(mom[!, 2:end], sigdigits=3)
    end
    mom[!, :Variable] = names(dat)
    return mom
end


function hamilton_filter(x::Vector; h=8)
    x_h = ShiftedArrays.lag(x, h)
    x_h1 = ShiftedArrays.lag(x_h, 1)
    x_h2 = ShiftedArrays.lag(x_h, 2)
    x_h3 = ShiftedArrays.lag(x_h, 3)
    X = DataFrame(x=x, x_h=x_h, x_h1=x_h1, x_h2=x_h2, x_h3=x_h3)  # Construct DataFrame directly
    # rename
    ols = lm(@formula(x ~ x_h + x_h1 + x_h2 + x_h3), X)
    missing_values = Vector(undef, h+3)
    missing_values .= missing
    return vcat(missing_values, residuals(ols))
end


function growth_filter(x)
    x = diff(x)
    # demean
    x = x .- mean(x)
    return vcat(missing, x)
end

function linear_filter(x; intercept=true)
    X = DataFrame()
    X[!, :x] = x
    T = range(1, size(X)[1], step=1)
    #X[!, :ones] = ones(size(x))
    X[!, :time] = T
    ols = lm(@formula(x~  time), X)
    cons = coef(ols)[1]
    if intercept
        return residuals(ols)
    else
        return residuals(ols) .+ cons
    end
end  

function polynomial_filter(x, n=1; intercept=true)
    T = range(1, size(x)[1], step=1)
    Time_vars = zeros(length(T), n+1)
    for i in 1:(n+1)
        Time_vars[:,i] = T.^(i-1)
    end
    ols = lm(Time_vars, x)
    cons = coef(ols)[1]
    if intercept
        return residuals(ols)
    else
        return residuals(ols) .+ cons
    end
end  

function bk_filter(y::DataFrames.DataFrame; wl=6, wu=32, K=12)
    ### Arguments
    #y: data to be filtered
    #wu: upper cutoff frequencies
    #wl: lower cutoff frequencies
    #K: number of leads and lags of moving average

    # Returns
    #y_cycle: cyclical component
    #y_trend: trend component
    T = length(y)
    #T = nrow(y)
    w1 = 2pi/wu
    w2 = 2pi/wl
    b = vcat((w2-w1)/pi, [(sin(w2*i)-sin(w1*i))/(i*pi) for i=1:K])
    theta = (b[1] + 2*sum(b[2:end]))/(2K+1)
    B = b .- theta
    #vals = y.value
    #y[!, :cycle] = Vector{Union{Missing, Float64}}(undef, T)
    cycle = Vector{Union{Missing, Float64}}(undef, T)
    for t = K+1:T-K
        y[t,:cycle] = vals[1]*B[1] + (vals[t-1:-1:t-K]'*B[2:end]) + (vals[t+1:t+K]'*B[2:end])
    end
    y[!,:trend] = y[:, :value] - y[:, :cycle]
    return dropmissing(y)
end

function bkfilter(y::Vector{<:Union{Missing, Float64}}; wl=6, wu=32, K=12)
    ### Arguments
    #y: data to be filtered
    #wu: upper cutoff frequencies
    #wl: lower cutoff frequencies
    #K: number of leads and lags of moving average

    # Returns
    #y_cycle: cyclical component
    #y_trend: trend component
    T = length(y)
    w1 = 2pi/wu
    w2 = 2pi/wl
    b = vcat((w2-w1)/pi, [(sin(w2*i)-sin(w1*i))/(i*pi) for i=1:K])
    theta = (b[1] + 2*sum(b[2:end]))/(2K+1)
    B = b .- theta
    cycle = Vector{Union{Missing, Float64}}(undef, T)
    for t = K+1:T-K
        cycle[t] = y[1]*B[1] + (y[t-1:-1:t-K]'*B[2:end]) + (y[t+1:t+K]'*B[2:end])
    end
    #trend = y - cycle
    return filter!(!ismissing, cycle)
end

function bkfilter(y; wl=6, wu=32, K=12)
    ### Arguments
    #y: data to be filtered
    #wu: upper cutoff frequencies
    #wl: lower cutoff frequencies
    #K: number of leads and lags of moving average

    # Returns
    #y_cycle: cyclical component
    #y_trend: trend component
    T = length(y)
    w1 = 2pi/wu
    w2 = 2pi/wl
    b = vcat((w2-w1)/pi, [(sin(w2*i)-sin(w1*i))/(i*pi) for i=1:K])
    theta = (b[1] + 2*sum(b[2:end]))/(2K+1)
    B = b .- theta
    cycle = Vector{Union{Missing, Float64}}(undef, T)
    for t = K+1:T-K
        cycle[t] = y[1]*B[1] + (y[t-1:-1:t-K]'*B[2:end]) + (y[t+1:t+K]'*B[2:end])
    end
    #trend = y - cycle
    return filter!(!ismissing, cycle)
end

# function hp_filter(y, λ) where T<:Real
#     ### Arguments
#     #y: data to be filtered
#      #λ: smoothing parameter (6.25 for annual, 1600 for quarterly, 129600 for monthly)
   
#      # Returns
#      #cycle: cyclical component

#    n = length(y)
#    if n <= 3
#     return zeros(n), y
#    end
   
#    # Setting up the matrix equation
#    A = zeros(n-2, n)
#    for i in 1:(n-2)
#         A[i, i:i+2] .= [1.0, -2.0, 1.0]
#     end
   
#    # Create sparse array and solve
#    D = sparse(A)
#    B = sparse(I(n)) + λ * (D' * D)
#    τ = B \ y
   
#    return y - τ
# end

"""
    hp_filter(y, lambda)

Hodrick-Prescott cyclical component: `y - trend`, where

    trend = argmin_τ  Σ(y_t - τ_t)² + λ Σ(τ_{t+1} - 2τ_t + τ_{t-1})²
          = (I + λ D'D)^{-1} y

with `D` the (n-2)×n **second**-difference operator. `I + λD'D` is pentadiagonal, so
the solve is sparse and O(n).

FIXED September 5, 2026. The previous implementation built the *tridiagonal*
matrix `Tridiagonal(-λ, 1+2λ, -λ)`, which is `I + λD'D` for the **first**-difference
operator — a Whittaker/random-walk smoother, not the HP filter. At λ=1,600 it returned
a "cycle" with ~5x the correct standard deviation and correlation ≈0.26 with the true
HP cycle (verified against `statsmodels.tsa.filters.hp_filter`). Any model-side second
moment computed with the old version is invalid; see context/pending_tasks.md.
"""
function hp_filter(y, lambda::Real)
    n = length(y)
    n < 5 && return y .- mean(y)   # HP is not identified for very short series

    # Second-difference operator D: (n-2) x n, rows [1 -2 1]
    D = spdiagm(n - 2, n, 0 => ones(n - 2), 1 => fill(-2.0, n - 2), 2 => ones(n - 2))

    A = sparse(I, n, n) + lambda * (D' * D)   # pentadiagonal
    trend = A \ collect(float.(y))

    return y .- trend
end

function time_series_object(out::Matrix, fields::Vector{Symbol})
    sim_length = size(out)[1]
    initial = Date(1800, 1, 1)
    final = initial + Month(sim_length-1)
    date_seq = initial:Month(1):final

    df = TS(out, date_seq);
    df.fields = fields
    return df
end

function binvec(x, n::Int,
    rng::AbstractRNG=Random.default_rng())
    n > 0 || throw(ArgumentError("number of bins must be positive"))
    l = length(x)

    # find bin sizes
    d, r = divrem(l, n)
    lens = fill(d, n)
    lens[1:r] .+= 1
    # randomly decide which bins should be larger
    shuffle!(rng, lens)

    # ensure that we have data sorted by x, but ties are ordered randomly
    df = DataFrame(id=axes(x, 1), x=x, r=rand(rng, l))
    sort!(df, [:x, :r])

    # assign bin ids to rows
    binids = reduce(vcat, [fill(i, v) for (i, v) in enumerate(lens)])
    df.binids = binids

    # recover original row order
    sort!(df, :id)
    return df.binids
end

function estimate_markov(X::Vector, nstates::Int64)
    P = zeros(nstates, nstates)
    n = length(X) - 1
    for t = 1:n
        # add one each count is observed
        P[X[t], X[t+1]] = P[X[t], X[t+1]] + 1
    end
    for i = 1:nstates
        # divide by total number of counts
        P[i, :] .= P[i, :]/sum(P[i, :])
    end
    return P
end



"""
Monthly-to-quarterly averages by time means
"""
function monthly_to_quarterly(df)
    quarterly_values = []
    num_cols = ncol(df)
    
    for col in names(df)
        quarterly_col_values = []
        for i in 1:3:size(df, 1)-2
            quarterly_value = mean(df[i:i+2, col])
            push!(quarterly_col_values, float(quarterly_value))
        end
        push!(quarterly_values, quarterly_col_values)
    end
    df_quarterly = hcat([quarterly_values[idx] for idx in 1:num_cols]...)
    df_quarterly = DataFrame(df_quarterly, names(df))
    return df_quarterly
end