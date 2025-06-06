%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMS - Parameter initialization function for nested DSGE model
%
% This function defines all model parameters and steady-state values for a
% dynamic stochastic general equilibrium (DSGE) model. Parameters are
% organized into four categories:
%   1. Predetermined parameters
%   2. Predetermined steady-state values
%   3. Target steady-state values
%   4. Estimated parameters
%
% Outputs: Returns all parameters needed for model simulation and estimation
%
% Version History:
% Created by Zhesheng Qiu at 15:48 2019-01-14
% Multiple updates through 14:28 2025-09-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ ...
    beta, eta_L, tau, x_v, ...                % predetermined parameters
    f_ss, q_ss, N_ss, w_ss, ls_ss, ...          % Target steady-state values
    sigma, b_ratio, xi_inv, delta, epsi, ...      % Estimated parameters
        rho_z, rho_s, rho_delta, sigma_z, sigma_s, sigma_delta ...
    ] = params


%% Predetermined parameters - Fixed model coefficients
beta       = 0.99673;   % Discount factor (4% annual interest rate)
eta_L      = 0.6;       % Elasticity of matching function wrt unemployment
tau        = 0.031;     % Aggregate separation rate %tau = s+ delta(1-s) \approx s + delta
x_v        = 1.00;      % share of sunk/non-fixed matching costs to total costs: κ = (1-x_v)/x_v*K/q


%% Predetermined steady state - Growth and inflation parameters


%% The steady state to be targeted - Macroeconomic ratios and levels
f_ss        = 0.41;     % Job finding rate (41% monthly)
q_ss        = 0.8;      % Vacancy filling rate (80% monthly)
N_ss        = 1;        % Mass of firms (normalization)
w_ss        = 1;        % Wage rate (normalization)
ls_ss       = 0.66;     % Labor share of income


%% Estimated parameters - estimated via Bayesian methods

% Labor market parameters
sigma      = 1.0;       % Inverse of intertemporal elast. of sub. (log utility)
b_ratio     = 0.71;     % Outside option of worker
xi_inv      = 1/0.265;  % inverse elasticity of entry to vacancy value

% Entry and monopolistic competition
delta       = 0.00874;  % Product destruction rate (set to 10% annual)
epsi        = 1.3/0.3;  % Prod. elasticity of substitution (set to 30% annual markups)

% Shock process parameters
rho_z        = 0.963;   % Persistence of tech shocks
rho_s        = 0.78;    % Persistence of separation shocks
rho_delta    = 0.935;   % Persistence of product destruction shocks
sigma_z      = 0.0055;  % Std. dev of tech shocks
sigma_s      = 0.051;   % std. dev of separation shocks
sigma_delta  = 0.0197;  % Std. dev of prod. dest. shocks

% Note, given these we obtain dependent parameters A, F, zbar, phi, f_e
% A: level parameter of matching function 
% F: Mass of potential recruiters 
% zbar: level of technology
% phi: Bargaining power of worker
% f_e: sunk entry cost of firm