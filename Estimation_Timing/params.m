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
    eta, F, epsi, zeta, xy, ...                 % Predetermined parameters
    r_ss, delta_ss, tau_ss, ...
    f_ss, q_ss, N_ss, w_ss, ...                 % Target steady-state values    
    sigma, b_ratio, xi_inv, x_v, ...            % Estimated parameters
        rho_z, rho_s, rho_delta, sigma_z, sigma_s, sigma_delta ...
    ] ...
    = params


%% Predetermined parameters
eta         = 0.6;      % Share of unemployment in matching function
F           = 1.0;      % Number of recruiters
epsi        = 4.3;      % Elasticity of substitution between varieties
zeta        = 1/3.3;    % Love-for-variety
xy          = 0.015;    % Share of hiring cost in GDP

%% Target steady state
r_ss        = 1.04^(1/12) - 1;      % real interest rate (set to 4% annual)
delta_ss    = 1-(1-0.0754)^(1/12);  % Product destruction rate (set to 7.54% annual)
tau_ss      = 0.031;    % Aggregate separation rate %tau = s+ delta(1-s) \approx s + delta
f_ss        = 0.41/(1-delta_ss);    % Job finding rate (41% monthly)
q_ss        = 0.8/(1-delta_ss);     % Vacancy filling rate (80% monthly)
N_ss        = 1;        % Mass of firms (normalization)
w_ss        = 1;        % Wage rate (normalization)


%% Estimated parameters - estimated via Bayesian methods

% Labor market parameters
sigma       = 2.0;      % Inverse of intertemporal elast. of sub. (log utility)
b_ratio     = 0.7;      % Outside option of worker
xi_inv      = 2.0;      % inverse elasticity of entry to vacancy value
x_v         = 0.5;      % share of hiring costs from vacancy creation

% Shock process parameters
rho_z        = 0.963;   % Persistence of tech shocks
rho_s        = 0.78;    % Persistence of separation shocks
rho_delta    = 0.935;   % Persistence of product destruction shocks
sigma_z      = 0.0055;  % Std. dev of tech shocks
sigma_s      = 0.051;   % std. dev of separation shocks
sigma_delta  = 0.0197;  % Std. dev of prod. dest. shocks