%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script initializes the estimation by setting up model parameters and
% saving them in a format suitable for Metropolis-Hastings algorithm
% 
% Inputs: None (loads parameters from params.m)
% Outputs: Saves model_mh_mode.mat containing:
%          - xparam1: Parameter vector
%          - hh: Initial scaling matrix for MCMC
%          - fval: Initial likelihood value
%          - parameter_names: Cell array of parameter names
%
% Version History:
% Created by Zhesheng Qiu at 16:31 2019-01-29
% Multiple updates through 12:07 2025-09-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Names of estimated parameters
parameter_names = [
    % Structural prameters 
    {'delta'}     % Product destruction rate
    {'b'}         % Outside option of unemployed worker
    {'epsi'}        % elasticity of substitution 
    {'x_v'}         % share of sunk/non-fixed matching costs to total costs: κ = (1-x_v)/x_v*K/q
    {'xi_inv'}     % inverse elasticity of entry to vacancy value

    % Shock processes
    {'rho_z'}     % Persistence of technology shock
    {'rho_delta'}  % Persistence of product destruction shock
    {'rho_s'}       % Persistence of idiosyncratic separation shock
    {'sigma_z'}    % std dev of technology shock
    {'sigma_delta'}  % std dev of destruction shock
    {'sigma_s'}     % std dev of idiosyncratic separation shock
    ];

%% Comments
% Bargaining power phi will be obtained indirectly from labor share of
% income=0.66
% Given tau and delta, we obtain s
% Normalized wage = 1 pins down ss tech level, also means that b represents
% ratio of unemployment benefits to wage


% Load all model parameters from the params.m file
% The params.m file should contain definitions for all these variables
[ ...
    ~, ~, ~, ~, ~, ...
    ~, ~, ~, ~, ...
    ~, ~, ~, ~, ~, ~, ~, ...
    delta, b, epsi, x_v, xi_inv, ...  % Parameters being estimated
    rho_z, rho_delta, rho_s, sigma_z, sigma_delta, sigma_s ...
    ] ...
    = params;

% Create vector of parameters to be estimated
% Note: The transpose operation (') is crucial for correct dimensionality
xparam1 = [delta, b, epsi, x_v, xi_inv, ...
         rho_z, rho_delta, rho_s, sigma_z, sigma_delta, sigma_s]';

% Initialize optimization/MCMC variables
fval = 500;                     % Initial value for the likelihood
hh = 1e4*eye(length(xparam1));  % Initial scaling matrix for MCMC

% Save all relevant variables to a mat file for Metropolis-Hastings algorithm
% This file will be loaded by the estimation routine
save model_mh_mode.mat xparam1 hh fval parameter_names
%%
% 
% 