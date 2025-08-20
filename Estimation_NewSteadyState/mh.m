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
    {'sigma'}       % risk aversion paramter
    {'b_ratio'}     % outside option of unemployed worker
    {'xi_inv'}      % inverse elasticity of entry to vacancy value
    {'x_v'}         % share of hiring costs from vacancy creation
    % Shock processes
    {'rho_z'}       % persistence of technology shock
    {'rho_s'}       % persistence of idiosyncratic separation shock
    {'rho_delta'}   % persistence of product destruction shock
    {'sigma_z'}     % std dev of technology shock
    {'sigma_s'}     % std dev of idiosyncratic separation shock
    {'sigma_delta'} % std dev of destruction shock
    ];

% Load all model parameters from the params.m file
% The params.m file should contain definitions for all these variables
[ ...
    ~, ~, ~, ~, ~, ...
    ~, ~, ~, ~, ~, ~, ~, ...
    sigma, b_ratio, xi_inv, x_v, ...  % Parameters being estimated
    rho_z, rho_s, rho_delta, sigma_z, sigma_s, sigma_delta ...
    ] ...
    = params;

% Create vector of parameters to be estimated
% Note: The transpose operation (') is crucial for correct dimensionality
xparam1 = [sigma, b_ratio, xi_inv, x_v, ...
         rho_z, rho_s, rho_delta, sigma_z, sigma_s, sigma_delta]';

% Initialize optimization/MCMC variables
fval = 500;                     % Initial value for the likelihood
hh = 1e5*eye(length(xparam1));  % Initial scaling matrix for MCMC

% Save all relevant variables to a mat file for Metropolis-Hastings algorithm
% This file will be loaded by the estimation routine
save model_mh_mode.mat xparam1 hh fval parameter_names