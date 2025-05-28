%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP_IRF_MATCHING - Prepares impulse response data for model estimation
%
% Purpose:
%   Loads and processes empirical impulse response functions (IRFs) and their
%   standard errors for use in DSGE model estimation via impulse response matching.
%
% Inputs:
%   Requires MAT-file 'SVAR_IRFsAndSEs' containing:
%     - IRFout: Empirical impulse responses (horizon × variables × shocks)
%     - IRFse: Standard errors for impulse responses
%     - namefig: Variable names for the IRFs
%
% Outputs:
%   Sets global variables containing processed IRF data for estimation:
%     - psihat: Vector of empirical IRFs to match
%     - Vhat: Diagonal covariance matrix of IRF standard errors
%     - inv_Vhat: Inverse of Vhat (for weighting in estimation)
%     - logdetVhat: Log determinant of Vhat (for likelihood calculation)
%
% Version History:
%   Created by Zhesheng Qiu at 00:10 2019-01-29
%   Multiple updates through 16:37 2025-05-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load empirical impulse response data
load SVAR_IRFsAndSEs;  % Contains IRFout, IRFse, and namefig variables


%% Specify observables and model variables
global horizon logdetVhat inv_Vhat psihat data_varnames mod_var_list mod_shock_list

% Set estimation horizon based on loaded IRF data
horizon = size(IRFout, 1);  % Number of time periods in the IRFs

% Create cell array of observable variable names from namefig matrix
data_varnames = arrayfun(@(i) {namefig(i,:)}, 1:size(namefig,1))';

% Define the complete lists of model variables and shocks to construct IRFs
mod_var_list = {'Y','C','I','L','u','w','Pi','R','Ra','lp','ls','Ical','tau','muy','mux'};
mod_shock_list = {'epsR','epsn','epsx'};  % Monetary, neutral tech, investment tech shocks

% Create selection vector for which observables to use in estimation
% 1 = include, 0 = exclude
select_vec = [
    1  % Real GDP (%)
    1  % Real consumption (%)
    1  % Real investment (%)
    1  % Hours worked (%)
    1  % Capacity utilization (APR)
    1  % Relative price of investment (%)
    1  % Real wage (%)
    1  % Inflation (GDP deflator, APR)
    1  % Fed fund rate (APR)
    1  % Labor productivity (%)
    1  % Labor share (%)
    ]';% Include all 11 variables by default

% Get indices of selected variables (global for use in other functions)
global select_idx
if isempty(select_idx)
    select_idx = find(select_vec == 1);
end


%% Process impulse response data
% Extract IRFs and their standard errors for each shock type:
% 1. Monetary policy shock (epsR)
% 2. Neutral technology shock (epsn) 
% 3. Investment technology shock (epsx)

% Note: Negative sign on IRFR for monetary policy shock by convention
IRFR = -IRFout(:,select_idx,3);  % Monetary policy shock responses
IRFn = IRFout(:,select_idx,2);   % Neutral tech shock responses
IRFx = IRFout(:,select_idx,1);   % Investment tech shock responses

% Corresponding standard errors
IRFRse = IRFse(:,select_idx,3);
IRFnse = IRFse(:,select_idx,2);
IRFxse = IRFse(:,select_idx,1);


%% Prepare data structures for estimation

% Stack all IRFs into single column vector
psihat = [IRFR, IRFn, IRFx];  % Combine all shock responses
psihat = psihat(:);           % Vectorize

% Stack all squared standard errors into single column vector
Vhat = [IRFRse, IRFnse, IRFxse].^2;  % Variance of each IRF point
Vhat = Vhat(:);                      % Vectorize variances

% Remove zero elements in Vhat and psihat (due to information restictions)
nonzero_idx = find(psihat ~= 0);
psihat = psihat(nonzero_idx);
Vhat = Vhat(nonzero_idx);

% Create diagonal covariance matrix
Vhat = diag(Vhat);

% Pre-compute matrix operations for estimation efficiency
inv_Vhat = inv(Vhat);          % Inverse covariance matrix for weighting
logdetVhat = logdet(Vhat);     % Log determinant for likelihood calculation