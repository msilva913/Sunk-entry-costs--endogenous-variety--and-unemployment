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

load moments_bf_empirical_35.mat;

global logdetVhat inv_Vhat psihat

% psihat = Values;
psihat = Values([1:5,7:11,13:16,18:20,22:23,25,27,29:33,35]);
Vhat = 0.01*diag(abs(psihat));

% Pre-compute matrix operations for estimation efficiency
inv_Vhat = inv(Vhat);          % Inverse covariance matrix for weighting
logdetVhat = logdet(Vhat);     % Log determinant for likelihood calculation