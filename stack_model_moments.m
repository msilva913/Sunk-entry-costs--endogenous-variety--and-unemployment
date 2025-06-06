%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS_MODEL_IRFS - Prepares model impulse responses for comparison with VAR
%
% This code processes model-generated impulse response functions (IRFs) for:
%   1. Monetary policy shocks
%   2. Neutral technology shocks
%   3. Investment technology shocks
%
% For each shock type, it:
%   1. Creates a 16-period × 15-variable response matrix
%   2. Adjusts variables for cumulative growth effects where needed
%   3. Annualizes rate variables (×4) 
%   4. Converts to percentage points (×100)
%
% Structure:
%   Each shock section follows the same variable ordering for consistency
%   Rows represent time periods (0-15 quarters)
%   Columns represent economic variables (see details below)
%
% Variable order (columns):
%   1-Y, 2-C, 3-I, 4-L, 5-u, 6-inv_price, 7-w, 8-Pi, 9-R/Ra, 
%   10-lp, 11-ls, 12-Ical, 13-tau, 14-composite, 15-uk
%
% Note: The 'composite' variable (col 14) combines wage, labor and rate effects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% u_obs, v_obs, theta_obs, z_obs, s_obs, bf_obs;

model_moments = nan(27,1);

model_moments( 1) = sqrt(oo_.var(1,1));
model_moments( 2) = sqrt(oo_.var(2,2));
model_moments( 3) = sqrt(oo_.var(3,3));
model_moments( 4) = sqrt(oo_.var(4,4));
model_moments( 5) = sqrt(oo_.var(5,5));
model_moments( 6) = sqrt(oo_.var(6,6));
model_moments( 7) = oo_.var(1,2)/(sqrt(oo_.var(1,1))*sqrt(oo_.var(2,2)));
model_moments( 8) = oo_.var(1,3)/(sqrt(oo_.var(1,1))*sqrt(oo_.var(3,3)));
model_moments( 9) = oo_.var(1,4)/(sqrt(oo_.var(1,1))*sqrt(oo_.var(4,4)));
model_moments(10) = oo_.var(1,5)/(sqrt(oo_.var(1,1))*sqrt(oo_.var(5,5)));
model_moments(11) = oo_.var(1,6)/(sqrt(oo_.var(1,1))*sqrt(oo_.var(6,6)));
model_moments(12) = oo_.var(2,3)/(sqrt(oo_.var(2,2))*sqrt(oo_.var(3,3)));
model_moments(13) = oo_.var(2,4)/(sqrt(oo_.var(2,2))*sqrt(oo_.var(4,4)));
model_moments(14) = oo_.var(2,5)/(sqrt(oo_.var(2,2))*sqrt(oo_.var(5,5)));
model_moments(15) = oo_.var(2,6)/(sqrt(oo_.var(2,2))*sqrt(oo_.var(6,6)));
model_moments(16) = oo_.var(3,4)/(sqrt(oo_.var(3,3))*sqrt(oo_.var(4,4)));
model_moments(17) = oo_.var(3,5)/(sqrt(oo_.var(3,3))*sqrt(oo_.var(5,5)));
model_moments(18) = oo_.var(3,6)/(sqrt(oo_.var(3,3))*sqrt(oo_.var(6,6)));
model_moments(19) = oo_.var(4,5)/(sqrt(oo_.var(4,4))*sqrt(oo_.var(5,5)));
model_moments(20) = oo_.var(4,6)/(sqrt(oo_.var(4,4))*sqrt(oo_.var(6,6)));
model_moments(21) = oo_.var(5,6)/(sqrt(oo_.var(5,5))*sqrt(oo_.var(6,6)));
model_moments(22) = oo_.autocorr{1,3}(1,1);
model_moments(23) = oo_.autocorr{1,3}(2,2);
model_moments(24) = oo_.autocorr{1,3}(3,3);
model_moments(25) = oo_.autocorr{1,3}(4,4);
model_moments(26) = oo_.autocorr{1,3}(5,5);
model_moments(27) = oo_.autocorr{1,3}(6,6);