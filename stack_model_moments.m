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

% u_obs, v_obs, theta_obs, z_obs, s_obs, delta_obs bf_obs;

n_obs = 7;
n_mom = n_obs*2 + n_obs*(n_obs-1)/2;
n_var = M_.orig_endo_nbr;

oo_.var = get_variance_of_endogenous_variables ...
    (oo_.dr, oo_.dr.inv_order_var(n_var-2*n_obs+1:n_var));

model_moments = nan(n_mom,1);

k = 0;
for i = 1:n_obs
    k = k+1;
    model_moments(k) ...
        = sqrt(oo_.var(i,i));
end
for i = 1:n_obs
    for j = i+1:n_obs
        k = k+1;
        model_moments(k) ...
            = oo_.var(i,j)/(sqrt(oo_.var(i,i))*sqrt(oo_.var(j,j)));
    end
end
for i = 1:n_obs
    k = k+1;
    model_moments(k) ...
        = oo_.var(i,i+n_obs)/(sqrt(oo_.var(i,i))*sqrt(oo_.var(i+n_obs,i+n_obs)));
end