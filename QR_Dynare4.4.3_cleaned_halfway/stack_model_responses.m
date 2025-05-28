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

%% Monetary Policy Shock Responses
model_resp_mon = zeros(16,15);  % Initialize 16q × 15var matrix

% Cumulative investment technology effect (for relative price calculation)
muilev_epsR = cumsum(mux_epsR);

% Fill response matrix (period 1 = impact response)
model_resp_mon(2:end,1)  = Y_epsR(1:end-1);         % Real GDP (%)
model_resp_mon(2:end,2)  = C_epsR(1:end-1);         % Real consumption (%)
model_resp_mon(2:end,3)  = I_epsR(1:end-1);         % Real investment (%)
model_resp_mon(2:end,4)  = L_epsR(1:end-1);         % Hours worked (%)
model_resp_mon(2:end,5)  = u_epsR(1:end-1);         % Capacity utilization (APR)
model_resp_mon(2:end,6)  = -muilev_epsR(1:end-1);   % Relative price of investment (%)
model_resp_mon(2:end,7)  = w_epsR(1:end-1);         % Real wage (%)
model_resp_mon(2:end,8)  = 4*Pi_epsR(1:end-1);      % Inflation (GDP deflator, APR)
model_resp_mon(:,9)      = 4*R_epsR;                % Fed fund rate (APR)
model_resp_mon(2:end,10) = lp_epsR(1:end-1);        % Labor productivity (%)
model_resp_mon(2:end,11) = ls_epsR(1:end-1);        % Labor share (%)
model_resp_mon(2:end,12) = Ical_epsR(1:end-1);      % Occupancy (%)
model_resp_mon(2:end,13) = tau_epsR(1:end-1);       % Markup (p.p.)
model_resp_mon(2:end,14) = model_resp_mon(2:end,7) + L_epsR(1:end-1) + Ra_epsR(1:end-1);
                                                    % Working capital (%)
model_resp_mon(2:end,15) = uk_epsR(1:end-1);        % Capital utilization (APR)
model_resp_mon = 100*model_resp_mon;                % Convert to percentage points


%% Neutral Technology Shock Responses
model_resp_ntech = zeros(16,15);  % Initialize matrix

% Cumulative growth effects
muylev_epsn = cumsum(muy_epsn);   % Output growth
muilev_epsn = cumsum(mux_epsn);   % Investment tech growth

% Fill response matrix (all periods included)
model_resp_ntech(:,1)  = Y_epsn + muylev_epsn;      % Real GDP (%)
model_resp_ntech(:,2)  = C_epsn + muylev_epsn;      % Real consumption (%)
model_resp_ntech(:,3)  = I_epsn + muylev_epsn + muilev_epsn;
                                                    % Real investment (%)
model_resp_ntech(:,4)  = L_epsn;                    % Hours worked (%)
model_resp_ntech(:,5)  = u_epsn;                    % Capacity utilization (APR)
model_resp_ntech(:,6)  = -muilev_epsn;              % Relative price of investment (%)
model_resp_ntech(:,7)  = w_epsn + muylev_epsn;      % Real wage (%)
model_resp_ntech(:,8)  = 4*Pi_epsn;                 % Inflation (GDP deflator, APR)
model_resp_ntech(:,9)  = 4*Ra_epsn;                 % Fed fund rate (APR)
model_resp_ntech(:,10) = lp_epsn + muylev_epsn;     % Labor productivity (%)
model_resp_ntech(:,11) = ls_epsn;                   % Labor share (%)
model_resp_ntech(:,12) = Ical_epsn;                 % Occupancy (%)
model_resp_ntech(:,13) = tau_epsn;                  % Markup (p.p.)
model_resp_ntech(:,14) = model_resp_ntech(:,7) + L_epsn + Ra_epsn;
                                                    % Working capital (%)
model_resp_ntech(:,15) = uk_epsn;                   % Capital utilization (APR)
model_resp_ntech = 100*model_resp_ntech;            % Convert to percentage points


%% Investment Technology Shock Responses 
model_resp_itech = zeros(16,15);    % Initialize matrix

% Cumulative growth effects
muylev_epsx = cumsum(muy_epsx);     % Output growth
muilev_epsx = cumsum(mux_epsx);     % Investment tech growth

% Fill response matrix
xs_SS = 1-1/(1+17.5/62.5);          % Investment-to-output ratio (17.5/62.5 = C/K ratio)
model_resp_itech(:,1)  = Y_epsx + muylev_epsx + xs_SS*muilev_epsx;
                                                    % Real GDP (%)
model_resp_itech(:,2)  = C_epsx + muylev_epsx + xs_SS*muilev_epsx;
                                                    % Real consumption (%)
model_resp_itech(:,3)  = I_epsx + muylev_epsx + xs_SS*muilev_epsx;
                                                    % Real investment (%)
model_resp_itech(:,4)  = L_epsx;                    % Hours worked (%)
model_resp_itech(:,5)  = u_epsx;                    % Capacity utilization (APR)
model_resp_itech(:,6)  = -muilev_epsx + xs_SS*muilev_epsx;
                                                    % Relative price of investment (%)
model_resp_itech(:,7)  = w_epsx + muylev_epsx + xs_SS*muilev_epsx;
                                                    % Real wage (%)
model_resp_itech(:,8)  = 4*Pi_epsx - mux_epsx;      % Inflation (GDP deflator, APR)
model_resp_itech(:,9)  = 4*Ra_epsx;                 % Fed fund rate (APR)
model_resp_itech(:,10) = lp_epsx + muylev_epsx + xs_SS*muilev_epsx;
                                                    % Labor productivity (%)
model_resp_itech(:,11) = ls_epsx;                   % Labor share (%)
model_resp_itech(:,12) = Ical_epsx;                 % Occupancy (%)
model_resp_itech(:,13) = tau_epsx;                  % Markup (p.p.)
model_resp_itech(:,14) = model_resp_itech(:,7) + L_epsx + Ra_epsx;
                                                    % Working capital (%)
model_resp_itech(:,15) = uk_epsx;                   % Capital utilization (APR)
model_resp_itech = 100*model_resp_itech;            % Convert to percentage points