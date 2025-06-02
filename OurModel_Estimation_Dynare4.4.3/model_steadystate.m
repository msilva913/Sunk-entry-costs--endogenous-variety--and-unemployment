%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function sets parameters and steady state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ys,check] = model_steadystate(ys,exo)
 
global M_


%% Read parameters

% if isempty(M_) == 1
%     M_.params = NaN;
% end

if isnan(M_.params)
    [ ...
    beta, eta_L, tau, ...                % predetermined parameters
    f_ss, q_ss, N_ss, w_ss, ls_ss, ...          % Target steady-state values
    sigma, b_ratio, x_v, xi_inv, delta, epsi, ...      % Estimated parameters
        rho_z, rho_delta, rho_s, sigma_z, sigma_delta, sigma_s ...
    ] ...
    = params;
else
    [ ...
    beta, eta_L, tau, ...                % predetermined parameters
    f_ss, q_ss, N_ss, w_ss, ls_ss, ...          % Target steady-state values
    ] ...
    = params;
    sigma       = M_.params(1);
    b_ratio     = M_.params(2);    
    x_v         = M_.params(3);
    xi_inv      = M_.params(4);    
    delta       = M_.params(5);
    epsi        = M_.params(6);
    rho_z       = M_.params(7);    
    rho_delta   = M_.params(8);    
    rho_s       = M_.params(9);
    sigma_z     = M_.params(10);
    sigma_delta = M_.params(11);
    sigma_s     = M_.params(12);
end


%% Solve steady state

mu_ss       = epsi/(epsi-1);
s_agg_ss    = tau;
s_ss        = (tau-delta)/(1-delta);

p_ss        = N_ss^(1/(epsi-1));
N_e_ss      = delta/(1-delta)*N_ss;
    
rho_ss      = (1-beta)/beta;
fbar        = f_ss;
qbar        = q_ss;
f_ss        = fbar/(1-delta);
q_ss        = qbar/(1-delta);
    
theta_ss    = f_ss/q_ss;
u_ss        = tau/(tau+(1-delta)*f_ss);
L_ss        = 1 - u_ss;
v_ss        = theta_ss*u_ss;
e_ss        = delta*(v_ss+1-u_ss);
    
rs_ss       = (delta+(rho_ss+delta)*(epsi-1))/(delta+(rho_ss+delta)*epsi);
w_wR_ss     = ls_ss/rs_ss;
w_R_ss      = w_ss/(w_wR_ss);
z_ss        = (mu_ss/p_ss)*w_R_ss;
    
sr_ss       = (rho_ss+tau)/(1-delta)*(1/(q_ss*x_v));
K_ss        = (w_R_ss-w_ss)/(1+sr_ss);
kappa_ss    = (1-x_v)/x_v*K_ss/q_ss;
    
f_e_ss      = (mu_ss-1)*z_ss*(L_ss/N_ss)*(1-delta)/(delta*mu_ss+rho_ss);
    
nu_f_ss     = p_ss*f_e_ss/mu_ss;
d_f_ss      = (rho_ss+delta)/(1-delta)*nu_f_ss;
    
L_e_ss      = N_e_ss*f_e_ss/z_ss;
L_c_ss      = L_ss - L_e_ss;
Y_c_ss      = p_ss*z_ss*L_c_ss;
Q_ss        = K_ss*(1+rho_ss)/(rho_ss+delta);
    
F_ss        = e_ss/(Q_ss^(1/xi_inv));
    
Omega_ss    = F_ss/(1+xi_inv)*(e_ss/F_ss)^(1+xi_inv);
C_ss        = Y_c_ss - Omega_ss - kappa_ss*v_ss*q_ss;
lam_ss      = C_ss^(-sigma);
Y_ss        = Y_c_ss + nu_f_ss*N_e_ss;

A           = f_ss/(theta_ss^(1-eta_L));
f_e         = (mu_ss-1)*z_ss*(L_ss/N_ss)*(1-delta)/(delta*mu_ss+rho_ss);
kappa       = (1-x_v)/x_v*K_ss/q_ss;
F           = F_ss;
phi         = (w_ss-b_ratio*w_ss)/(w_R_ss-K_ss+theta_ss/(1-delta)*(K_ss+q_ss*kappa)-b_ratio*w_ss);

% A           = f_ss/(theta_ss^(1-eta_L));
% f_e         = nu_f_ss*mu_ss/p_ss;
% kappa       = (1-x_v)/x_v*K_ss/q_ss;
% F           = e_ss/(Q_ss^(1/xi_inv));
% phi         = (w_ss-b_ratio*w_ss)/(w_R_ss-K_ss+theta_ss/(1-delta)*(K_ss+q_ss*kappa)-b_ratio*w_ss);

log_z_ss    = log(z_ss);
log_delta_ss= log(delta);
log_s_ss    = log(s_ss);

u_obs_ss    = log(u_ss);
v_obs_ss    = log(v_ss);
theta_obs_ss= log(v_ss/u_ss);
lp_obs_ss   = z_ss;
s_obs_ss    = log(3*s_ss);
bf_obs_ss   = log(3*N_e_ss);


%% Put steady states into ys vector

ys          = nan(40,1);

ys(1)       = theta_ss;
ys(2)       = f_ss;
ys(3)       = q_ss;

ys(4)       = p_ss;

ys(5)       = lam_ss;
ys(6)       = nu_f_ss;
ys(7)       = d_f_ss;
ys(8)       = N_e_ss;

ys(9)       = K_ss;
ys(10)      = Q_ss;
ys(11)      = e_ss;
ys(12)      = w_R_ss;
ys(13)      = w_ss;

ys(14)      = N_ss;
ys(15)      = u_ss;
ys(16)      = v_ss;

ys(17)      = L_ss;
ys(18)      = L_e_ss;
ys(19)      = L_c_ss;
ys(20)      = Omega_ss;
ys(21)      = C_ss;
ys(22)      = Y_c_ss;
ys(23)      = Y_ss;
ys(24)      = ls_ss;
ys(25)      = s_agg_ss;

ys(26)      = log_z_ss;
ys(27)      = log_delta_ss;
ys(28)      = log_s_ss;

ys(29)      = u_obs_ss;
ys(30)      = v_obs_ss;
ys(31)      = theta_obs_ss;
ys(32)      = lp_obs_ss;
ys(33)      = s_obs_ss;
ys(34)      = bf_obs_ss;

ys(35)      = u_obs_ss;
ys(36)      = v_obs_ss;
ys(37)      = theta_obs_ss;
ys(38)      = lp_obs_ss;
ys(39)      = s_obs_ss;
ys(40)      = bf_obs_ss;


%% Save parameters

if isnan(M_.params)
    
    M_.params       = nan(20,1);
    M_.params(1)    = sigma;   
    M_.params(2)    = b_ratio;
    M_.params(3)    = x_v;
    M_.params(4)    = xi_inv;
    M_.params(5)    = delta;
    M_.params(6)    = epsi;
    M_.params(7)    = rho_z;
    M_.params(8)    = rho_delta;
    M_.params(9)    = rho_s;
    M_.params(10)    = sigma_z;
    M_.params(11)   = sigma_delta;
    M_.params(12)   = sigma_s;
        
end

M_.params(13)   = beta;
M_.params(14)   = eta_L;
%M_.params(14)   = sigma;
M_.params(15)   = tau;

M_.params(16)   = A;
M_.params(17)   = f_e;
M_.params(18)   = kappa;
M_.params(19)   = F;
M_.params(20)   = phi;

% disp('------ Parameter Values ------');
% disp(['beta = ', num2str(beta)]);
 disp(['sigma = ', num2str(sigma)]);
% disp(['delta = ', num2str(delta)]);  % Check delta is in (0,1)
% disp(['tau = ', num2str(tau)]);      % Must be > delta
% disp(['epsi = ', num2str(epsi)]);    % Must be > 1
% disp(['xi_inv = ', num2str(xi_inv)]); % Must allow real roots

% disp('------ Parameter Values ------');
% disp(M_.params);
% disp('------ Steady State Values ------');
% disp(ys);

check = 0;