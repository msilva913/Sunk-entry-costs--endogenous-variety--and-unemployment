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
    eta, F, epsi, zeta, xy, ...                 % Predetermined parameters
    r_ss, delta_ss, tau_ss, ...
    f_ss, q_ss, N_ss, w_ss, ...                 % Target steady-state values
    sigma, b_ratio, xi_inv, x_v, ...            % Estimated parameters
        rho_z, rho_s, rho_delta, sigma_z, sigma_s, sigma_delta ...
    ] ...
    = params;
else
    [ ...
    eta, F, epsi, zeta, xy, ...                 % Predetermined parameters
    r_ss, delta_ss, tau_ss, ...
    f_ss, q_ss, N_ss, w_ss, ...                 % Target steady-state values
    ] ...
    = params;
    
    sigma       = M_.params(1);
    b_ratio     = M_.params(2);
    xi_inv      = M_.params(3);
    x_v         = M_.params(4);
    rho_z       = M_.params(5);    
    rho_s       = M_.params(6);    
    rho_delta   = M_.params(7);
    sigma_z     = M_.params(8);
    sigma_s     = M_.params(9);
    sigma_delta = M_.params(10);
end


%% Solve steady state

% Abbreviation
mu          = epsi/(epsi-1);
hc          = (r_ss+tau_ss)/(1-delta_ss);
kq          = (r_ss+delta_ss)/(1+r_ss);
lsg         = (delta_ss+(epsi-1)*(r_ss+delta_ss)) ...
            /(delta_ss+epsi*(r_ss+delta_ss));

% Labor market variables
theta_ss    = f_ss/q_ss;
u_ss        = tau_ss/(tau_ss+f_ss*(1-delta_ss));
L_ss        = 1 - u_ss;
v_ss        = theta_ss*u_ss;
e_ss        = delta_ss*(v_ss+1-u_ss);
s_ss        = (tau_ss-delta_ss)/(1-delta_ss);

ws          = (1+hc/(x_v*q_ss))/(e_ss/(1+xi_inv)/kq+(1-x_v)/x_v*v_ss) ...
            * xy/(1+xy)* L_ss/lsg;

% Pin down the levels
wint_ss     = w_ss/(1-ws);
Ygross_ss   = L_ss*wint_ss/lsg;
X_ss        = xy/(1+xy)*Ygross_ss;
Y_ss        = Ygross_ss - X_ss;
Q_ss        = X_ss/(e_ss/(1+xi_inv)+(1-x_v)/x_v*kq*v_ss);
K_ss        = kq*Q_ss;

% business formation
p_ss        = N_ss^zeta;
z_ss        = (mu/p_ss)*wint_ss;
f_e         = (mu-1)*z_ss*(L_ss/N_ss)*(1-delta_ss)/(r_ss+mu*delta_ss);
nu_ss       = wint_ss*f_e/z_ss;
d_ss        = (r_ss+delta_ss)/(1-delta_ss)*nu_ss;

% National accounting
Ne_ss       = delta_ss/(1-delta_ss)*N_ss;
Le_ss       = Ne_ss*f_e/z_ss;
Lc_ss       = L_ss - Le_ss;
% Lc_ss       = (r_ss+delta_ss)*L_ss/(r_ss+mu*delta_ss);
% Le_ss       = (mu-1)*delta_ss*L_ss/(r_ss+mu*delta_ss);
Yc_ss       = p_ss*z_ss*Lc_ss;
C_ss        = Yc_ss - X_ss;
lam_ss      = C_ss^(-sigma);
ls_ss       = w_ss*L_ss/Y_ss;
lp_ss       = Y_ss/(p_ss*L_ss);

% Other implied value of parameters
beta        = 1/(1+r_ss);
kappa       = (1-x_v)/(x_v*q_ss)*K_ss;
A           = f_ss*theta_ss^(eta-1);
x_m         = Q_ss*(e_ss/F)^(-xi_inv);
phi         = (1-b_ratio) ...
            *w_ss/(wint_ss-K_ss+f_ss*(kappa+K_ss/q_ss)-b_ratio*w_ss);

% Observables
u_obs_ss    = log(u_ss);
v_obs_ss    = log(v_ss);
tau_obs_ss  = log(3*tau_ss);
f_obs_ss    = log(3*(1-delta_ss)*f_ss);
delta_obs_ss= log(3*delta_ss);
bf_obs_ss   = log(3*Ne_ss);
lp_obs_ss   = log(lp_ss);


%% Put steady states into ys vector

ys          = nan(M_.orig_endo_nbr,1);

ys(1)       = theta_ss;
ys(2)       = f_ss;
ys(3)       = q_ss;

ys(4)       = p_ss;

ys(5)       = lam_ss;
ys(6)       = nu_ss;
ys(7)       = d_ss;
ys(8)       = Ne_ss;

ys(9)       = K_ss;
ys(10)      = Q_ss;
ys(11)      = e_ss;
ys(12)      = wint_ss;
ys(13)      = w_ss;

ys(14)      = N_ss;
ys(15)      = u_ss;
ys(16)      = v_ss;
ys(17)      = tau_ss;

ys(18)      = L_ss;
ys(19)      = Le_ss;
ys(20)      = Lc_ss;
ys(21)      = X_ss;
ys(22)      = C_ss;
ys(23)      = Yc_ss;
ys(24)      = Ygross_ss;
ys(25)      = Y_ss;
ys(26)      = ls_ss;
ys(27)      = lp_ss;

ys(28)      = z_ss;
ys(29)      = s_ss;
ys(30)      = delta_ss;

ys(31)      = u_obs_ss;
ys(32)      = v_obs_ss;
ys(33)      = tau_obs_ss;
ys(34)      = f_obs_ss;
ys(35)      = delta_obs_ss;
ys(36)      = bf_obs_ss;
ys(37)      = lp_obs_ss;

ys(38:44)   = ys(31:37);


%% Save parameters

if isnan(M_.params)
    
    M_.params       = nan(M_.param_nbr,1);
    
    M_.params(1)    = sigma;
    M_.params(2)    = b_ratio;
    M_.params(3)    = xi_inv;
    M_.params(4)    = x_v;
    M_.params(5)    = rho_z;
    M_.params(6)    = rho_s;
    M_.params(7)    = rho_delta;
    M_.params(8)    = sigma_z;
    M_.params(9)    = sigma_s;
    M_.params(10)   = sigma_delta;
        
end

M_.params(11)   = eta;
M_.params(12)   = F;
M_.params(13)   = epsi;
M_.params(14)   = zeta;
M_.params(15)   = xy;

M_.params(16)   = beta;
M_.params(17)   = kappa;
M_.params(18)   = f_e;
M_.params(19)   = A;
M_.params(20)   = x_m;
M_.params(21)   = phi;

% disp('------ Parameter Values ------');
% disp(['beta = ', num2str(beta)]);
% disp(['delta = ', num2str(delta)]);  % Check delta is in (0,1)
% disp(['tau = ', num2str(tau)]);      % Must be > delta
% disp(['epsi = ', num2str(epsi)]);    % Must be > 1
% disp(['xi_inv = ', num2str(xi_inv)]); % Must allow real roots

% disp('------ Parameter Values ------');
% disp(M_.params);
% disp('------ Steady State Values ------');
% disp(ys);

check = 0;