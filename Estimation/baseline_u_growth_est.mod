

// define variables

var //x           ${x}$ (long_name='technology')
    Q         ${Q}$ (long_name='value of vacancy')
    K           ${K}$ (long_name='expected discounted difference in vacancy value')
    w        ${w}$ (long_name = 'wage')
    e           ${e}$ (long_name = 'new vacancies')
    v         ${v}$ (long_name = 'vacancy stock')
    u         ${u}$ (long_name = 'unemployment')
    theta     ${theta}$ (long_name='market tightness')
    f         ${f}$ (long_name = 'job finding rate')
    q         ${q}$ (long_name = 'vacancy filling rate')

    N         ${N}$ (long_name='number of firms')
    N_e       ${N_e}$ (long_name='number of new businesses')
    p         ${p}$ (long_name='relative price')
    w_R       ${w_R}$ (long_name='recruiter payment')
    nu_f      ${nu_f}$ (long_name='retailer firm value')
    d_f       ${d_f}$ (long_name='dividend')

    Y_c       ${Y_c}$ (long_name='retail output')
    C         ${C}$ (long_name= 'consumption')
    Omega     ${\Omega}$ (long_name= 'sunk entry costs of investment')
    lam       ${\lambda}$ (long_name='Marginal utility of consumption')
    Y         
    Y_R         ${Y}$ (long_name='data-consistent aggregate output')

    L         ${L}$ (long_name='total labor')
    L_e       ${L_e}$ (long_name='Labor in entry')
    L_c       ${L_c}$ (long_name='Labor in consumption')
    labor_share_obs (long_name='labor share')
    cons_share_obs_m  (long_name='consumptions share of output: monthly')
    cons_share_obs  (long_name='consumptions share of output: quarterly')
    
    //s_agg
    theta_z     ${\theta_z}$ (long_name='technology shock')
    theta_b     ${\theta_b}$ (long_name='discount factor shock')
    theta_phi     ${\theta_{\phi}}$ (long_name='wage markup shock')
    theta_F     ${\theta_F}$ (long_name='shock to vacancy supply')
    theta_A     ${\theta_A}$ (long_name='matching efficiency shocks')
    

    
    % Observation equations
    x      % labor productivity
    x_obs_m % monthly growth rate of labor productivity
    x_obs % quarterly growth rate of labor productivity
    w_obs_m % montly growth rate of real wage
    w_obs % quarterly growth rate of real wage

    u_obs % demeaned unemployment rate
    v_obs % demeaned vacancy rate
    //s_obs % demeaned aggregate separation rate

    u_imp
    v_imp
    theta_imp
    e_imp
    N_imp
    N_e_imp
    C_imp;
    

varexo e_z ${e_z}$ (long_name= 'technology shock')
       e_b ${e_{b}}$ (long_name= 'discount factor shock')
       e_phi ${e_{phi,L}}$ (long_name= 'wage markup shock')
       e_F   ${e_F}$ (long_name = 'vacancy supply shock')
       e_A   ${e_F}$ (long_name = 'matching efficiency shock')
    ;
    
parameters 
% A, F,phi, and s are dependent parameters and are thus commented out
    z     ${x}$        (long_name='Steady-state productivity')
    beta  ${\beta}$    (long_name='Discount factor')
    delta ${\delta}$   (long_name='Product destruction rate')
    mu_net  ${\mu-1}$ (long_name='net markup')
    sigma  ${\sigma}$ (long_name='inverse of intertemporal elasticity of substitution')
    //s ${S}$            (long_name='Worker separation rate')
    tau ${\tau}$       (long_name='Aggregate separation rate')

    fbar
    qbar ${\overline{q}}$ (long_name ='Mean vacancy filling rate')

    //F ${F}$            (long_name = 'Level parameter in vacancy creation')
    xi_inv  ${1/\xi}$        (long_name = 'Inverse elasticity of entrants to vacancy value')
    b_ratio $b$              (long_name = 'Replacement ratio: unemployment benefits to wages')
    //A $A$              (long_name='Matching function level parameter')
    eta_L $\eta_L$         (long_name='Elasticity of matching function with respect to unemployment')
    
    labor_share // replaces phi
    N_ss ${N_{ss}}$ // replaces f_e

    rho_z    ${\rho_z}$               (long_name='persistence technology shock')
    rho_b    ${\rho_b}$               (long_name='persistence discount-factor shock')
    rho_phi    ${\rho_phi}$               (long_name='persistence wage markup shock')
    rho_F    ${\rho_F}$               (long_name='persistence cost vacancy supply shock')
    rho_A    ${\rho_A}$               (long_name='persistence cost vacancy supply shock')
    ;

%----------------------------------------------------------------
% set parameter values 
%----------------------------------------------------------------
z = 1.0;
beta = 0.99673; % monthly discount factor
xi_inv = 0.5; % 0 corresponds to infinitely elastic vacancy creation--standard DMP free entry
//r_ann = 0.04; % annual interest rate 

delta = 0.00514;
tau = 0.031;
//fbar = 0.41;
fbar = 0.41;
qbar = 0.8;

labor_share = 0.66;
N_ss = 1.0;

mu_net = 0.3;
sigma = 1.5;
b_ratio = 0.71;
eta_L = 0.6;


rho_z = 0.979;
rho_b = 0.979;
rho_phi = 0.979;
rho_F = 0.979;
rho_A = 0.979;

%----------------------------------------------------------------
% enter model equations
%----------------------------------------------------------------

model;

% Dependent parameters
# rho=(1-beta)/beta;
# xi = 1/xi_inv;

% Correct job finding and vacancy finding probabilities
#f_ss = fbar/(1-delta);
#q_ss = qbar/(1-delta);
#theta_ss = f_ss/q_ss;
#u_ss = tau/(tau+(1-delta)*f_ss);
#L_ss = 1-u_ss;

% Vacancies
#v_ss = theta_ss*u_ss;
%
#A = f_ss/(theta_ss^(1-eta_L));
% Separation rate
#s = (tau-delta)/(1-delta);

% Entrants
#e_ss = delta*(v_ss+1-u_ss);

#epsi = (mu_net+1)/mu_net;
#mu = epsi/(epsi-1);
#p_ss = N_ss^(1/(epsi-1));
#N_e_ss = delta/(1-delta)*N_ss;
#w_R_ss = p*z/mu;

#recruiter_share= (delta+(rho+delta)*(epsi-1))/(delta+(rho+delta)*epsi);
#w_ss = w_R_ss*labor_share/recruiter_share;

#b = b_ratio*w_ss;

#surplus_ratio = (rho+tau)/(1-delta)*(1/q);
#K_ss = (w_R_ss-w_ss)/(1+surplus_ratio);

// From wage equation find ϕ
#phi = (w_ss-b)/(w_R_ss-K_ss+theta_ss*K_ss/(1-delta)-b);

#Q_ss = K_ss*(1+rho)/(rho+delta);
#F = e_ss/(Q_ss^(1/xi_inv));

#f_e = (mu-1)*(L_ss/N_ss)*(1-delta)/(delta*mu+rho);
    
[name= 'Job creation condition']
K/q = beta*exp(theta_b)*lam(+1)/lam*(1-delta)*(w_R(+1)-w(+1)-K(+1)+(1-s)*K(+1)/q(+1));

[name ='Marginal revenue product']
w_R = p*exp(theta_z)/mu;

[name = 'Wage equation']
w = phi*exp(theta_phi)*(w_R - K + theta/(1-delta)*K) + (1-phi*exp(theta_phi))*b;

[name = 'Value of a vacancy']
Q = (e/F*exp(theta_F))^(xi_inv);

[name = 'Expected discounted difference in vacancy value']
K = (Q-beta*exp(theta_b)*lam(+1)/lam*(1-delta)*Q(+1));

[name = 'Market tightness']
theta = v/u;

[name = 'Job finding probability']
f = A*exp(theta_A)*theta^(1-eta_L);

[name = 'Vacancy filling probability']
q = A*exp(theta_A)*theta^(-eta_L);

[name='Aggregate labor']
L = 1-u;

[name = 'Composition of labor']
L = L_c + L_e;

[name='Relative price']
p = N^(1/(epsi-1));

[name='Lagrangian multiplier']
lam = C^(-sigma);

[name='Firm value']
nu_f = beta*exp(theta_b)*(1-delta)*lam(+1)/lam*(nu_f(+1)+d_f(+1));

[name ='Retail output: resources']
Y_c = p*exp(theta_z)*L_c;

[name ='Retail output: expenditure']
Y_c = C + Omega;

[name = 'Sunk vacancy posting costs']
Omega = F*exp(theta_F)/(1+xi_inv)*(e/(F*exp(theta_F)))^(1+xi_inv);

[Name = 'Business entrants']
N_e = L_e*exp(theta_z)/f_e;

[Name = 'Firm value relative to price']
nu_f = p*f_e/mu;

[Name= 'Expenditure = income']
Y_c + nu_f*N_e = w_R*L + N*d_f;

[Name = 'Output']
Y = Y_c+nu_f*N_e;

[Name = 'Data-consistent output']
Y_R = Y/p;

[name = 'LOM of vacancies']
v = (1-delta)*((1-q(-1))*v(-1) + s*(1-u(-1))) + e;

[name = 'LOM of unemployment']
u = (1-(1-delta)*f(-1))*u(-1) + tau*(1-u(-1)); % u is predetermined

[name = 'LOM of firms']
N = (1-delta)*(N(-1) + N_e(-1)); % N is predetermined

[name = 'labor productivity']
x = Y_R/L;

[name = 'labor share']
labor_share_obs = w*L/Y;

[name = 'consumption share']
cons_share_obs_m = C/Y;

cons_share_obs = (cons_share_obs_m + cons_share_obs_m(-1) + cons_share_obs_m(-2))/3;

[name = 'Aggregate separation rate']
//s_agg = (1-(1-delta*exp(theta_delta))*(1-s));

% Observation variables: first differences (demeaned) -> link to data in first differences (p. 58 of Pfeifer's Observation Equations)
//z_obs_m = theta_z - theta_z(-1); % Monthly growth rate of productivity (Delta \log x_t = theta_xt - theta_{x,t-1})
x_obs_m = log(x) - log(x(-1));
w_obs_m = log(w/p) - log(w(-1)/p(-1)); % Monthly growth rate of data-consistent real wages

% Quarterly growth rate of labor productivity and wages is sum of past 3 months
x_obs = x_obs_m + x_obs_m(-1) + x_obs_m(-2);
w_obs = w_obs_m + w_obs_m(-1) + w_obs_m(-2);

u_obs = log(u/u(-1)); % Demeaned unemployment rate
//u_obs = u - steady_state(u);
v_obs = log(v/v(-1)); % Demeaned vacancy rate
//v_obs = v - steady_state(v);
//s_obs = s_agg - tau;

% Variables for impulse responses (% deviation from steady state)
u_imp = 100*log(u);
v_imp = 100*log(v);
theta_imp = 100*log(theta);
e_imp = 100*log(e);
N_imp = 100*log(N);
N_e_imp = 100*log(N_e);
C_imp   = 100*log(C);
//Equivalent of x_imp is theta_x

% Exogenous processes
[name ='Technology shock process']
theta_z = rho_z*theta_z(-1) - e_z;

[name ='Discount factor process']
theta_b = rho_b*theta_b(-1) + e_b;

[name ='Bargaining power process']
theta_phi = rho_phi*theta_phi(-1) + e_phi;

[name ='Vacancy supply process']
theta_F = rho_F*theta_F(-1) + e_F;

[name ='Matching efficiency process']
theta_A = rho_A*theta_A(-1) + e_A;

end;

steady_state_model;
    //Do Calibration
    //Steady state imposes fbar, qbar, r_ss, tau
    N = 1;
    epsi_ss = (mu_net+1)/mu_net;
    mu_ss = epsi_ss/(epsi_ss-1);

    p = N^(1/(epsi_ss-1));
    w_R = p*z/mu_ss;
    N_e = delta/(1-delta)*N;

    rho_ss = (1-beta)/beta;
    f = fbar/(1-delta);
    q = qbar/(1-delta);

    theta = f/q;
    u = tau/(tau+(1-delta)*f);
    L = 1 - u;
    v = theta*u;
    e = delta*(v+1-u);

    recruiter_share_ss= (delta+(rho_ss+delta)*(epsi_ss-1))/(delta+(rho_ss+delta)*epsi_ss);
    w = w_R*labor_share/recruiter_share_ss;

    surplus_ratio_ss = (rho_ss+tau)/(1-delta)*(1/q);
    K = (w_R-w)/(1+surplus_ratio_ss);

    //K = (rho_ss+delta)/(1+rho_ss)*(e/F)^(xi_inv);

    //b = b_ratio*w;
    //w = phi*(w_R-K+theta/(1-delta)*K) + (1-phi)*b;

    f_e_ss = (mu_ss-1)*(L/N)*(1-delta)/(delta*mu_ss+rho_ss);

    nu_f = p*f_e_ss/mu_ss;
    d_f = (rho_ss+delta)/(1-delta)*nu_f;

    L_e = N_e*f_e_ss/z;
    L_c = L - L_e;
    Y_c = p*z*L_c;
    Q = K*(1+rho_ss)/(rho_ss+delta);

    F_ss = e/(Q^(1/xi_inv));

    Omega = F_ss/(1+xi_inv)*(e/F_ss)^(1+xi_inv);
    C = Y_c - Omega;
    lam = C^(-sigma);
    Y = (Y_c + nu_f*N_e);
    Y_R = Y/p;
    x = Y_R/L;

    labor_share_obs = labor_share;
    cons_share_obs_m = C/Y;
    cons_share_obs = C/Y;

    theta_z = 0;
    theta_b = 0;
    theta_phi = 0;
    theta_F = 0;
    theta_A = 0;

    x_obs = 0;
    x_obs_m = 0;
    w_obs_m = 0;
    w_obs = 0;

    u_obs = 0;
    //s_obs = 0;
    v_obs = 0;

    u_imp = 100*(log(u));
    v_imp = 100*(log(v));
    theta_imp = 100*(log(theta));
    e_imp = 100*(log(e));
    N_imp = 100*(log(N));
    N_e_imp = 100*(log(N_e));
    C_imp   = 100*(log(C));

end;

//set shock variances
shocks;
    var e_z=0.007^2;
    var e_b = 0.0072^2;
    var e_phi=0.0072^2;
end;

// local identification
//identification(ar=1);
//check the starting values for the steady state
resid;

// compute steady state given the starting values
steady;
// check Blanchard-Kahn-conditions
check;


estimated_params;
//x, init_value, lower bound, upper bound, prior shape, prior mean, prior std

//phi, 0.5, 0.0,0.99,             BETA_PDF, 0.5, 0.25;
b_ratio, 0.71, 0.4, 0.95,         BETA_PDF, 0.71, 0.1; %b_ratio = b/w
eta_L, 0.6, 0.5, 0.8,             BETA_PDF, 0.6, 0.2;
delta, 0.00514, 0.0, tau,         BETA_PDF, 0.00514, 0.001, 0.0, tau;
xi_inv, 1.0, 0.0, 20,             GAMMA_PDF, 1.0, 0.5;
//epsi, 4.0, 1.1, 100,              GAMMA_PDF, 4.3, 1.0;
mu_net, 0.3, 0.0, 2.0,            GAMMA_PDF, 0.3, 0.1;
sigma, 1.5, 1.0, 4.0,             GAMMA_PDF, 1.5, 0.5;


% Persistence parameters
rho_z,  0.9, 0.0001, 0.99999,       BETA_PDF, 0.6, 0.2;
rho_b, 0.9, 0.01, 0.999999,         BETA_PDF, 0.6, 0.2;
rho_phi,  0.9, 0.01, 0.999999,      BETA_PDF, 0.6, 0.2;
rho_F,  0.9, 0.01, 0.999999,      BETA_PDF, 0.6, 0.2;
rho_A,  0.9, 0.01, 0.999999,      BETA_PDF, 0.6, 0.2;

% Standard errors
stderr e_z, 0.01, 0.0000001, 0.2,      INV_GAMMA_PDF, 0.01, 0.1;
stderr e_b, 0.01, 0.00001, 0.2,        INV_GAMMA_PDF, 0.01, 0.1;
stderr e_phi, 0.01, 0.00001, 0.2,      INV_GAMMA_PDF, 0.01, 0.1;
stderr e_F, 0.01, 0.00001, 0.2,      INV_GAMMA_PDF, 0.01, 0.1;
stderr e_A, 0.01, 0.00001, 0.2,      INV_GAMMA_PDF, 0.01, 0.1;

end;

options_.TeX=1;
% Estimate on monthly unemployment rate, monthly vacancy rate, quarterly labor productivity growth

varobs u_obs, v_obs, x_obs, w_obs;
//varobs u_obs, x_obs;
estimation(tex, optim=('MaxIter', 200), 
datafile=observables_u_growth, 
//mode_file=baseline_est_mode, %With _mh option uses mode after MCM run
nograph,
//load_mh_file, 
//mh_recover,
mcmc_jumping_covariance=prior_variance,

mode_compute=4,
presample=0, 
lik_init=2,
mh_jscale=0.001, 
mh_init_scale =0.0001,
//mh_jscale=0.1,
mode_check, 
mh_replic=100000, 
//mh_replic=0,
mh_nblocks=2, 
//bayesian_irf,
//irf=100,
mh_drop=0.3, 
//moments_varendo,
prior_trunc=0)
u_obs, v_obs, x_obs, w_obs, cons_share_obs, theta_z, u_imp, v_imp, theta_imp, e_imp, N_imp, N_e_imp, C_imp, labor_share_obs, L_c, L_e;


% Stochastic simulation -> for conditional FEVD and IRF
stoch_simul (order=1, nofunctions, irf=80, periods=0)
//conditional_variance_decomposition=[1 4 8 40])
u_obs, v_obs, x_obs, w_obs, cons_share_obs, theta_z, 
u_imp, v_imp, theta_imp, e_imp, N_imp, N_e_imp, C_imp, labor_share_obs, L_c, L_e;

