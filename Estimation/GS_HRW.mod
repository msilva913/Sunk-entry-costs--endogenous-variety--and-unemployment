

// define variables

var //x           ${x}$ (long_name='technology')
    K           ${K}$ (long_name='expected discounted difference in vacancy value')
    w        ${w}$ (long_name = 'wage')
    e           ${e}$ (long_name = 'new vacancies')
    v         ${v}$ (long_name = 'vacancy stock')
    u         ${u}$ (long_name = 'unemployment')
    theta     ${theta}$ (long_name='market tightness')
    f         ${f}$ (long_name = 'job finding rate')
    q         ${q}$ (long_name = 'vacancy filling rate')
    
    theta_x     ${\theta_x}$ (long_name='technology shock')
    theta_b     ${\theta_b}$ (long_name='discount factor shock')
    theta_alphaL     ${\theta_{\alpha_L}}$ (long_name='bargaining power shock')

    
    % Observation equations
    x_obs_m % monthly growth rate of labor productivity
    x_obs % quarterly growth rate of labor productivity

    w_obs_m % monthly growth rate of real wages
    w_obs   % quarterly growth rate of real wages
    u_obs % demeaned unemployment rate
    v_obs % demeaned vacancy rate


    u_imp
    v_imp
    theta_imp
    e_imp;
    //predetermined_variables u, u_obs, u_imp; % adjust to beginning-of-stock convention
    

varexo e_x ${e_x}$ (long_name= 'technology shock')
       e_b ${e_{b}}$ (long_name= 'discount factor shock')
       e_alphaL ${e_{alpha,L}}$ (long_name= 'bargaining power shock')
    ;
    
parameters 
% A, F, and s are dependent parameters and are thus commented out
    x     ${x}$        (long_name='Steady-state productivity')
    beta  ${\beta}$    (long_name='Discount factor')
    delta ${\delta}$   (long_name='Product destruction rate')
    //s ${S}$            (long_name='Worker separation rate')
    tau ${\tau}$       (long_name='Aggregate separation rate')
    fbar
    qbar ${\overline{q}}$ (long_name ='Mean vacancy filling rate')
    //F ${F}$            (long_name = 'Level parameter in vacancy creation')
    xi_inv  ${1/\xi}$        (long_name = 'Inverse elasticity of entrants to vacancy value')
    //xi  ${\xi}$        (long_name = 'Elasticity of entrants to vacancy value')
    alpha_L $\alpha_L$ (long_name = 'Worker bargaining power')
    b $b$              (long_name = 'Unemployment value')
    //A $A$              (long_name='Matching function level parameter')
    nu_L $\phi$         (long_name='Elasticity parameter of HRW')



    rho_x    ${\rho_x}$               (long_name='persistence technology shock')
    rho_b   ${\rho_Z}$                (long_name='persistence discount factor shock')
    rho_alphaL    ${\rho_{\alpha_L}}$ (long_name='persistence bargaining power shock')
    ;

%----------------------------------------------------------------
% set parameter values 
%----------------------------------------------------------------
x = 1.0;
beta = 0.99673; % monthly discount factor
xi_inv = 1.0; % 0 corresponds to infinitely elastic vacancy creation--standard DMP free entry
//r_ann = 0.04; % annual interest rate 
alpha_L = 0.566;
b=0.9;
nu_L = 1.5857;
delta = 0.00514;
tau = 0.034;
//fbar = 0.41;
fbar = 1/2.2;
qbar = 1 - (1-1/3)^4;


rho_x = 0.979;
rho_b = 0.979;
rho_alphaL = 0.979;

%----------------------------------------------------------------
% enter model equations
%----------------------------------------------------------------

model;

% Dependent parameters
# r=(1-beta)/beta;
# xi = 1/xi_inv;

% Correct job finding and vacancy finding probabilities
#f_ss = fbar/(1-delta);
#q_ss = qbar/(1-delta);
#theta_ss = f_ss/q_ss;
#u_ss = tau/(tau+(1-delta)*f_ss);

% Vacancies
#v_ss = theta_ss*u_ss;
%
//#A = f_ss/(theta_ss^(1-phi));
% Separation rate
#s = (tau-delta)/(1-delta);

% Entrants
#e_ss = v_ss - (1-delta)*((1-q_ss)*v_ss+s*(1-u_ss));

#K_ss = q_ss*(1-delta)*(1-alpha_L)*(x-b)/(r+tau+alpha_L*f_ss+q_ss*(1-delta)*(1-alpha_L));
    
#F = e_ss*((r+delta)/((1+r)*K_ss))^xi;

[name= 'Job creation condition']
K/q = beta*exp(theta_b)*(1-delta)*(exp(theta_x(+1))-w(+1)-K(+1)+(1-s)*K(+1)/q(+1));

[name = 'Expected discounted difference in vacancy value']
K = (e^(1/xi) - beta*exp(theta_b)*(1-delta)*e(+1)^(1/xi))/(F^(1/xi));

[name = 'Wage equation']
w = alpha_L*exp(theta_alphaL)*(exp(theta_x) - K + theta/(1-delta)*K) + (1-alpha_L*exp(theta_alphaL))*b;

[name = 'LOM of vacancies']
v = (1-delta)*((1-q(-1))*v(-1) + s*(1-u(-1))) + e;

[name = 'LOM of unemployment']
u = (1-(1-delta)*f(-1))*u(-1) + tau*(1-u(-1));

[name = 'Market tightness']
theta = v/u;

[name = 'Job finding probability']
f = theta*q;

[name = 'Vacancy filling probability']
//q = A*theta^(-phi);
q = (1+theta^(nu_L))^(-1/nu_L);

% Observation variables: first differences (demeaned) -> link to data in first differences (p. 58 of Pfeifer's Observation Equations)
x_obs_m = theta_x - theta_x(-1); % Monthly growth rate of productivity (Delta \log x_t = theta_xt - theta_{x,t-1})

% Quarterly growth rate of labor productivity is sum of past 3 months
x_obs = x_obs_m + x_obs_m(-1) + x_obs_m(-2);

w_obs_m = log(w/w(-1)); % Monthly growth rate of real wages

% Quarterly growth rate of wages is sum of past 3 months
w_obs = w_obs_m + w_obs_m(-1) + w_obs_m(-2);

%

//u_obs = log(u/u(-1)); % Growth rate of unemployment
u_obs = u - steady_state(u);
//v_obs = log(v/v(-1)); % Growth rate of vacancies
v_obs = v - steady_state(v);

% Variables for impulse responses (% deviation from steady state)
u_imp = 100*log(u/steady_state(u));
v_imp = 100*log(v/steady_state(v));
theta_imp = 100*log(theta/steady_state(theta));
e_imp = 100*log(e/steady_state(e));
//Equivalent of x_imp is theta_x

% Exogenous processes
[name ='Labor productivity process']
theta_x = rho_x*theta_x(-1) - e_x;

[name ='Discount factor process']
theta_b = rho_b*theta_b(-1) + e_b;

[name ='Bargaining power process']
theta_alphaL = rho_alphaL*theta_alphaL(-1) + e_alphaL;

end;

steady_state_model;
    //Do Calibration
    //Steady state imposes fbar, qbar, r_ss, tau
    
    r_ss = (1-beta)/beta;
    f = fbar/(1-delta);
    q = qbar/(1-delta);

    theta = f/q;
    u = tau/(tau+(1-delta)*f);
    v = theta*u;
    e = delta*(v+1-u);
    
    //K = (e/F)^(1/xi)*(r_ss+delta)/(1+r_ss);
    K = q*(1-delta)*(1-alpha_L)*(x-b)/(r_ss+tau+alpha_L*f+q*(1-delta)*(1-alpha_L));
    w = alpha_L*(x-K+theta/(1-delta)*K) + (1-alpha_L)*b;

    theta_x = 0;
    theta_b = 0;
    theta_alphaL = 0;

    x_obs_m = 0;
    x_obs = 0;
    w_obs_m = 0;
    w_obs = 0;
    u_obs = 0;
    v_obs = 0;

    x_imp = 0;
    v_imp = 0;
    theta_imp = 0;

end;

//set shock variances
shocks;
    var e_x=0.007^2;
    //var e_b = 0.0072^2;
    //var e_alphaL=0.0072^2;
end;

// local identification
//identification(ar=1);
//check the starting values for the steady state
resid;

// compute steady state given the starting values
steady;
// check Blanchard-Kahn-conditions
check;

/*
estimated_params;
//x, init_value, lower bound, upper bound, prior shape, prior mean, prior std

alpha_L, 0.5, 0.0,0.99,             BETA_PDF, 0.5, 0.25;
b, 0.71, 0.4, 0.95,                 BETA_PDF, 0.71, 0.1;
phi, 0.6, 0.5, 0.8,                 BETA_PDF, 0.6, 0.2;
delta, 0.00514, 0.0, tau,           BETA_PDF, 0.00514, 0.001, 0.0, tau;
xi_inv, 1.0, 0.0, 20,               GAMMA_PDF, 1.0, 0.5;

//xi, 0.85, 0.5, 2.0,        GAMMA_PDF, 0.85, 0.1;


% Persistence parameters
rho_x,  0.9, 0.0001, 0.99999,             BETA_PDF, 0.6, 0.2;
rho_b, 0.9, 0.01, 0.999999,          BETA_PDF, 0.6, 0.2;
rho_alphaL,  0.9, 0.01, 0.999999,    BETA_PDF, 0.6, 0.2;

% Standard errors
stderr e_x, 0.01, 0.0000001, 0.2,      INV_GAMMA_PDF, 0.01, 0.1;
stderr e_b, 0.01, 0.00001, 0.2,        INV_GAMMA_PDF, 0.01, 0.1;
stderr e_alphaL, 0.01, 0.00001, 0.2,   INV_GAMMA_PDF, 0.01, 0.1;

% news



end;

options_.TeX=1;
% Estimate on growth rates of monthly unemployment rate, monthly vacancy rate, quarterly labor productivity
varobs u_obs, x_obs, w_obs;
//varobs u_obs, x_obs;


estimation(tex, optim=('MaxIter', 200), 
datafile=observables_u_level, 
//mode_file=GS_mh_mode, %With _mh option uses mode after MCM run
nograph,
//load_mh_file, 
//mh_recover,
mcmc_jumping_covariance=prior_variance,

mode_compute=4,
presample=0, 
lik_init=2,
mh_jscale=0.006, 
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
theta_x, u_imp, v_imp, theta_imp, e_imp;



%----------------------------------------------------------------
% generate LaTeX output
%----------------------------------------------------------------

write_latex_dynamic_model;
write_latex_parameter_table;
write_latex_definitions;
write_latex_prior_table;
//generate_trace_plots(1);
collect_latex_files;
% if system(['pdflatex -halt-on-error -interaction=batchmode ' M_.fname '_TeX_binder.tex'])
%     error('TeX-File did not compile.')
% end
*/

% Stochastic simulation -> for conditional FEVD and IRF
stoch_simul (order=1, pruning, nofunctions, irf=80, periods=0)
//conditional_variance_decomposition=[1 4 8 40])
 theta_x, u_imp, v_imp, theta_imp, e_imp;

