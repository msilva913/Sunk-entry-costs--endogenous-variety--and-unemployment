

// define variables

var //x           ${x}$ (long_name='technology')
    w        ${w}$ (long_name = 'wage')
    v         ${v}$ (long_name = 'vacancy stock')
    u         ${u}$ (long_name = 'unemployment')
    theta     ${theta}$ (long_name='market tightness')
    f         ${f}$ (long_name = 'job finding rate')
    q         ${q}$ (long_name = 'vacancy filling rate')
    
    theta_x     ${\theta_x}$ (long_name='technology shock')


    u_imp
    v_imp
    theta_imp
;
    

varexo e_x ${e_x}$ (long_name= 'technology shock')
    ;
    
parameters 
    x     ${x}$        (long_name='Steady-state productivity')
    beta  ${\beta}$    (long_name='Discount factor')
    delta ${\delta}$   (long_name='Product destruction rate')
    s ${s}$            (long_name='Worker separation rate')
    fbar
    qbar ${\overline{q}}$ (long_name ='Mean vacancy filling rate')
    //k ${k}$            (long_name = 'Vacancy posting cost')
    alpha_L $\alpha_L$ (long_name = 'Worker bargaining power')
    b $b$              (long_name = 'Unemployment value')
    //A $A$              (long_name='Matching function level parameter')
    phi $\phi$         (long_name='Elasticity of matching function with respect to unemployment')



    rho_x    ${\rho_x}$               (long_name='persistence technology shock')
    rho_b   ${\rho_Z}$                (long_name='persistence discount factor shock')
    rho_alphaL    ${\rho_{\alpha_L}}$ (long_name='persistence bargaining power shock')
    //rho_delta    ${\rho_N}$           (long_name='persistence product destruction rate shock')
    ;

%----------------------------------------------------------------
% set parameter values 
%----------------------------------------------------------------
x = 1.0;
beta = exp(-4/1200); % monthly discount factor
alpha_L = 0.5;
b=0.92;
phi=0.5;
s = 0.035;
fbar = 1/2.2;
qbar = 1 - (1-1/3)^4;


rho_x = 0.979;
stdx = 0.007;

%----------------------------------------------------------------
% enter model equations
%----------------------------------------------------------------

model;

% Dependent parameters
# r=(1-beta)/beta;

#theta_ss = fbar/qbar;
#A =fbar/theta_ss^(1-phi);
#k = (1-alpha_L)*(x-b)/((r+s)/qbar+alpha_L*theta_ss);

[name= 'Job creation condition']
k/q = beta*(exp(theta_x(+1))-w(+1)+(1-s)*k/q(+1));

[name = 'Wage equation']
w = alpha_L*(exp(theta_x)+k*theta) + (1-alpha_L)*b;

[name = 'LOM of unemployment']
u = (1-f(-1))*u(-1) + s*(1-u(-1));

[name = 'Market tightness']
theta = v/u;

[name = 'Job finding probability']
f = A*theta^(1-phi);

[name = 'Vacancy filling probability']
q = A*theta^(-phi);


% Variables for impulse responses (% deviation from steady state)
u_imp = 100*log(u/steady_state(u));
v_imp = 100*log(v/steady_state(v));
theta_imp = 100*log(theta/steady_state(theta));
//Equivalent of x_imp is theta_x

% Exogenous processes
[name ='Labor productivity process']
theta_x = rho_x*theta_x(-1) - e_x;

end;

steady_state_model;
    //Do Calibration
    //Steady state imposes fbar, qbar, r_ss, tau
    
    r_ss = (1-beta)/beta;
    f = fbar;
    q = qbar;
    theta = f/q;
    u = s/(s+f);
    v = theta*u;

    v = theta*u;
    k_ss = (1-alpha_L)*(x-b)/((r_ss+s)/qbar+alpha_L*theta);
    w = alpha_L*(x+k_ss*theta) + (1-alpha_L)*b;
   
    theta_x = 0;
  
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

xi, 0.85, 0.5, 2.0,        GAMMA_PDF, 0.85, 0.1;


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
varobs u_obs, v_obs, x_obs;


estimation(tex, optim=('MaxIter', 200), 
datafile=observables_sectoral, 
mode_file=BRS_sectoral_mh_mode, %With _mh option uses mode after MCM run
//nograph,
load_mh_file, 
//mh_recover,
mcmc_jumping_covariance=prior_variance,

mode_compute=0,
presample=0, 
lik_init=2,
mh_jscale=0.006, 
mh_init_scale =0.0001,
//mh_jscale=0.1,
mode_check, 
//mh_replic=75000, 
mh_replic=0,
mh_nblocks=2, 
//bayesian_irf,
//irf=100,
mh_drop=0.3, 
//moments_varendo,
prior_trunc=0)
Y_obs, Y_N_obs, I_obs, p_I_obs, C_obs, NC_obs, NI_obs, util_ND_obs, util_D_obs, SR_obs, util_obs, D_obs, h_obs;
//log_Y, log_Y_N, log_I, log_p_I, log_C, log_N, log_NC, log_NI, util;



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
stoch_simul (order=2, nofunctions, irf=80, periods=0, pruning)
//conditional_variance_decomposition=[1 4 8 40])
theta_x, u_imp, v_imp, theta_imp;

