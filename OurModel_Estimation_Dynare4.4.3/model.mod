///////////////////////////////////////////////////////////////////////////
//
// Medium-scale DSGE model with shopping friction and monetary policy
//
// Key Features:
//   - New Keynesian framework with price and wage rigidities
//   - Shopping friction mechanism affecting consumption
//   - Three exogenous shocks:
//       1. Monetary policy shock
//       2. Neutral technology shock
//       3. Investment-specific technology shock
//   - Bayesian estimation using impulse response matching
//   - Requires modified Dynare functions for IRF matching
//
// Version History:
//   Original version created by Zhesheng Qiu at 08:19 2019-01-05
//   Multiple updates by Zhesheng Qiu through 18:40 2025-05-27
//
// Note:
//   All Variables have been denominated (if not mentioned) and detrended.
//   The syntax is consistent with Dynare 4.4.3.
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
// Define variables（in the same order as equations
///////////////////////////////////////////////////////////////////////////

var
// Group 1: search friction
    theta       ${theta}$       (long_name='Market tightness')
    f           ${f}$           (long_name='Job finding rate')
    q           ${q}$           (long_name='Vacancy filling rate')
// Group 2: relative price
    p           ${p}$           (long_name='Relative price')
// Group 3: business formation
    lam         ${\lambda}$     (long_name='Marginal utility of consumption')
    nu_f        ${nu_f}$        (long_name='Retailer firm value')
    d_f         ${d_f}$         (long_name='Dividend')
    N_e         ${N_e}$         (long_name='Number of new businesses')
// Group 4: vacancy creation
    K           ${K}$           (long_name='Expected discounted difference in vacancy value')
    Q           ${Q}$           (long_name='Value of vacancy')
    e           ${e}$           (long_name='New vacancies')
    w_R         ${w_R}$         (long_name='Recruiter payment')
    w           ${w}$           (long_name='Wage')
// Group 5: law of motion
    N           ${N}$           (long_name='Number of firms')
    u           ${u}$           (long_name='Unemployment')
    v           ${v}$           (long_name='Vacancy stock')
// Group 6: aggregate identity
    L           ${L}$           (long_name='Total labor')
    L_e         ${L_e}$         (long_name='Labor in entry')
    L_c         ${L_c}$         (long_name='Labor in consumption')
    Omega       ${\Omega}$      (long_name='Sunk entry costs of investment')
    C           ${C}$           (long_name='Consumption')
    Y_c         ${Y_c}$         (long_name='Retail output')
    Y           ${Y}$           (long_name='Aggregate output')
    ls          ${\ell s}$      (long_name='Labor share')
    s_agg       ${s_{agg}$      (long_name='Aggregate separation rate')
// Group 7: exogenous process
    log_z       ${log(z)}$      (long_name='Tech shocks')
    log_delta   ${log(\delta)}$ (long_name='Product destruction shocks')
    log_s       ${log(s)}$      (long_name='Separation shocks')
// Group 8: observables
    u_obs       
    v_obs      
    theta_obs   
    lp_obs    
    s_obs      
    bf_obs   
    u_obs_lag   
    v_obs_lag  
    theta_obs_lag
    lp_obs_lag   
    s_obs_lag   
    bf_obs_lag
;

// Define exogenous shocks
varexo 
    e_z         ${e_z}$         (long_name='Technology shock')
    e_delta     ${e_\delta}$    (long_name='Firm destruction shock')
    e_s         ${e_s}$         (long_name='Separation shock');


///////////////////////////////////////////////////////////////////////////
// Model parameters
///////////////////////////////////////////////////////////////////////////

parameters
    // Labor market parameters
    b_ratio     ${b}$           (long_name='Replacement ratio: unemployment benefits to wages')
    x_v         ${x_v}$         (long_name='Share of hiring costs from vacancy creation')
    xi_inv      ${1/\xi}$       (long_name='Inverse elasticity of entrants to vacancy value')
    // Entry and monopolistic competition
    delta       ${\delta}$      (long_name='Product destruction rate')
    epsi        ${\varepsilon}$ (long_name='Elasticity of substitution')
    % Shock process parameters
    rho_z       ${\rho_z}$      (long_name='Persistence of tech shocks')
    rho_delta   ${\rho_delta}$  (long_name='Persistence of product destruction shocks')
    rho_s       ${\rho_s}$      (long_name='Persistence of separation shocks')
    sigma_z     ${\sigma_z}$    (long_name='Std. dev of tech shocks')
    sigma_delta ${\sigma_delta}$(long_name='Std. dev of prod. dest. shocks')
    sigma_s     ${\sigma_s}$    (long_name='std. dev of separation shocks')
    // Fixed parameters
    beta        ${\beta}$       (long_name='Discount factor')
    eta_L       ${\eta_L}$      (long_name='Elasticity of matching function with respect to unemployment')
    sigma       ${\sigma}$      (long_name='Inverse of intertemporal elasticity of substitution')
    tau         ${\tau}$        (long_name='Aggregate separation rate')
    // Implied parameters from the steady state
    A           ${A}$           (long_name='Match efficiency')
    f_e         ${f_e}$         (long_name='Cost of business creation')
    kappa       ${\kappa}$      (long_name='Cost of job posting')
    F           ${F}$           (long_name='Parameter of vacancy creation')
    phi         ${\phi}$        (long_name='Bargaining power')
;


///////////////////////////////////////////////////////////////////////////
// Model equations
///////////////////////////////////////////////////////////////////////////

model;

// Group 1: search friction

    [name='Market tightness']
    theta = v/u;
    
    [name='Job finding probability']
    f = A*theta^(1-eta_L);
    
    [name='Vacancy filling probability']
    q = A*theta^(-eta_L);

// Group 2: relative price

    [name='Relative price']
    p = N^(1/(epsi-1));

// Group 3: business formation

    [name='Lagrangian multiplier']
    lam = C^(-sigma);
    
    [name='Firm value relative to price']
    nu_f = p*f_e*(epsi-1)/epsi;

    [name='Firm value']
    nu_f = beta*(1-exp(log_delta))*lam(+1)/lam*(nu_f(+1)+d_f(+1));

    [name='Business entrants']
    N_e = L_e*exp(log_z)/f_e;

// Group 4: Vacancy creation

    [name='Job creation condition']
    (kappa+K/q) = beta*lam(+1)/lam*(1-exp(log_delta))*(w_R(+1)-w(+1)-K(+1)+(1-exp(log_s))*(kappa+K(+1)/q(+1)));
    
    [name='Expected discounted difference in vacancy value']
    K = (Q-beta*lam(+1)/lam*(1-exp(log_delta))*Q(+1));
    
    [name='Value of a vacancy']
    Q = (e/F)^(xi_inv);
    
    [name='Marginal revenue product']
    w_R = p*exp(log_z)*(epsi-1)/epsi;
    
    [name='Wage equation']
    w = phi*(w_R - K + theta/(1-exp(log_delta))*(K+q*kappa)) + (1-phi)*b_ratio*steady_state(w);
    
// Group 5: laws of motion
    
    [name='LOM of firms']
    N = (1-exp(log_delta))*(N(-1) + N_e(-1));
    
    [name='LOM of unemployment']
    u = (1-(1-exp(log_delta))*f(-1))*u(-1) + s_agg*(1-u(-1));
    
    [name='LOM of vacancies']
    v = (1-exp(log_delta))*((1-q(-1))*v(-1) + exp(log_s)*(1-u(-1))) + e;

// Group 6: aggregate identity
    
    [name='Aggregate labor']
    L = 1-u;
    
    [name='Composition of labor']
    L = L_c + L_e;    
    
    [name='Retail output: resources']
    Y_c = p*exp(log_z)*L_c;
    
    [name='Sunk vacancy posting costs']
    Omega = F/(1+xi_inv)*(e/F)^(1+xi_inv);
    
    [name='Retail output: expenditure']
    Y_c = C+Omega+kappa*v*q;
    
    [name='Output: expenditure']
    Y = Y_c + nu_f*N_e;
    
    [name='Output: income']
    Y = w_R*L + N*d_f;

    [name='Labor share']
    ls = w*L/Y;

    [name='Aggregate separation rate']
    s_agg = (1-(1-exp(log_delta))*(1-exp(log_s)));

// Group 7: exogenous process

    [name='Labor productivity process']
    log_z = rho_z*log_z(-1) + (1-rho_z)*steady_state(log_z) - sigma_z*e_z;

    [name='Product line destruction process']
    log_delta = rho_delta*log_delta(-1) + (1-rho_delta)*steady_state(log_delta) + sigma_delta*e_delta;

    [name='Separation rate process']
    log_s = rho_s*log_s(-1) + (1-rho_s)*steady_state(log_s) - sigma_s*e_s;

// Group 8: observables

    [name='Unemployment observation']
    u_obs = 1/3*(u+u(-1)+u(-2));

    [name='Vacancy observation']
    v_obs = 1/3*(v+v(-1)+v(-2));

    [name='Tightness observation']
    theta_obs = v_obs/u_obs;

    [name='Labor productivity observation']
    lp_obs = 1/3*(exp(log_z)+exp(log_z(-1))+exp(log_z(-2)));

    [name='Job separation observation'] //
    s_obs = 1/3*(exp(log_s)+exp(log_s(-1))+exp(log_s(-2)));

    [name='Business formation observation']
    bf_obs = (N_e+N_e(-1)+N_e(-2))/steady_state(N);

// Group 9: lag observables

    [name='Lag unemployment observation']
    u_obs_lag = u_obs(-3);

    [name='Lag vacancy observation']
    v_obs_lag = v_obs(-3);

    [name='Lag tightness observation']
    theta_obs_lag = theta_obs(-3);

    [name='Lag labor productivity observation']
    lp_obs_lag = lp_obs(-3);

    [name='Lag job separation observation'] //
    s_obs_lag = s_obs(-3);

    [name='Lag business formation observation']
    bf_obs_lag = bf_obs(-3);

end;
    

///////////////////////////////////////////////////////////////////////////

resid;
steady;
check;

shocks;
var e_z = 1;
var e_delta = 1;
var e_s = 1;
end;


///////////////////////////////////////////////////////////////////////////
// Observed variables (irrelevant for impulse response matching
///////////////////////////////////////////////////////////////////////////

varobs Y;   // check out data.m. We just pass-on random data to prevent Dynare
            // from getting stuck here. VAR impulse response data is loaded and 
            // handled below.

// add folder with modified Dynare functions and utilities
addpath('DynareUtilites');


///////////////////////////////////////////////////////////////////////////
// Priors for estimation
///////////////////////////////////////////////////////////////////////////

estimated_params;
b_ratio,    beta_pdf,       0.71, 0.2;
x_v,        beta_pdf,       0.5, 0.25;
xi_inv,     gamma_pdf,      1.0, 2;
delta,      beta_pdf,       0.0083, 0.005; // 0.0083 monthly corresponds to 10% annual
epsi,       gamma_pdf,      4.2, 1.5;
rho_z,      beta_pdf,       0.8, 0.2;
rho_delta,  beta_pdf,       0.93, 0.1;
rho_s,      beta_pdf,       0.8, 0.2;
sigma_z,    inv_gamma_pdf,  0.01, 1.0;
sigma_delta,inv_gamma_pdf,  0.02, 0.01; // prior can be influenced by product destruction data, not in dataset
sigma_s,    inv_gamma_pdf,  0.01, 1.0;
end;


///////////////////////////////////////////////////////////////////////////

// Begin of estimation

do_imp_resp_matching=1;                  %if =1, model estimated based on VAR impulse responses
if do_imp_resp_matching==1,
   addpath('DynareImpRespMatching');     %path with adjusted DYNARE functions dsge_likelihood.m and dynare_estimation_1.m  
   get_empirical_moments;                %get VAR impulse responses: load VAR impulse responses and precision matrix
                                         %from .mat file; store them in global variables psihat and inv_Vhat
                                         %The model-VAR impulse response matching is done in dsge_likelihood.m (in folder DynareImpRespMatching). 
                                         %Check out lines 822-950 for the details. 
end;

% The dsge_likelihood.m needs to be modified for our own model.


///////////////////////////////////////////////////////////////////////////

// Begin of estimation

options_.mh_conf_sig    = 0.95;
options_.prior_interval = 0.95;
options_.conf_sig       = 0.95;

estimation(first_obs=1,
           datafile=dummy_data,
           nograph,
           order=1,
           lik_init=1,
           prior_trunc=0,
           mode_compute=6,
           mode_file=model_mh_mode,
           mh_replic=100000,%910000,%1500000,%720000,           
           mh_nblocks=1,%11,
           mh_init_scale=0.5,
           mh_jscale=0.5,
           mh_drop=0.50);
    
// evaluate at MCMC joint mode

load model_mh_mode.mat; % need to reload this for each estimation?
for ii=1:1:length(xparam1)
    para_idx=loc(M_.param_names,char(parameter_names(ii)));
    M_.params(para_idx)=xparam1(ii);
end

// call steady state
model_steadystate([],[]);

// Stochastic simulation
stoch_simul(order=1,nofunctions,nograph,noprint)
    u_obs, v_obs, theta_obs, lp_obs, s_obs, bf_obs;
stack_model_moments;
save results;
delete *.eps;