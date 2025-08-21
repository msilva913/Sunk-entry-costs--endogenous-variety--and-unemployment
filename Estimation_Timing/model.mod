///////////////////////////////////////////////////////////////////////////
//
// The model with business formation and labor market frictions
//
// Key Features:
//   - Labor market search combined with business formation
//   - Vacancy cost removed from GDP defintion
//   - Distribution fo sunk cost normalized
//   - Three exogenous shocks:
//       1. Technology shock
//       2. Job separation shock
//       3. Firm destruction shock
//   - Bayesian estimation using method of simulated moments
//   - Requires modified Dynare functions for moment matching
//
// Version History:
//   Original version created by Zhesheng Qiu at xx:xx 2025-xx-xx
//   Multiple updates by Zhesheng Qiu through 18:24 2025-08-12
//
// Note:
//   The syntax is consistent with Dynare 4.4.3.
//
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
// Define variables（in the same order as equations)
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
    nu          ${nu}$          (long_name='Retailer firm value')
    d           ${d}$           (long_name='Dividend')
    Ne          ${N^{e}}$       (long_name='Number of new businesses')
// Group 4: vacancy creation
    K           ${K}$           (long_name='Expected discounted difference in vacancy value')
    Q           ${Q}$           (long_name='Value of vacancy')
    e           ${e}$           (long_name='New vacancies')
    wint        ${w^(int)}$     (long_name='Recruiter payment')
    w           ${w}$           (long_name='Wage')
// Group 5: law of motion
    N           ${N}$           (long_name='Number of firms')
    u           ${u}$           (long_name='Unemployment')
    v           ${v}$           (long_name='Vacancy stock')
    tau         ${\tau$         (long_name='Aggregate separation rate')
// Group 6: aggregate identity
    L           ${L}$           (long_name='Total labor')
    Le          ${L^{e}}$       (long_name='Labor in entry')
    Lc          ${L^{c}}$       (long_name='Labor in consumption')
    X           ${X}$           (long_name='Hiring cost')
    C           ${C}$           (long_name='Consumption')
    Yc          ${Y^{c}}$       (long_name='Retail output')
    Ygross      ${Y^{Gross}}$   (long_name='Gross output')
    Y           ${Y}$           (long_name='Net output')
    ls          ${\ell s}$      (long_name='Labor share')
    lp          ${\ell p}$      (long_name='Labor productivity')
// Group 7: exogenous process
    z           ${z}$           (long_name='Tech shocks')
    s           ${s}$           (long_name='Separation shocks')
    delta       ${\delta}$      (long_name='Product destruction shocks')
// Group 8: observables
    u_obs
    v_obs
    tau_obs
    f_obs
    delta_obs
    bf_obs
    lp_obs
// Group 9: log observables
    u_obs_lag
    v_obs_lag
    tau_obs_lag
    f_obs_lag
    delta_obs_lag
    bf_obs_lag
    lp_obs_lag
;

// Define exogenous shocks
varexo 
    e_z         ${e_z}$         (long_name='Technology shock')
    e_s         ${e_s}$         (long_name='Job separation shock')
    e_delta     ${e_\delta}$    (long_name='Firm destruction shock')
;


///////////////////////////////////////////////////////////////////////////
// Model parameters
///////////////////////////////////////////////////////////////////////////

parameters
    // Labor market parameters
    sigma       ${\sigma}$      (long_name='Inverse of intertemporal elasticity of substitution')
    b_ratio     ${b_{ratio}}$   (long_name='Replacement ratio: unemployment benefits to wages')
    xi_inv      ${\xi_{inv}}$   (long_name='Inverse elasticity of entrants to vacancy value')
    x_v         ${x_{v}}$       (long_name='Share of flow costs from vacancy creation')
    // Shock process parameters
    rho_z       ${\rho_z}$      (long_name='Persistence of tech shocks')
    rho_s       ${\rho_s}$      (long_name='Persistence of separation shocks')
    rho_delta   ${\rho_delta}$  (long_name='Persistence of product destruction shocks')
    sigma_z     ${\sigma_z}$    (long_name='Std. dev of tech shocks')
    sigma_s     ${\sigma_s}$    (long_name='std. dev of separation shocks')
    sigma_delta ${\sigma_delta}$(long_name='Std. dev of prod. dest. shocks')
    // Predetermined constant parameters
    eta         ${\eta}$        (long_name='Elasticity of matching function with respect to unemployment')
    F           ${F}$           (long_name='Parameter of vacancy creation')
    epsi        ${\varepsilon}$ (long_name='Elasticity of substitution')
    zeta        $\zeta$         (long_name='Love for variety')
    XY          $XY$            (long_name='Share of hiring cost in GDP')
    // Implied parameters from the steady state
    beta        ${\beta}$       (long_name='Discount factor')
    kappa       ${\kappa}$      (long_name='Cost of job posting')
    f_e         ${f_e}$         (long_name='Cost of business creation')
    A           ${A}$           (long_name='Match efficiency')
    x_m         ${x_{m}}$       (long_name='Cost of job creation')
    phi         ${\phi}$        (long_name='Bargaining power')
;


///////////////////////////////////////////////////////////////////////////
// Model equations
///////////////////////////////////////////////////////////////////////////

model;

// Abbreviations

    #mu = epsi/(epsi-1);
    #w_ss = STEADY_STATE(w);

// Group 1: search friction

    [name='Market tightness']
    theta = v/u(-1);
    
    [name='Job finding probability']
    f = A*theta^(1-eta);
    
    [name='Vacancy filling probability']
    q = A*theta^(-eta);

// Group 2: relative price

    [name='Relative price']
    p = N(-1)^zeta;

// Group 3: business formation

    [name='Lagrangian multiplier']
    lam = C^(-sigma);

    [name='Firm value']
    nu = beta*(1-delta)*lam(+1)/lam*(nu(+1)+d(+1));

    [name='Free entry']
    Ne*nu = wint*Le;
    
    [name='Entry technology']
    Ne*f_e = z*Le;

// Group 4: Vacancy creation

    [name='Job creation condition']
    (kappa+K/q) = beta*(1-delta)*lam(+1)/lam*(wint(+1)-w(+1)-K(+1)+(1-s(+1))*(kappa+K(+1)/q(+1)));

    [name='Expected discounted difference in vacancy value']
    K = Q - beta*(1-delta)*lam(+1)/lam*Q(+1);
    
    [name='Value of a vacancy']
    Q = x_m*(e/F)^xi_inv;
    
    [name='Marginal cost of production']
    wint = p*z/mu;
    
    [name='Wage equation']
    w = phi*(wint-K+f*(kappa+K/q)) + (1-phi)*b_ratio*w_ss;
    
// Group 5: laws of motion
    
    [name='LOM of firms']
    N = (1-delta)*(N(-1)+Ne);
    
    [name='LOM of unemployment']
    u = (1-(1-delta)*f)*u(-1) + tau*(1-u(-1));
    
    [name='LOM of vacancies']
    v = (1-delta(-1))*((1-q(-1))*v(-1)+s(-1)*(1-u(-1))) + e;

    [name='Aggregate separation rate']
    tau = 1-(1-delta)*(1-s);

// Group 6: aggregate identity
    
    [name='Aggregate labor']
    L = 1 - u(-1);
    
    [name='Composition of labor']
    L = Lc + Le;    
    
    [name='Retail output: resources']
    Yc = p*z*Lc;
    
    [name='Sunk vacancy posting costs']
    X = Q*e/(1+xi_inv) + kappa*v*q;
    
    [name='Retail output: expenditure']
    Yc = C + X;    

    [name='Gross output: expenditure']
    Ygross = Yc + Ne*nu;
    
    [name='Gross output: income']
    Ygross = wint*L + N(-1)*d;

    [name='GDP']
    Y = Ygross - X;

    [name='Labor share']
    ls = w*L/Y;

    [name='Labor productivity (data consistent)']
    lp = Y/(p*L);

// Group 7: exogenous process

    [name='Labor productivity process']
    log(z) = rho_z*log(z(-1)) + (1-rho_z)*log(steady_state(z)) - sigma_z*e_z;

    [name='Separation rate process']
    log(s) = rho_s*log(s(-1)) + (1-rho_s)*log(steady_state(s)) + sigma_s*e_s;

    [name='Product line destruction process']
    log(delta) = rho_delta*log(delta(-1)) + (1-rho_delta)*log(steady_state(delta)) + sigma_delta*e_delta;

// Group 8: observables

    [name='Unemployment observation']
    u_obs = log(1/3*(u+u(-1)+u(-2)));

    [name='Vacancy observation']
    v_obs = log(1/3*(v+v(-1)+v(-2)));

    [name='Job separation observation']
    tau_obs = log(1/3*(tau+tau(-1)+tau(-2)));

    [name='Job finding rate']
    f_obs = log(1/3*((1-delta)*f+(1-delta(-1))*f(-1)+(1-delta(-2))*f(-2)));

    [name='Business destruction observation']
    delta_obs = log(delta+delta(-1)+delta(-2));

    [name='Business formation observation']
    bf_obs = log(1/3*(Ne+Ne(-1)+Ne(-2)));

    [name='Labor productivity observation']
    lp_obs = log(1/3*(lp+lp(-1)+lp(-2)));

// Group 9: lag observables

    [name='Lag unemployment observation']
    u_obs_lag = u_obs(-3);

    [name='Lag vacancy observation']
    v_obs_lag = v_obs(-3);

    [name='Job separation observation']
    tau_obs_lag = tau_obs(-3);

    [name='Job finding rate']
    f_obs_lag = f_obs(-3);

    [name='Lag business destruction observation']
    delta_obs_lag = delta_obs(-3);

    [name='Lag business formation observation']
    bf_obs_lag = bf_obs(-3);

    [name='Lag labor productivity observation']
    lp_obs_lag = lp_obs(-3);

end;
    

///////////////////////////////////////////////////////////////////////////

resid;
steady;
check;

shocks;
var e_z = 1;
var e_s = 1;
var e_delta = 1;
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
sigma,      gamma_pdf,      1.0, 0.5;
b_ratio,    beta_pdf,       0.7, 0.1;
xi_inv,     gamma_pdf,      2.0, 1.0;
x_v,        beta_pdf,       0.5, 0.2;
rho_z,      beta_pdf,       0.8, 0.10;
rho_s,      beta_pdf,       0.6, 0.20;
rho_delta,  beta_pdf,       0.9, 0.05;
sigma_z,    inv_gamma_pdf,  0.01, 0.005;
sigma_s,    inv_gamma_pdf,  0.04, 0.02;
sigma_delta,inv_gamma_pdf,  0.04, 0.02;
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
           mode_compute=0,
           mode_file=model_mh_mode,
           mh_replic=10000,%910000,%1500000,%720000,           
           mh_nblocks=1,%11,
           mh_init_scale=0.5,
           mh_jscale=0.5,
           mh_drop=0.50);
    
// evaluate at MCMC joint mode

load model_mh_mode.mat;
for ii=1:1:length(xparam1)
    para_idx = find(strcmp(M_.param_names, char(parameter_names(ii))));
    // para_idx=loc(M_.param_names,char(parameter_names(ii)));
    M_.params(para_idx)=xparam1(ii);
end

// call steady state
model_steadystate([],[]);

// Stochastic simulation
stoch_simul(order=1,nofunctions,noprint)
    u_obs, v_obs, tau_obs, f_obs, delta_obs, bf_obs, lp_obs;
stack_model_moments;
Data_vs_Model = [Values', model_moments];
save results;
delete *.eps;