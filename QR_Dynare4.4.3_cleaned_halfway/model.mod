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
//
///////////////////////////////////////////////////////////////////////////


// Model declaration and variables

var

    // Calvo wage variables
    h_1         ${h_1}$         (long_name = 'Calvo wage numerator (auxiliary)')
    h_2         ${h_2}$         (long_name = 'Calvo wage denominator (auxiliary)')
    wp          ${w^p}$         (long_name = 'adjusted real wage (log)')
    w           ${w}$           (long_name = 'average real wage (log)')
    
    // Shopping friction variables
    Q           ${Q}$           (long_name = 'market tightness (log)')
    Z           ${Z}$           (long_name = 'shopping disutility scaling (log)')
    CA          ${CA}$          (long_name = 'consumption aggregate (log)')
    
    // Shadow prices
    lambdax     ${\lambda_x}$   (long_name = 'shadow price of investment')
    lambdak     ${\lambda_k}$   (long_name = 'shadow price of capital')
    lambdae     ${\lambda_e}$   (long_name = 'shadow price of wealth')
    
    // Investment and capital
    Xk          ${X_k}$         (long_name = 'investment (log)')
    uk          ${u_k}$         (long_name = 'capacity utilization rate (log)')
    rk          ${r^k}$         (long_name = 'real rental rate of capital (log)')
    K           ${K}$           (long_name = 'capital stock (log)')
    
    // Matching and production
    Ical        ${\mathcal{I}}$ (long_name = 'occupancy rate (log)')
    Ycal        ${\mathcal{Y}}$ (long_name = 'production capacity (log)')
    L           ${L}$           (long_name = 'labor input (log)')
    
    // Prices and policy rates
    mc          ${mc}$          (long_name = 'real marginal cost (log)')
    tau         ${\tau}$        (long_name = 'price-cost markup (log)')
    Pi          ${\Pi}$         (long_name = 'gross inflation rate (log)')
    R           ${R}$           (long_name = 'gross nominal interest rate (log)')
    Ra          ${R^a}$         (long_name = 'pre-shock nominal rate (log)')
    
    // Technology trends
    mun         ${\mu_n}$       (long_name = 'neutral technology growth (log)')
    mux         ${\mu_x}$       (long_name = 'investment tech growth (log)')
    muy         ${\mu_y}$       (long_name = 'output growth (log)')
    muk         ${\mu_k}$       (long_name = 'capital growth (log)')
        
    // Abbreviations (levels)
    U           ${U}$           (long_name = 'utility')
    Gx          ${G_x}$         (long_name = 'investment growth')
    S           ${S}$           (long_name = 'investment adjustment cost')
    Sp          ${S\prime}$     (long_name = 'derivative of adjustment cost')
    a           ${a}$           (long_name = 'capital utilization cost')
    ap          ${a\prime}$     (long_name = 'derivative of utilization cost')
    chi         ${\chi}$        (long_name = 'price adjustment cost')
    Mcal        ${\mathcal{M}}$ (long_name = 'stochastic discount factor')
    
    // Abbreviations (logs)
    muu         ${\mu_u}$       (long_name = 'marginal utility growth rate')
    C           ${C}$           (long_name = 'consumption (log)')
    I           ${I}$           (long_name = 'investment (log)')
    Y           ${Y}$           (long_name = 'output (log)')
    lp          ${lp}$          (long_name = 'labor productivity (log)')
    ls          ${ls}$          (long_name = 'labor share (log)')
    u           ${u}$           (long_name = 'composite capital utilization (log)')

;

// Exogenous shocks
varexo
    epsR        ${\epsilon_R}$  (long_name = 'monetary policy shock (log)')
    epsn        ${\epsilon_n}$  (long_name = 'neutral technology shock (log)')
    epsx        ${\epsilon_x}$  (long_name = 'investment-specific tech shock (log)')
;


///////////////////////////////////////////////////////////////////////////

// Model parameters

parameters

    // Estimated parameters
    sigma_a     ${\sigma_a}$     (long_name = 'curvature of utilization cost')
    ssigma      ${\sigma}$       (long_name = 'inverse IES')
    xi          ${\xi}$          (long_name = 'inverse Frisch elasticity')
    h           ${h}$            (long_name = 'habit persistence')
    Spp         ${S\prime\prime}$(long_name = 'curvature of investment adjustment cost')
    kappawt     ${\kappa_w}$     (long_name = 'curvature of price adjustment cost')
    theta_w     ${\theta_w}$     (long_name = 'Calvo wage rigidity')
    phi_pi      ${\phi_\pi}$     (long_name = 'Taylor inflation coefficient')
    phi_y       ${\phi_y}$       (long_name = 'Taylor output coefficient')
    rho_R       ${\rho_R}$       (long_name = 'interest rate smoothing')
    sigma_R     ${\sigma_R}$     (long_name = 'std dev of FFR shock')
    sigma_n     ${\sigma_n}$     (long_name = 'std dev of neutral tech shock')
    sigma_x     ${\sigma_x}$     (long_name = 'std dev of investment tech shock')
    rho_x       ${\rho_x}$       (long_name = 'persistent of investment tech shock')
    
    // Calibrated and implied parameters
    rho         ${\rho}$         (long_name = 'benefit of search')
    nu          ${\nu}$          (long_name = 'cost of search')
    varphi      ${\varphi}$      (long_name = 'bottleneck of search')
    ggamma      ${\gamma}$       (long_name = 'complementarity of search')
    Psi         ${\Psi}$         (long_name = 'macro elasticity of occupancy')
    bbeta       ${\beta}$        (long_name = 'discount factor')
    delta_k     ${\delta_k}$     (long_name = 'capital depreciation')
    aalpha      ${\alpha}$       (long_name = 'capital share')
    vartheta    ${\vartheta}$    (long_name = 'fixed cost of production')
    zeta        ${\zeta}$        (long_name = 'shopping disutility level')
    eta         ${\eta}$         (long_name = 'working disutility level')
    sigma_b     ${\sigma_b}$     (long_name = 'utilization cost level')
    
    // Fixed parameters
    omega       ${\omega}$       (long_name = 'technology diffusion speed')
    rho_w       ${\rho_w}$       (long_name = 'gross wage markup')
    varrho      ${\varrho}$      (long_name = 'fraction of searched varieties')
    varpi       ${\varpi}$       (long_name = 'fraction of reinstalled inputs')
    iota        ${\iota}$        (long_name = 'CRRA coefficient')

;


///////////////////////////////////////////////////////////////////////////

model;
    
// Steady-state abbreviations

    #Pi_SS      = STEADY_STATE(Pi);
    #Ical_SS    = STEADY_STATE(Ical);
    #Ycal_SS    = STEADY_STATE(Ycal);
    #tau_SS     = STEADY_STATE(tau);
    #Y_SS       = STEADY_STATE(Y);
    #R_SS       = STEADY_STATE(R);
    #mun_SS     = STEADY_STATE(mun);
    #mux_SS     = STEADY_STATE(mux);
    #muk_SS     = STEADY_STATE(muk);
    #muy_SS     = STEADY_STATE(muy);
    #ls_SS      = STEADY_STATE(ls);


///////////////////////////////////////////////////////////////////////////

// Wage setting equations (Calvo-style)

[name = 'eq01 Calvo wage']
-h_1 + eta*exp((1+xi)*(rho_w/(rho_w-1)*w+L)) 
        + bbeta*theta_w*exp((rho_w/(rho_w-1)*(1+xi))*(Pi(+1)+muy(+1)))*h_1(+1);

[name = 'eq02 Calvo wage']
-h_2 + lambdae*exp(rho_w/(rho_w-1)*w+L) 
        + bbeta*theta_w*exp(1/(rho_w-1)*(Pi(+1)+muy(+1))+muu(+1))*h_2(+1);

[name = 'eq03 Calvo wage']
-exp((1+rho_w/(rho_w-1)*xi)*wp) + rho_w*h_1/h_2;

[name = 'eq04 Calvo wage']
-exp(1/(1-rho_w)*w) 
        + (1-theta_w)*exp(1/(1-rho_w)*wp) 
        + theta_w*exp(1/(1-rho_w)*(w(-1)-Pi-muy));

// eq05. to be revised
//- (rho-1)*exp(rho*Ical+Ycal) + zeta*exp(Z+(1+nu)*Q);
//- exp(Z)*(exp(Q)*( 1-(1-varrho)*exp(Ical(-1))) )^(1+nu)/(1-(1-varrho)*exp(Ical(-1)-Ical)) + bbeta*(1-varrho)*exp(lambdae(+1)-lambdae)*exp((Ical(+1)-Ical)*(-rho))*exp(Z(+1))*(exp(Q(+1))*( 1-(1-varrho)*exp(Ical)) )^(1+nu)/(1-(1-varrho)*exp(Ical-Ical(+1))) + (rho-1)/zeta*exp(rho*Ical+Ycal);
- exp(Z)*(exp(Q)*( 1-(1-varrho)*exp(Ical_SS)) )^(1+nu)/(1-(1-varrho)*exp(Ical_SS-Ical)) + (rho-1)/zeta*exp(rho*Ical+Ycal);

// eq06.
- lambdae + exp((rho-1)*Ical)
    *((U-h*U(-1)/exp(muy))^(-ssigma)-bbeta*h*(U(+1)*exp(muy(+1))-h*U)^(-ssigma));
// eq07.
- lambdae + exp((rho-1)*Ical)*lambdax;
// eq08.
// - lambdax + lambdak*(1-S-Sp*Gx) + bbeta*lambdak(+1)*exp(muu(+1)-muk(+1))*Sp(+1)*Gx(+1)^2;
- lambdax + lambdak*(1-S-Sp*Gx) + bbeta*lambdak(+1)*exp(muu(+1)-muk(+1))*Sp(+1)*Gx(+1)^2*exp(mux(+1));
// eq09.
- lambdak + bbeta*(lambdae(+1)*exp(rk(+1)+uk(+1))-lambdax(+1)*a(+1)+lambdak(+1)*(1-delta_k))*exp(muu(+1)-muk(+1));
// eq10.
- lambdae*exp(rk) + lambdax*ap;
// eq11.
- lambdae + bbeta*lambdae(+1)*exp(muu(+1)-muy(+1)+R(+1)-Pi(+1));
// eq12.
- exp(rho*Ical+Ycal) + exp(CA) + exp(Xk) + a*exp(K(-1)-muk);
// eq13.
- exp(K) + (1-delta_k)*exp(K(-1)-muk) + (1-S)*exp(Xk);

// eq14. to be revised
//- exp(Ical-Ical_SS) + (1-varphi+varphi*exp(-ggamma*Q))^(-1/ggamma);
//- (exp(Ical)-(1-varrho)*exp(Ical(-1)))/(varrho*exp(Ical_SS)) * (1-(1-varrho)*exp(Ical_SS))/(1-(1-varrho)*exp(Ical(-1))) + (1-varphi+varphi*exp(-ggamma*Q))^(-1/ggamma);
- (exp(Ical)-(1-varrho)*exp(Ical_SS))/(varrho*exp(Ical_SS)) + (1-varphi+varphi*exp(-ggamma*Q))^(-1/ggamma);

// eq15.
- (1+chi*exp(Ical))*exp(Ycal) + exp(aalpha*(uk+K(-1)-muk)+(1-aalpha)*L) - exp(Ycal_SS)*vartheta;
//((1+tau_SS)*exp(ls_SS)/(1-aalpha)-1);
// eq16.
- aalpha*exp(Ra+w+L) + (1-aalpha)*exp(rk+uk+K(-1)-muk);
// eq17
- exp(mc) + (exp(rk)/aalpha)^aalpha*(exp(Ra+w)/(1-aalpha))^(1-aalpha);

// eq18. to be revised
- (exp(Pi)-exp(Pi_SS))*exp(Pi)
    + kappawt*(1+tau_SS)/rho
        *( rho*exp(mc-Ical)
            - (1-varrho)*exp(Ical_SS-Ical)
                - (1 - (1-varrho)*exp(Ical_SS-Ical))
                    *(1+varphi/(1-varphi)*exp(-ggamma*Q))
                        + chi*(rho*exp(mc)-1) )
                            + 1.00*bbeta*Mcal(+1)*(exp(Pi(+1))-exp(Pi_SS))*exp(Pi(+1));
//- (exp(Pi)-exp(Pi_SS))*exp(Pi)
//    + kappawt*(1+tau_SS)/rho
//        *( rho*exp(mc-Ical)
//            - 1/(1- (varphi*exp(-ggamma*Q))/(1-varphi+varphi*exp(-ggamma*Q))*(1-(1-varrho)*exp(Ical(-1)-Ical)) )
//                + chi*(rho*exp(mc)-1) )
//                    + bbeta*Mcal(+1)*(exp(Pi(+1))-exp(Pi_SS))*exp(Pi(+1));
//- (exp(Pi)-exp(Pi_SS))*exp(Pi)
//    + kappawt*(1+tau_SS)/rho
//        *( rho*exp(mc-Ical)-exp(Ical-Ical_SS)^(-ggamma)/(1-varphi)
//            + chi*(rho*exp(mc)-1) )
//                + bbeta*Mcal(+1)*(exp(Pi(+1))-exp(Pi_SS))*exp(Pi(+1));
// - (exp(Pi)-exp(Pi_SS))*exp(Pi)
//     + kappawt*(1+tau_SS)/rho
//         *( rho*exp(mc-Ical)-exp(Ical)^(-ggamma)/(1-varphi)
//             + chi*(rho*exp(mc)-1) )
//                 + bbeta*Mcal(+1)*(exp(Pi(+1))-exp(Pi_SS))*exp(Pi(+1));
// - (exp(Pi)-exp(Pi_SS))*exp(Pi)
//     + kappawt*(1-varphi)*exp(Ical_SS)^(ggamma-1)
//         *(rho*exp(mc)-(chi+(1-chi)*exp(Ical)^(1-ggamma)/(1-varphi)))
//             + bbeta*Mcal(+1)*(exp(Pi(+1))-exp(Pi_SS))*exp(Pi(+1));
// - (exp(Pi)-exp(Pi_SS))*exp(Pi)
//     + kappawt*exp(Ical_SS)^(ggamma_PC-1)
//         *(rho*exp(mc)-(chi+(1-chi)*exp(Ical)^(1-ggamma_PC)))
//             + bbeta*Mcal(+1)*(exp(Pi(+1))-exp(Pi_SS))*exp(Pi(+1));

// eq19.
- (Ra-R_SS) + rho_R*(R-R_SS)
    + (1-rho_R)*(phi_pi*(Pi-Pi_SS)+phi_y*(Y-Y_SS));
//- (Ra-R_SS) + rho_R*(R-R_SS)
//    + (1-rho_R)*(phi_pi*(Pi-Pi_SS-0.25*(mux-mux_SS))+phi_y*(Y-Y_SS));
// eq20.
- R + Ra(-1) - sigma_R*epsR/400;

///////////////////////////////////////////////////////////////////////////

// Exogenous processes

[name = 'eq21 technology']
- (mun-mun_SS) + sigma_n*epsn/100;

[name = 'eq22 technology']
- (mux-mux_SS) + rho_x*(mux(-1)-mux_SS) + sigma_x*epsx/100;

[name = 'eq23 technology']
- muy + aalpha/(1-aalpha)*mux + mun;

[name = 'eq24 technology']
- muk + mux + muy;

[name = 'eq25 technology']
- Z + (1-omega)*(Z(-1)-muy);


// eq26.
- tau + exp(Ical-mc) - 1;


///////////////////////////////////////////////////////////////////////////

// abbreviations in levels

[name = 'a01 utility']
-U + exp(CA) - zeta*exp(Z)*(exp(Q)*(1-(1-varrho)*exp(Ical_SS)))^(1+nu)/(1+nu);

[name = 'a01 utility']
//- Gx + exp(Xk-Xk(-1)+muk);
- Gx + exp(Xk-Xk(-1)+muy);
// a03.
//- S + (exp(sqrt(Spp)*(Gx-exp(muk_SS)))+exp(-sqrt(Spp)*(Gx-exp(muk_SS))))/2 - 1;
- S + (exp(sqrt(Spp)*(Gx-exp(muy_SS)))+exp(-sqrt(Spp)*(Gx-exp(muy_SS))))/2 - 1;
// a04.
//- Sp + sqrt(Spp)*(exp(sqrt(Spp)*(Gx-exp(muk_SS)))-exp(-sqrt(Spp)*(Gx-exp(muk_SS))))/2;
- Sp + sqrt(Spp)*(exp(sqrt(Spp)*(Gx-exp(muy_SS)))-exp(-sqrt(Spp)*(Gx-exp(muy_SS))))/2;
// a05.
- a + sigma_a*sigma_b/2*exp(uk)^2 + sigma_b*(1-sigma_a)*exp(uk) + sigma_b*(sigma_a/2-1);
// a06.
- ap + sigma_a*sigma_b*exp(uk) + sigma_b*(1-sigma_a);
// a07.
- chi + rho/(rho-1)/(1+tau_SS/kappawt)/2*(exp(Pi-Pi_SS)-1)^2;
// - chi + (exp(Ical_SS)^(1-ggamma)/(2*(rho-1)*kappawt))*(exp(Pi-Pi_SS)-1)^2;
// a08.
- Mcal + lambdae/lambdae(-1)*exp(Ical-Ical(-1)+Ycal-Ycal(-1)+muu);


// abbreviations in logs

// a09.
- muu + (1-ssigma)*muy;
// a10.
- C - (rho-1)*Ical + CA;
// a11.
- I - (rho-1)*Ical + Xk;
// - exp(I) + exp(-(rho-1)*Ical+Xk) + a*exp(K(-1)-muk));
// a12.
- exp(Y) + exp(C) + exp(I);
// a13.
- lp + Y - L;
// a14.
- ls + w + L - Y;
// a15.
// - u + Ical + uk;
- u + Ical + uk*(((1+tau_SS)*ls_SS/(1-aalpha)-1)*Ical_SS+1)*aalpha;

end;


///////////////////////////////////////////////////////////////////////////

resid;
steady;
check;

shocks;
var epsR = 1;
var epsn = 1;
var epsx = 1;
end;


///////////////////////////////////////////////////////////////////////////

// Observed variables (irrelevant for impulse response matching

varobs Y;   // check out data.m. We just pass-on random data to prevent Dynare
            // from getting stuck here. VAR impulse response data is loaded and 
            // handled below.

// add folder with modified Dynare functions and utilities
addpath('DynareUtilites');


///////////////////////////////////////////////////////////////////////////

// Priors for estimation

estimated_params;
// vartheta,   gamma_pdf,  0.2,0.02;
sigma_a,    gamma_pdf,  0.5,0.3;
ssigma,     gamma_pdf,  1.0,0.15;
xi,         gamma_pdf,  1.0,0.05;
theta_w,    beta_pdf,   0.75,0.1;
kappawt,    gamma_pdf,  0.1,0.01;
h,          beta_pdf,   0.5,0.15;
Spp,        gamma_pdf,  8,2;
phi_pi,     gamma_pdf,  2.0,0.20;
phi_y,      gamma_pdf,  0.1,0.05;
rho_R,      beta_pdf,   0.7,0.15;
sigma_R,    gamma_pdf,  0.65,0.05;
sigma_n,    gamma_pdf,  0.1,0.05;
sigma_x,    gamma_pdf,  0.1,0.05;
rho_x,      beta_pdf,   0.75,0.1;
end;


///////////////////////////////////////////////////////////////////////////

// Begin of estimation

do_imp_resp_matching=1;                  %if =1, model estimated based on VAR impulse responses
if do_imp_resp_matching==1,
   addpath('DynareImpRespMatching');     %path with adjusted DYNARE functions dsge_likelihood.m and dynare_estimation_1.m  
   get_VAR_resp;                         %get VAR impulse responses: load VAR impulse responses and precision matrix
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
           %prior_trunc=0,
           prior_trunc=0,
           mode_compute=6,
           %mode_check,
           mode_file=model_mh_mode,
           %load_mh_file,
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

// plot VAR vs. Model resp.
stoch_simul(order=1,irf=16,nomoments,nocorr,nograph,noprint);
stack_model_responses;
plot_var_model_irf; 
save results;
delete *.eps;