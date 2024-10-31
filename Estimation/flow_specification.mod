

// define variables

var //x           ${x}$ (long_name='technology')
    w        ${w}$ (long_name = 'wage')
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
    lam       ${\lambda}$ (long_name='Marginal utility of consumption')
    Y         ${Y}$ (long_name='aggregate output')

    L         ${L}$ (long_name='total labor')
    L_e       ${L_e}$ (long_name='Labor in entry')
    L_c       ${L_c}$ (long_name='Labor in consumption')
    labor_share_obs (lone_name='Labor share')
    
    //s_agg
    theta_z     ${\theta_z}$ (long_name='technology shock')
    

    
    % Observation equations
    z_obs_m % monthly growth rate of labor productivity
    z_obs % quarterly growth rate of labor productivity
    u_obs % demeaned unemployment rate
    v_obs % demeaned vacancy rate
    //s_obs % demeaned aggregate separation rate

    u_imp
    v_imp
    theta_imp
    N_imp
    N_e_imp
    C_imp;
    

varexo e_z ${e_z}$ (long_name= 'technology shock')
    ;
    
parameters 
% A, F, and s are dependent parameters and are thus commented out
    w_ss     ${w_{ss}}$        (long_name='Steady-state wage')
    beta  ${\beta}$    (long_name='Discount factor')
    delta ${\delta}$   (long_name='Product destruction rate')
    //epsi  ${\varepsilon}$ (long_name='elasticity of substitution')
     mu_net  ${\mu-1}$ (long_name='net markup')
    sigma  ${\sigma}$ (long_name='inverse of intertemporal elasticity of substitution')
    //s ${S}$            (long_name='Worker separation rate')
    tau ${\tau}$       (long_name='Aggregate separation rate')

    fbar
    qbar ${\overline{q}}$ (long_name ='Mean vacancy filling rate')

    //alpha_L $\alpha_L$ (long_name = 'Worker bargaining power')      
    b_ratio $b$              (long_name = 'Replacement ratio: unemployment benefits to wages')
    //A $A$              (long_name='Matching function level parameter')
    eta_L $\eta_L$         (long_name='Elasticity of matching function with respect to unemployment')
    
    labor_share // replaces alpha_L
    N_ss ${N_{ss}}$ // replaces f_e

    rho_z    ${\rho_z}$               (long_name='persistence technology shock')
    ;

%----------------------------------------------------------------
% set parameter values 
%----------------------------------------------------------------
w_ss = 1.0;
beta = 0.99673; % monthly discount factor
xi_inv = 1; % 0 corresponds to infinitely elastic vacancy creation--standard DMP free entry
//r_ann = 0.04; % annual interest rate 

delta = 0.00514;
tau = 0.031;
//fbar = 0.41;
fbar = 0.41;
qbar = 0.8;

labor_share = 0.66;
N_ss = 1.0;
//mu_net = 0.3;
mu_net = 0.05;
sigma = 1.5;
b_ratio = 0.71;
eta_L = 0.6;


rho_z = 0.979;

%----------------------------------------------------------------
% enter model equations
%----------------------------------------------------------------

model;

% Dependent parameters
# rho=(1-beta)/beta;

% Correct job finding and vacancy finding probabilities
#f_ss = fbar/(1-delta);
#q_ss = qbar/(1-delta);
#theta_ss = f_ss/q_ss;
#u_ss = tau/(tau+(1-delta)*f_ss);
#L_ss = 1-u_ss;

#b = b_ratio*w_ss;

% Vacancies
#v_ss = theta_ss*u_ss;
%
#A = f_ss/(theta_ss^(1-eta_L));
% Separation rate
#s = (tau-delta)/(1-delta);

#epsi = (mu_net+1)/mu_net;
#mu = epsi/(epsi-1);
#p_ss = N_ss^(1/(epsi-1));
#N_e_ss = delta/(1-delta)*N_ss;

#recruiter_share= (delta+(rho+delta)*(epsi-1))/(delta+(rho+delta)*epsi);
#w_wR = labor_share/recruiter_share;
#w_R_ss = w_ss/(w_wR);
#z = (mu/p_ss)*w_R_ss;

#gam = q_ss*(1-delta)/(rho+tau)*(w_R_ss-w_ss);

// From wage equation find ϕ
#phi = (w_ss-b)/(w_R_ss+theta_ss*gam/(1-delta)-b);


#f_e = (mu-1)*z*(L_ss/N_ss)*(1-delta)/(delta*mu+rho);
    
[name= 'Job creation condition']
gam/q = beta*lam(+1)/lam*(1-delta)*(w_R(+1)-w(+1)+(1-s)*gam/q(+1));

[name ='Marginal revenue product']
w_R = p*z*exp(theta_z)/mu;

[name = 'Wage equation']
w = phi*(w_R +theta*gam/(1-delta)) + (1-phi)*b;

[name = 'Market tightness']
theta = v/u;

[name = 'Job finding probability']
f = A*theta^(1-eta_L);

[name = 'Vacancy filling probability']
q = A*theta^(-eta_L);

[name='Aggregate labor']
L = 1-u;

[name = 'Composition of labor']
L = L_c + L_e;

[name='Relative price']
p = N^(1/(epsi-1));

[name='Lagrangian multiplier']
lam = C^(-sigma);

[name='Firm value']
nu_f = beta*(1-delta)*lam(+1)/lam*(nu_f(+1)+d_f(+1));

[name ='Retail output: resources']
Y_c = p*z*exp(theta_z)*L_c;

[name ='Retail output: expenditure']
Y_c = C+gam*v;

[Name = 'Business entrants']
N_e = L_e*z*exp(theta_z)/f_e;

[Name = 'Firm value relative to price']
nu_f = p*f_e/mu;

[Name= 'Output: expenditure']
Y = Y_c + nu_f*N_e;

[Name = 'Output: income']
Y = w_R*L + N*d_f;

[name = 'LOM of unemployment']
u = (1-(1-delta)*f(-1))*u(-1) + tau*(1-u(-1)); % u is predetermined

[name = 'LOM of firms']
N = (1-delta)*(N(-1) + N_e(-1)); % N is predetermined

[name = 'labor share']
labor_share_obs = w*L/Y;


[name = 'Aggregate separation rate']
//s_agg = (1-(1-delta*exp(theta_delta))*(1-s));

% Observation variables: first differences (demeaned) -> link to data in first differences (p. 58 of Pfeifer's Observation Equations)
z_obs_m = theta_z - theta_z(-1); % Monthly growth rate of productivity (Delta \log x_t = theta_xt - theta_{x,t-1})

% Quarterly growth rate of labor productivity is sum of past 3 months
z_obs = z_obs_m + z_obs_m(-1) + z_obs_m(-2);

//u_obs = log(u/u(-1)); % Growth rate of unemployment
u_obs = u - steady_state(u);
//v_obs = log(v/v(-1)); % Growth rate of vacancies
v_obs = v - steady_state(v);
//s_obs = s_agg - tau;

% Variables for impulse responses (% deviation from steady state)
u_imp = 100*log(u);
v_imp = 100*log(v);
theta_imp = 100*log(theta);
N_imp = 100*log(N);
N_e_imp = 100*log(N_e);
C_imp   = 100*log(C);
//Equivalent of x_imp is theta_x

% Exogenous processes
[name ='Labor productivity process']
theta_z = rho_z*theta_z(-1) - e_z;

end;

steady_state_model;
    //Do Calibration
    //Steady state imposes fbar, qbar, r_ss, tau
    w = w_ss;
    N = N_ss;
    mu_ss = 1+mu_net;
    epsi_ss = (mu_net+1)/mu_net;
    b_ss = b_ratio*w_ss;

    p = N^(1/(epsi_ss-1));
    N_e = delta/(1-delta)*N;

    rho_ss = (1-beta)/beta;
    f = fbar/(1-delta);
    q = qbar/(1-delta);

    theta = f/q;
    u = tau/(tau+(1-delta)*f);
    L = 1 - u;
    v = theta*u;

    recruiter_share_ss= (delta+(rho_ss+delta)*(epsi_ss-1))/(delta+(rho_ss+delta)*epsi_ss);
    w_wR_ss = labor_share/recruiter_share_ss;
    w_R = w/(w_wR_ss);
    z_ss = (mu_ss/p)*w_R;

    f_e_ss = (mu_ss-1)*z_ss*(L/N)*(1-delta)/(delta*mu_ss+rho_ss);

    nu_f = p*f_e_ss/mu_ss;
    d_f = (rho_ss+delta)/(1-delta)*nu_f;

    L_e = N_e*f_e_ss/z_ss;
    L_c = L - L_e;
    Y_c = p*z_ss*L_c;
   
    gam_ss = q*(1-delta)/(rho_ss+tau)*(w_R-w);
    C = Y_c - gam_ss*v;
    lam = C^(-sigma);
    Y = Y_c + nu_f*N_e;

    labor_share_obs = labor_share;

    theta_z = 0;

    z_obs = 0;
    z_obs_m = 0;
    u_obs = 0;
    //s_obs = 0;
    v_obs = 0;

    z_imp = 0;
    u_imp = 100*(log(u));
    v_imp = 100*(log(v));
    theta_imp = 100*(log(theta));
    N_imp = 100*(log(N));
    N_e_imp = 100*(log(N_e));
    C_imp   = 100*(log(C));

end;

//set shock variances
shocks;
    var e_z=0.007^2;
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


% Stochastic simulation -> for conditional FEVD and IRF
stoch_simul (order=1, nofunctions, irf=80, periods=0)
//conditional_variance_decomposition=[1 4 8 40])
u_obs, v_obs, z_obs, theta_z, 
u_imp, v_imp, theta_imp, N_imp, N_e_imp, C_imp, labor_share_obs, L_c, L_e;

