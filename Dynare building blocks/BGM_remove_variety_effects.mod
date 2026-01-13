var C N_e N d nu w L L_C L_E Y Y_r rho 
    log_Y log_Y_r log_L log_L_C log_L_E
    Z log_Z;

varexo e_z;

parameters beta r delta chi psi epsi f_e sigma rho_z;

//****************************************************************************
//Set parameter values
//**********

beta = 0.99;
r = (1-beta)/beta;
delta = 0.025;
chi = 0.9241;
psi = 4.0;
epsi = 3.8;
f_e = 1.0;
sigma = 0.9;
//Z = 1;

rho_z = 0.979;
//sigma_z = 0.0072;

//****************************************************************************
//enter the model equations (model-block)
//*********
model;
#mu = epsi/(epsi-1);
[name='Euler equation']
C^(-sigma) = beta*(1-delta)*C(+1)^(-sigma)*(nu(+1)+d(+1))/nu;

[name='Variety effects']
//rho = N^(1/(epsi-1)); 
rho = 1.0;

[name='Consumption sector'] 
C = Z*rho*L_C;

[name ='Labor in entry']
L_E = L - L_C;

[name = 'Wage from pricing']
w = rho*Z/mu;

[name= 'Labor supply']
L = (w*C^(-sigma)/chi)^(psi);

[name = 'Free entry']
nu = w*f_e/Z;

[name = 'Profits']
d = C/(N*epsi);

[name = 'Aggregate accounting']
nu*N_e + C = w*L+N*d;

[name = 'Output']
Y = C + nu*N_e;
[name = 'Data-consistent output']
Y_r = Y/rho;

[name='Firm law of motion']
N = (1-delta)*(N(-1) +N_e(-1));

% Exogenous processes 
log(Z) = rho_z*log(Z(-1)) + e_z;

log_Y = 100*log(Y);
log_Y_r = 100*log(Y_r);
log_L = 100*log(L);
log_L_C = 100*log(L_C);
log_L_E = 100*log(L_E);
log_Z = 100*log(Z);
end;

initval; // closed-form steady state inside the .mod 7 8
Z = 1;
L = 0.87;
N = (1-delta)*Z*L/(f_e*(r+delta)*(epsi-1)+delta);
//rho = N^(1/(epsi-1));
rho = 1.0;
w = rho*Z/(epsi/(epsi-1));
L_C = (epsi-1)*(r+delta)*L/((epsi-1)*(r+delta)+delta);
L_E = L - L_C;
C = Z*rho*L_C;
N_e = delta*N/(1-delta);
d = C/(epsi*N);
nu = w*f_e/Z;
Y = C+nu*N_e;
Y_r = (Y)/rho;

log_Y = 100*log(Y);
log_Y_r = 100*log(Y_r);
log_L = 100*log(L);
log_L_C = 100*log(L_C);
log_L_E = 100*log(L_E);
log_Z = 100*log(Z);
end;

shocks;
    //var e_z=0.0072^2;
    var e_z = 0.01^2;
end;

//check the starting values for the steady state
resid;
// compute steady state given the starting values
steady;
// check Blanchard Kahn conditions 
check;


stoch_simul(order=2, irf=40) log_Z log_Y log_Y_r log_L log_L_C log_L_E N rho;