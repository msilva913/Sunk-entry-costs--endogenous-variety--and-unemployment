%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function sets parameters and steady state for the nested model

% created by Zhesheng Qiu at 10:47 2019-01-15, Hong Kong
% updated by Zhesheng Qiu at 17:21 2019-01-31, Hong Kong
% updated by Zhesheng Qiu at 23:12 2019-01-31, Hong Kong
% updated by Zhesheng Qiu at 22:21 2019-02-03, Yangxin
% updated by Zhesheng Qiu at 14:59 2019-02-10, Hong Kong

% updated by Zhesheng Qiu at 12:40 2021-01-12, Hong Kong

% updated by Zhesheng Qiu at 18:39 2021-03-13, Hong Kong

% updated by Zhesheng Qiu at 23:40 2021-03-16, Hong Kong

% updated by Zhesheng Qiu at 16:40 2021-04-07, Hong Kong
% updated by Zhesheng Qiu at 13:50 2021-04-14, Hong Kong

% updated by Zhesheng Qiu at 12:18 2021-09-23, Hong Kong

% updated by Zhesheng Qiu at 17:41 2021-09-26, Hong Kong

% updated by Zhesheng Qiu at 17:01 2021-09-27, Hong Kong

% updated by Zhesheng Qiu at 18:26 2021-10-03, Hong Kong

% updated by Zhesheng Qiu at 19:02 2024-04-27, Hong Kong

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ys,check] = model_steadystate(ys,exo)
 
global M_


%% Read parameters

% if isempty(M_) == 1
%     M_.params = NaN;
% end

if isnan(M_.params)
    [ ...
    omega, rho_w, varrho, varpi, iota ...
    Pi_SS, muy_SS, muk_SS, mux_SS, ...
    ds_SS, xs_SS, ls_SS, R_SS, uk_SS, L_SS, Ical_SS, ...
    sigma_a, ssigma, xi, h, Spp, kappawt, theta_w, ...
        phi_pi, phi_y, rho_R, sigma_R, sigma_n, sigma_x, rho_x ...
    ] ...
    = params;
else
    [ ...
    omega, rho_w, varrho, varpi, iota ...
    Pi_SS, muy_SS, muk_SS, mux_SS, ...
    ds_SS, xs_SS, ls_SS, R_SS, uk_SS, L_SS, Ical_SS, ...
    ] ...
    = params;
    
    sigma_a     = M_.params(1);    
    ssigma      = M_.params(2);
    xi          = M_.params(3);    
    h           = M_.params(4);
    Spp         = M_.params(5);
    kappawt     = M_.params(6);    
    theta_w     = M_.params(7);    
    phi_pi      = M_.params(8);
    phi_y       = M_.params(9);
    rho_R       = M_.params(10);
    sigma_R     = M_.params(11);
    sigma_n     = M_.params(12);
    sigma_x     = M_.params(13);
    rho_x       = M_.params(14);
end


%% Solve steady state

Ra_SS       = R_SS;
muu_SS      = muy_SS^(1-ssigma);
bbeta       = Pi_SS/Ra_SS*muy_SS/muu_SS;

delta_k     = (muk_SS-1)/(xs_SS/ds_SS-1);
aalpha      = 1/(ls_SS/ds_SS*delta_k/(muk_SS/muy_SS*R_SS/Pi_SS-(1-delta_k))+1);
% Gamma       = (1-aalpha)/ls_SS*(1+vartheta);
Gamma       = 1.34;
vartheta    = Gamma*ls_SS/(1-aalpha)-1;

beta_miche  = (1/1.58+1.92/4.72)/2;%1/1.58;%1.92/4.72;
beta_stick  = (0.642+0.316)/2;
% beta_miche  = 0.010;
% beta_stick  = 0.001;
% beta_petro  = 107.95/100;
% % beta_stick  = 0.5;
% % beta_petro  = 2;
% 
% rho         = Gamma/(1-beta_stick);
% 
% a           = 1+(2-rho)*beta_stick/(1-beta_stick);
% b           = - (1-1/beta_petro+(2*(2-rho)-1/beta_miche)*beta_stick/(1-beta_stick));
% c           = - (1/beta_miche+rho-2)*beta_stick/(1-beta_stick);
% 
% varphi      = (-b+sqrt(b^2-4*a*c))/(2*a);
% varrho      = beta_stick/(1-beta_stick)*(1-varphi)/varphi;
% nu          = varrho*(1/beta_miche + rho - 2);
% 
% Psi         = varrho*varphi*beta_petro;
% 
% ggamma      = (gGamma-1)*varrho*varphi/beta_stick + 1;

varrho      = 1.00;
% varrho      = 0.25;
varphi      = beta_stick/(varrho+(1-varrho)*beta_stick);
rho         = Gamma/(1-beta_stick);
nu          = varrho*(1/beta_miche + rho - 2);
% nu          = 0;

% rho    = 2.54;
% varrho = 0.25;
% nu     = 0.74;
% varphi = 0.78;

Psi         = 1/((1+nu-varphi)/varphi/varrho+2-rho);
ggamma      = 0.509;

Q_SS        = 1;
mc_SS       = Ical_SS/Gamma;
sigma_b     = bbeta^(-1)*(muk_SS/muu_SS) - (1-delta_k);
rk_SS       = Ical_SS^(1-rho)*sigma_b;

w_SS        = (1-aalpha)*Ra_SS^(-1)*mc_SS^(1/(1-aalpha))*(rk_SS/aalpha)^(-aalpha/(1-aalpha));
K_SS        = (aalpha/(1-aalpha))*(Ra_SS*w_SS/rk_SS)*L_SS*muk_SS;
Ycal_SS     = 1/(1+vartheta)*(K_SS/muk_SS)^aalpha*L_SS^(1-aalpha);
Xk_SS       = (1-(1-delta_k)/muk_SS)*K_SS;
CA_SS       = Ical_SS^rho*Ycal_SS - Xk_SS;
Z_SS        = muy_SS^(-(1-omega)/omega);
%zeta        = (rho-1)*Ical_SS^rho*Ycal_SS/(Z_SS*Q_SS^(1+nu));
%zeta        = (rho-1)*Ical_SS^rho*Ycal_SS/(Z_SS*(Q_SS*(1-(1-varrho)*Ical_SS))^(1+nu)) * varrho/(1-bbeta*(1-varrho));
zeta        = (rho-1)*Ical_SS^rho*Ycal_SS/(Z_SS*(Q_SS*(1-(1-varrho)*Ical_SS))^(1+nu)) * varrho;

%U_SS        = CA_SS - zeta*Z_SS*(Q_SS^(1+nu))/(1+nu);
U_SS        = CA_SS - zeta*Z_SS*(Q_SS*(1-(1-varrho)*Ical_SS))^(1+nu)/(1+nu);

lambdae_SS  = Ical_SS^(rho-1)*(1-bbeta*h*muy_SS^(-ssigma))*((1-h/muy_SS)*U_SS)^(-ssigma);

wp_SS       = ((1-theta_w)/(1-theta_w*(Pi_SS*muy_SS)^(1/(rho_w-1))))^(rho_w-1)*w_SS;
h_2SS       = lambdae_SS*w_SS^(rho_w/(rho_w-1))*L_SS ...
            / (1-bbeta*theta_w*(Pi_SS*muy_SS)^(1/(rho_w-1))*muu_SS);
h_1SS       = (1/rho_w)*(wp_SS)^(1+rho_w/(rho_w-1)*xi)*h_2SS;
eta         = (1-bbeta*theta_w*(Pi_SS*muy_SS)^(rho_w/(rho_w-1)*(1+xi)))*h_1SS ...
            / (w_SS^(rho_w/(rho_w-1))*L_SS)^(1+xi);

lambdax_SS  = Ical_SS^(1-rho)*lambdae_SS;
lambdak_SS  = lambdax_SS;

mun_SS      = muy_SS/mux_SS^(aalpha/(1-aalpha));

Gx_SS       = muy_SS;
% Gx_SS       = muk_SS;
S_SS        = 0;
Sp_SS       = 0;
a_SS        = 0;
ap_SS       = sigma_b;
chi_SS      = 0;
Mcal_SS     = muu_SS;

C_SS        = Ical_SS^(1-rho)*CA_SS;
I_SS        = Ical_SS^(1-rho)*Xk_SS;
Y_SS        = Ical_SS*Ycal_SS;
lp_SS       = Y_SS/L_SS;
ls_SS       = w_SS*L_SS/Y_SS;

u_SS        = Ical_SS*uk_SS;


%% Put steady states into ys vector

ys          = NaN(41,1);

ys(1)       = h_1SS;
ys(2)       = h_2SS;
ys(3)       = log(wp_SS);
ys(4)       = log(w_SS);

ys(5)       = log(Q_SS);
ys(6)       = log(Z_SS);
ys(7)       = log(CA_SS);

ys(8)       = lambdax_SS;
ys(9)       = lambdak_SS;
ys(10)      = lambdae_SS;

ys(11)      = log(Xk_SS);
ys(12)      = log(uk_SS);
ys(13)      = log(rk_SS);
ys(14)      = log(K_SS);

ys(15)      = log(Ical_SS);
ys(16)      = log(Ycal_SS);
ys(17)      = log(L_SS);

ys(18)      = log(mc_SS);
ys(19)      = Gamma - 1;
ys(20)      = log(Pi_SS);
ys(21)      = log(R_SS);
ys(22)      = log(Ra_SS);

ys(23)      = log(mun_SS);
ys(24)      = log(mux_SS);
ys(25)      = log(muy_SS);
ys(26)      = log(muk_SS);

ys(27)      = U_SS;
ys(28)      = Gx_SS;
ys(29)      = S_SS;
ys(30)      = Sp_SS;
ys(31)      = a_SS;
ys(32)      = ap_SS;
ys(33)      = chi_SS;
ys(34)      = Mcal_SS;

ys(35)      = log(muu_SS);
ys(36)      = log(C_SS);
ys(37)      = log(I_SS);
ys(38)      = log(Y_SS);
ys(39)      = log(lp_SS);
ys(40)      = log(ls_SS);
ys(41)      = log(u_SS);

% steady state is transformed into log form if the original variable is
% also in log.


%% Save parameters

if isnan(M_.params)
    
    M_.params       = nan(29,1);
       
    M_.params(1)    = sigma_a;
    M_.params(2)    = ssigma;
    M_.params(3)    = xi;
    M_.params(4)    = h;
    M_.params(5)    = Spp;
    M_.params(6)    = kappawt;
    M_.params(7)    = theta_w;
    M_.params(8)    = phi_pi;
    M_.params(9)    = phi_y;
    M_.params(10)   = rho_R;
    M_.params(11)   = sigma_R;
    M_.params(12)   = sigma_n;
    M_.params(13)   = sigma_x;
    M_.params(14)   = rho_x;
        
end

M_.params(15)   = rho;
M_.params(16)   = nu;
M_.params(17)   = varphi;
M_.params(18)   = ggamma;
M_.params(19)   = Psi;
M_.params(20)   = bbeta;
M_.params(21)   = delta_k;
M_.params(22)   = aalpha;
M_.params(23)   = vartheta;
M_.params(24)   = zeta;
M_.params(25)   = eta;
M_.params(26)   = sigma_b;

M_.params(27)   = omega;
M_.params(28)   = rho_w;
M_.params(29)   = varrho;
M_.params(30)   = varpi;
M_.params(31)   = iota;

check = 0;