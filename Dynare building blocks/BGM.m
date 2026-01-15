%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'BGM';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('BGM.log');
M_.exo_names = 'e_z';
M_.exo_names_tex = 'e\_z';
M_.exo_names_long = 'e_z';
M_.endo_names = 'C';
M_.endo_names_tex = 'C';
M_.endo_names_long = 'C';
M_.endo_names = char(M_.endo_names, 'N_e');
M_.endo_names_tex = char(M_.endo_names_tex, 'N\_e');
M_.endo_names_long = char(M_.endo_names_long, 'N_e');
M_.endo_names = char(M_.endo_names, 'N');
M_.endo_names_tex = char(M_.endo_names_tex, 'N');
M_.endo_names_long = char(M_.endo_names_long, 'N');
M_.endo_names = char(M_.endo_names, 'd');
M_.endo_names_tex = char(M_.endo_names_tex, 'd');
M_.endo_names_long = char(M_.endo_names_long, 'd');
M_.endo_names = char(M_.endo_names, 'nu');
M_.endo_names_tex = char(M_.endo_names_tex, 'nu');
M_.endo_names_long = char(M_.endo_names_long, 'nu');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names_long = char(M_.endo_names_long, 'w');
M_.endo_names = char(M_.endo_names, 'L');
M_.endo_names_tex = char(M_.endo_names_tex, 'L');
M_.endo_names_long = char(M_.endo_names_long, 'L');
M_.endo_names = char(M_.endo_names, 'L_C');
M_.endo_names_tex = char(M_.endo_names_tex, 'L\_C');
M_.endo_names_long = char(M_.endo_names_long, 'L_C');
M_.endo_names = char(M_.endo_names, 'L_E');
M_.endo_names_tex = char(M_.endo_names_tex, 'L\_E');
M_.endo_names_long = char(M_.endo_names_long, 'L_E');
M_.endo_names = char(M_.endo_names, 'Y');
M_.endo_names_tex = char(M_.endo_names_tex, 'Y');
M_.endo_names_long = char(M_.endo_names_long, 'Y');
M_.endo_names = char(M_.endo_names, 'Y_r');
M_.endo_names_tex = char(M_.endo_names_tex, 'Y\_r');
M_.endo_names_long = char(M_.endo_names_long, 'Y_r');
M_.endo_names = char(M_.endo_names, 'rho');
M_.endo_names_tex = char(M_.endo_names_tex, 'rho');
M_.endo_names_long = char(M_.endo_names_long, 'rho');
M_.endo_names = char(M_.endo_names, 'log_Y');
M_.endo_names_tex = char(M_.endo_names_tex, 'log\_Y');
M_.endo_names_long = char(M_.endo_names_long, 'log_Y');
M_.endo_names = char(M_.endo_names, 'log_Y_r');
M_.endo_names_tex = char(M_.endo_names_tex, 'log\_Y\_r');
M_.endo_names_long = char(M_.endo_names_long, 'log_Y_r');
M_.endo_names = char(M_.endo_names, 'log_L');
M_.endo_names_tex = char(M_.endo_names_tex, 'log\_L');
M_.endo_names_long = char(M_.endo_names_long, 'log_L');
M_.endo_names = char(M_.endo_names, 'log_L_C');
M_.endo_names_tex = char(M_.endo_names_tex, 'log\_L\_C');
M_.endo_names_long = char(M_.endo_names_long, 'log_L_C');
M_.endo_names = char(M_.endo_names, 'log_L_E');
M_.endo_names_tex = char(M_.endo_names_tex, 'log\_L\_E');
M_.endo_names_long = char(M_.endo_names_long, 'log_L_E');
M_.endo_names = char(M_.endo_names, 'Z');
M_.endo_names_tex = char(M_.endo_names_tex, 'Z');
M_.endo_names_long = char(M_.endo_names_long, 'Z');
M_.endo_names = char(M_.endo_names, 'log_Z');
M_.endo_names_tex = char(M_.endo_names_tex, 'log\_Z');
M_.endo_names_long = char(M_.endo_names_long, 'log_Z');
M_.param_names = 'beta';
M_.param_names_tex = 'beta';
M_.param_names_long = 'beta';
M_.param_names = char(M_.param_names, 'r');
M_.param_names_tex = char(M_.param_names_tex, 'r');
M_.param_names_long = char(M_.param_names_long, 'r');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names_long = char(M_.param_names_long, 'delta');
M_.param_names = char(M_.param_names, 'chi');
M_.param_names_tex = char(M_.param_names_tex, 'chi');
M_.param_names_long = char(M_.param_names_long, 'chi');
M_.param_names = char(M_.param_names, 'psi');
M_.param_names_tex = char(M_.param_names_tex, 'psi');
M_.param_names_long = char(M_.param_names_long, 'psi');
M_.param_names = char(M_.param_names, 'epsi');
M_.param_names_tex = char(M_.param_names_tex, 'epsi');
M_.param_names_long = char(M_.param_names_long, 'epsi');
M_.param_names = char(M_.param_names, 'f_e');
M_.param_names_tex = char(M_.param_names_tex, 'f\_e');
M_.param_names_long = char(M_.param_names_long, 'f_e');
M_.param_names = char(M_.param_names, 'sigma');
M_.param_names_tex = char(M_.param_names_tex, 'sigma');
M_.param_names_long = char(M_.param_names_long, 'sigma');
M_.param_names = char(M_.param_names, 'rho_z');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_z');
M_.param_names_long = char(M_.param_names_long, 'rho_z');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 19;
M_.param_nbr = 9;
M_.orig_endo_nbr = 19;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('BGM_static');
erase_compiled_function('BGM_dynamic');
M_.lead_lag_incidence = [
 0 4 23;
 1 5 0;
 2 6 0;
 0 7 24;
 0 8 25;
 0 9 0;
 0 10 0;
 0 11 0;
 0 12 0;
 0 13 0;
 0 14 0;
 0 15 0;
 0 16 0;
 0 17 0;
 0 18 0;
 0 19 0;
 0 20 0;
 3 21 0;
 0 22 0;]';
M_.nstatic = 13;
M_.nfwrd   = 3;
M_.npred   = 3;
M_.nboth   = 0;
M_.nsfwrd   = 3;
M_.nspred   = 3;
M_.ndynamic   = 6;
M_.equations_tags = {
  1 , 'name' , 'Euler equation' ;
  2 , 'name' , 'Variety effects' ;
  3 , 'name' , 'Consumption sector' ;
  4 , 'name' , 'Labor in entry' ;
  5 , 'name' , 'Wage from pricing' ;
  6 , 'name' , 'Labor supply' ;
  7 , 'name' , 'Free entry' ;
  8 , 'name' , 'Profits' ;
  9 , 'name' , 'Aggregate accounting' ;
  10 , 'name' , 'Output' ;
  11 , 'name' , 'Data-consistent output' ;
  12 , 'name' , 'Firm law of motion' ;
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(19, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(9, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 56;
M_.NNZDerivatives(2) = 54;
M_.NNZDerivatives(3) = -1;
M_.params( 1 ) = 0.99;
beta = M_.params( 1 );
M_.params( 2 ) = (1-M_.params(1))/M_.params(1);
r = M_.params( 2 );
M_.params( 3 ) = 0.025;
delta = M_.params( 3 );
M_.params( 4 ) = 0.9241;
chi = M_.params( 4 );
M_.params( 5 ) = 4.0;
psi = M_.params( 5 );
M_.params( 6 ) = 3.8;
epsi = M_.params( 6 );
M_.params( 7 ) = 1.0;
f_e = M_.params( 7 );
M_.params( 8 ) = 0.9;
sigma = M_.params( 8 );
M_.params( 9 ) = 0.979;
rho_z = M_.params( 9 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 18 ) = 1;
oo_.steady_state( 7 ) = 0.87;
oo_.steady_state( 3 ) = (1-M_.params(3))*oo_.steady_state(18)*oo_.steady_state(7)/(M_.params(3)+M_.params(7)*(M_.params(3)+M_.params(2))*(M_.params(6)-1));
oo_.steady_state( 12 ) = oo_.steady_state(3)^(1/(M_.params(6)-1));
oo_.steady_state( 6 ) = oo_.steady_state(18)*oo_.steady_state(12)/(M_.params(6)/(M_.params(6)-1));
oo_.steady_state( 8 ) = oo_.steady_state(7)*(M_.params(3)+M_.params(2))*(M_.params(6)-1)/(M_.params(3)+(M_.params(3)+M_.params(2))*(M_.params(6)-1));
oo_.steady_state( 9 ) = oo_.steady_state(7)-oo_.steady_state(8);
oo_.steady_state( 1 ) = oo_.steady_state(18)*oo_.steady_state(12)*oo_.steady_state(8);
oo_.steady_state( 2 ) = M_.params(3)*oo_.steady_state(3)/(1-M_.params(3));
oo_.steady_state( 4 ) = oo_.steady_state(1)/(M_.params(6)*oo_.steady_state(3));
oo_.steady_state( 5 ) = M_.params(7)*oo_.steady_state(6)/oo_.steady_state(18);
oo_.steady_state( 10 ) = oo_.steady_state(1)+oo_.steady_state(5)*oo_.steady_state(2);
oo_.steady_state( 11 ) = oo_.steady_state(10)/oo_.steady_state(12);
oo_.steady_state( 13 ) = 100*log(oo_.steady_state(10));
oo_.steady_state( 14 ) = 100*log(oo_.steady_state(11));
oo_.steady_state( 15 ) = 100*log(oo_.steady_state(7));
oo_.steady_state( 16 ) = 100*log(oo_.steady_state(8));
oo_.steady_state( 17 ) = 100*log(oo_.steady_state(9));
oo_.steady_state( 19 ) = 100*log(oo_.steady_state(18));
if M_.exo_nbr > 0;
	oo_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];
end;
if M_.exo_det_nbr > 0;
	oo_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];
end;
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 0.0001;
resid;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 40;
options_.order = 2;
var_list_=[];
var_list_ = 'log_Z';
var_list_ = char(var_list_, 'log_Y');
var_list_ = char(var_list_, 'log_Y_r');
var_list_ = char(var_list_, 'log_L');
var_list_ = char(var_list_, 'log_L_C');
var_list_ = char(var_list_, 'log_L_E');
var_list_ = char(var_list_, 'N');
var_list_ = char(var_list_, 'rho');
info = stoch_simul(var_list_);
save('BGM_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('BGM_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('BGM_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('BGM_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('BGM_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
