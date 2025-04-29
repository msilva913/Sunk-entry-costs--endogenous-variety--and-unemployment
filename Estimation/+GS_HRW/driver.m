%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'GS_HRW';
M_.dynare_version = '5.4';
oo_.dynare_version = '5.4';
options_.dynare_version = '5.4';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(3,1);
M_.exo_names_tex = cell(3,1);
M_.exo_names_long = cell(3,1);
M_.exo_names(1) = {'e_x'};
M_.exo_names_tex(1) = {'{e_x}'};
M_.exo_names_long(1) = {'technology shock'};
M_.exo_names(2) = {'e_b'};
M_.exo_names_tex(2) = {'{e_{b}}'};
M_.exo_names_long(2) = {'discount factor shock'};
M_.exo_names(3) = {'e_alphaL'};
M_.exo_names_tex(3) = {'{e_{alpha,L}}'};
M_.exo_names_long(3) = {'bargaining power shock'};
M_.endo_names = cell(23,1);
M_.endo_names_tex = cell(23,1);
M_.endo_names_long = cell(23,1);
M_.endo_names(1) = {'K'};
M_.endo_names_tex(1) = {'{K}'};
M_.endo_names_long(1) = {'expected discounted difference in vacancy value'};
M_.endo_names(2) = {'w'};
M_.endo_names_tex(2) = {'{w}'};
M_.endo_names_long(2) = {'wage'};
M_.endo_names(3) = {'e'};
M_.endo_names_tex(3) = {'{e}'};
M_.endo_names_long(3) = {'new vacancies'};
M_.endo_names(4) = {'v'};
M_.endo_names_tex(4) = {'{v}'};
M_.endo_names_long(4) = {'vacancy stock'};
M_.endo_names(5) = {'u'};
M_.endo_names_tex(5) = {'{u}'};
M_.endo_names_long(5) = {'unemployment'};
M_.endo_names(6) = {'theta'};
M_.endo_names_tex(6) = {'{theta}'};
M_.endo_names_long(6) = {'market tightness'};
M_.endo_names(7) = {'f'};
M_.endo_names_tex(7) = {'{f}'};
M_.endo_names_long(7) = {'job finding rate'};
M_.endo_names(8) = {'q'};
M_.endo_names_tex(8) = {'{q}'};
M_.endo_names_long(8) = {'vacancy filling rate'};
M_.endo_names(9) = {'theta_x'};
M_.endo_names_tex(9) = {'{\theta_x}'};
M_.endo_names_long(9) = {'technology shock'};
M_.endo_names(10) = {'theta_b'};
M_.endo_names_tex(10) = {'{\theta_b}'};
M_.endo_names_long(10) = {'discount factor shock'};
M_.endo_names(11) = {'theta_alphaL'};
M_.endo_names_tex(11) = {'{\theta_{\alpha_L}}'};
M_.endo_names_long(11) = {'bargaining power shock'};
M_.endo_names(12) = {'x_obs_m'};
M_.endo_names_tex(12) = {'x\_obs\_m'};
M_.endo_names_long(12) = {'x_obs_m'};
M_.endo_names(13) = {'x_obs'};
M_.endo_names_tex(13) = {'x\_obs'};
M_.endo_names_long(13) = {'x_obs'};
M_.endo_names(14) = {'w_obs_m'};
M_.endo_names_tex(14) = {'w\_obs\_m'};
M_.endo_names_long(14) = {'w_obs_m'};
M_.endo_names(15) = {'w_obs'};
M_.endo_names_tex(15) = {'w\_obs'};
M_.endo_names_long(15) = {'w_obs'};
M_.endo_names(16) = {'u_obs'};
M_.endo_names_tex(16) = {'u\_obs'};
M_.endo_names_long(16) = {'u_obs'};
M_.endo_names(17) = {'v_obs'};
M_.endo_names_tex(17) = {'v\_obs'};
M_.endo_names_long(17) = {'v_obs'};
M_.endo_names(18) = {'u_imp'};
M_.endo_names_tex(18) = {'u\_imp'};
M_.endo_names_long(18) = {'u_imp'};
M_.endo_names(19) = {'v_imp'};
M_.endo_names_tex(19) = {'v\_imp'};
M_.endo_names_long(19) = {'v_imp'};
M_.endo_names(20) = {'theta_imp'};
M_.endo_names_tex(20) = {'theta\_imp'};
M_.endo_names_long(20) = {'theta_imp'};
M_.endo_names(21) = {'e_imp'};
M_.endo_names_tex(21) = {'e\_imp'};
M_.endo_names_long(21) = {'e_imp'};
M_.endo_names(22) = {'AUX_ENDO_LAG_11_1'};
M_.endo_names_tex(22) = {'AUX\_ENDO\_LAG\_11\_1'};
M_.endo_names_long(22) = {'AUX_ENDO_LAG_11_1'};
M_.endo_names(23) = {'AUX_ENDO_LAG_13_1'};
M_.endo_names_tex(23) = {'AUX\_ENDO\_LAG\_13\_1'};
M_.endo_names_long(23) = {'AUX_ENDO_LAG_13_1'};
M_.endo_partitions = struct();
M_.param_names = cell(13,1);
M_.param_names_tex = cell(13,1);
M_.param_names_long = cell(13,1);
M_.param_names(1) = {'x'};
M_.param_names_tex(1) = {'{x}'};
M_.param_names_long(1) = {'Steady-state productivity'};
M_.param_names(2) = {'beta'};
M_.param_names_tex(2) = {'{\beta}'};
M_.param_names_long(2) = {'Discount factor'};
M_.param_names(3) = {'delta'};
M_.param_names_tex(3) = {'{\delta}'};
M_.param_names_long(3) = {'Product destruction rate'};
M_.param_names(4) = {'tau'};
M_.param_names_tex(4) = {'{\tau}'};
M_.param_names_long(4) = {'Aggregate separation rate'};
M_.param_names(5) = {'fbar'};
M_.param_names_tex(5) = {'fbar'};
M_.param_names_long(5) = {'fbar'};
M_.param_names(6) = {'qbar'};
M_.param_names_tex(6) = {'{\overline{q}}'};
M_.param_names_long(6) = {'Mean vacancy filling rate'};
M_.param_names(7) = {'xi_inv'};
M_.param_names_tex(7) = {'{1/\xi}'};
M_.param_names_long(7) = {'Inverse elasticity of entrants to vacancy value'};
M_.param_names(8) = {'alpha_L'};
M_.param_names_tex(8) = {'\alpha_L'};
M_.param_names_long(8) = {'Worker bargaining power'};
M_.param_names(9) = {'b'};
M_.param_names_tex(9) = {'b'};
M_.param_names_long(9) = {'Unemployment value'};
M_.param_names(10) = {'nu_L'};
M_.param_names_tex(10) = {'\phi'};
M_.param_names_long(10) = {'Elasticity parameter of HRW'};
M_.param_names(11) = {'rho_x'};
M_.param_names_tex(11) = {'{\rho_x}'};
M_.param_names_long(11) = {'persistence technology shock'};
M_.param_names(12) = {'rho_b'};
M_.param_names_tex(12) = {'{\rho_Z}'};
M_.param_names_long(12) = {'persistence discount factor shock'};
M_.param_names(13) = {'rho_alphaL'};
M_.param_names_tex(13) = {'{\rho_{\alpha_L}}'};
M_.param_names_long(13) = {'persistence bargaining power shock'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 23;
M_.param_nbr = 13;
M_.orig_endo_nbr = 21;
M_.aux_vars(1).endo_index = 22;
M_.aux_vars(1).type = 1;
M_.aux_vars(1).orig_index = 12;
M_.aux_vars(1).orig_lead_lag = -1;
M_.aux_vars(1).orig_expr = 'x_obs_m(-1)';
M_.aux_vars(2).endo_index = 23;
M_.aux_vars(2).type = 1;
M_.aux_vars(2).orig_index = 14;
M_.aux_vars(2).orig_lead_lag = -1;
M_.aux_vars(2).orig_expr = 'w_obs_m(-1)';
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(3, 3);
M_.Correlation_matrix = eye(3, 3);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.orig_eq_nbr = 21;
M_.eq_nbr = 23;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 2;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 2;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 2;
M_.lead_lag_incidence = [
 0 13 36;
 1 14 37;
 0 15 38;
 2 16 0;
 3 17 0;
 0 18 0;
 4 19 0;
 5 20 39;
 6 21 40;
 7 22 0;
 8 23 0;
 9 24 0;
 0 25 0;
 10 26 0;
 0 27 0;
 0 28 0;
 0 29 0;
 0 30 0;
 0 31 0;
 0 32 0;
 0 33 0;
 11 34 0;
 12 35 0;]';
M_.nstatic = 9;
M_.nfwrd   = 2;
M_.npred   = 9;
M_.nboth   = 3;
M_.nsfwrd   = 5;
M_.nspred   = 12;
M_.ndynamic   = 14;
M_.dynamic_tmp_nbr = [7; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'Job creation condition' ;
  2 , 'name' , 'Expected discounted difference in vacancy value' ;
  3 , 'name' , 'Wage equation' ;
  4 , 'name' , 'LOM of vacancies' ;
  5 , 'name' , 'LOM of unemployment' ;
  6 , 'name' , 'Market tightness' ;
  7 , 'name' , 'Job finding probability' ;
  8 , 'name' , 'Vacancy filling probability' ;
  9 , 'name' , 'x_obs_m' ;
  10 , 'name' , 'x_obs' ;
  11 , 'name' , 'w_obs_m' ;
  12 , 'name' , 'w_obs' ;
  13 , 'name' , 'u_obs' ;
  14 , 'name' , 'v_obs' ;
  15 , 'name' , 'u_imp' ;
  16 , 'name' , 'v_imp' ;
  17 , 'name' , 'theta_imp' ;
  18 , 'name' , 'e_imp' ;
  19 , 'name' , 'Labor productivity process' ;
  20 , 'name' , 'Discount factor process' ;
  21 , 'name' , 'Bargaining power process' ;
};
M_.mapping.K.eqidx = [1 2 3 ];
M_.mapping.w.eqidx = [1 3 11 ];
M_.mapping.e.eqidx = [2 4 18 ];
M_.mapping.v.eqidx = [4 6 14 16 ];
M_.mapping.u.eqidx = [4 5 6 13 15 ];
M_.mapping.theta.eqidx = [3 6 7 8 17 ];
M_.mapping.f.eqidx = [5 7 ];
M_.mapping.q.eqidx = [1 4 7 8 ];
M_.mapping.theta_x.eqidx = [1 3 9 19 ];
M_.mapping.theta_b.eqidx = [1 2 20 ];
M_.mapping.theta_alphaL.eqidx = [3 21 ];
M_.mapping.x_obs_m.eqidx = [9 10 ];
M_.mapping.x_obs.eqidx = [10 ];
M_.mapping.w_obs_m.eqidx = [11 12 ];
M_.mapping.w_obs.eqidx = [12 ];
M_.mapping.u_obs.eqidx = [13 ];
M_.mapping.v_obs.eqidx = [14 ];
M_.mapping.u_imp.eqidx = [15 ];
M_.mapping.v_imp.eqidx = [16 ];
M_.mapping.theta_imp.eqidx = [17 ];
M_.mapping.e_imp.eqidx = [18 ];
M_.mapping.e_x.eqidx = [19 ];
M_.mapping.e_b.eqidx = [20 ];
M_.mapping.e_alphaL.eqidx = [21 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [2 4 5 7 8 9 10 11 12 14 22 23 ];
M_.exo_names_orig_ord = [1:3];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(23, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(13, 1);
M_.endo_trends = struct('deflator', cell(23, 1), 'log_deflator', cell(23, 1), 'growth_factor', cell(23, 1), 'log_growth_factor', cell(23, 1));
M_.NNZDerivatives = [71; -1; -1; ];
M_.static_tmp_nbr = [7; 1; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(1) = 1.0;
x = M_.params(1);
M_.params(2) = 0.99673;
beta = M_.params(2);
M_.params(7) = 1.0;
xi_inv = M_.params(7);
M_.params(8) = 0.566;
alpha_L = M_.params(8);
M_.params(9) = 0.9;
b = M_.params(9);
M_.params(10) = 1.5857;
nu_L = M_.params(10);
M_.params(3) = 0.00514;
delta = M_.params(3);
M_.params(4) = 0.034;
tau = M_.params(4);
M_.params(5) = 0.4545454545454545;
fbar = M_.params(5);
M_.params(6) = 0.8024691358024691;
qbar = M_.params(6);
M_.params(11) = 0.979;
rho_x = M_.params(11);
M_.params(12) = 0.979;
rho_b = M_.params(12);
M_.params(13) = 0.979;
rho_alphaL = M_.params(13);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 4.900000000000001e-05;
resid;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 80;
options_.nofunctions = true;
options_.order = 1;
options_.periods = 0;
options_.pruning = true;
var_list_ = {'theta_x';'u_imp';'v_imp';'theta_imp';'e_imp'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'GS_HRW_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'GS_HRW_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'GS_HRW_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'GS_HRW_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'GS_HRW_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'GS_HRW_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'GS_HRW_results.mat'], 'oo_recursive_', '-append');
end
disp('Note: 2 warning(s) encountered in the preprocessor')
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
