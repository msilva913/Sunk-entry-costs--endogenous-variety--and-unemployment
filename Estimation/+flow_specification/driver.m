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
M_.fname = 'flow_specification';
M_.dynare_version = '5.4';
oo_.dynare_version = '5.4';
options_.dynare_version = '5.4';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'e_z'};
M_.exo_names_tex(1) = {'{e_z}'};
M_.exo_names_long(1) = {'technology shock'};
M_.endo_names = cell(32,1);
M_.endo_names_tex = cell(32,1);
M_.endo_names_long = cell(32,1);
M_.endo_names(1) = {'w'};
M_.endo_names_tex(1) = {'{w}'};
M_.endo_names_long(1) = {'wage'};
M_.endo_names(2) = {'v'};
M_.endo_names_tex(2) = {'{v}'};
M_.endo_names_long(2) = {'vacancy stock'};
M_.endo_names(3) = {'u'};
M_.endo_names_tex(3) = {'{u}'};
M_.endo_names_long(3) = {'unemployment'};
M_.endo_names(4) = {'theta'};
M_.endo_names_tex(4) = {'{theta}'};
M_.endo_names_long(4) = {'market tightness'};
M_.endo_names(5) = {'f'};
M_.endo_names_tex(5) = {'{f}'};
M_.endo_names_long(5) = {'job finding rate'};
M_.endo_names(6) = {'q'};
M_.endo_names_tex(6) = {'{q}'};
M_.endo_names_long(6) = {'vacancy filling rate'};
M_.endo_names(7) = {'N'};
M_.endo_names_tex(7) = {'{N}'};
M_.endo_names_long(7) = {'number of firms'};
M_.endo_names(8) = {'N_e'};
M_.endo_names_tex(8) = {'{N_e}'};
M_.endo_names_long(8) = {'number of new businesses'};
M_.endo_names(9) = {'p'};
M_.endo_names_tex(9) = {'{p}'};
M_.endo_names_long(9) = {'relative price'};
M_.endo_names(10) = {'w_R'};
M_.endo_names_tex(10) = {'{w_R}'};
M_.endo_names_long(10) = {'recruiter payment'};
M_.endo_names(11) = {'nu_f'};
M_.endo_names_tex(11) = {'{nu_f}'};
M_.endo_names_long(11) = {'retailer firm value'};
M_.endo_names(12) = {'d_f'};
M_.endo_names_tex(12) = {'{d_f}'};
M_.endo_names_long(12) = {'dividend'};
M_.endo_names(13) = {'Y_c'};
M_.endo_names_tex(13) = {'{Y_c}'};
M_.endo_names_long(13) = {'retail output'};
M_.endo_names(14) = {'C'};
M_.endo_names_tex(14) = {'{C}'};
M_.endo_names_long(14) = {'consumption'};
M_.endo_names(15) = {'lam'};
M_.endo_names_tex(15) = {'{\lambda}'};
M_.endo_names_long(15) = {'Marginal utility of consumption'};
M_.endo_names(16) = {'Y'};
M_.endo_names_tex(16) = {'{Y}'};
M_.endo_names_long(16) = {'aggregate output'};
M_.endo_names(17) = {'L'};
M_.endo_names_tex(17) = {'{L}'};
M_.endo_names_long(17) = {'total labor'};
M_.endo_names(18) = {'L_e'};
M_.endo_names_tex(18) = {'{L_e}'};
M_.endo_names_long(18) = {'Labor in entry'};
M_.endo_names(19) = {'L_c'};
M_.endo_names_tex(19) = {'{L_c}'};
M_.endo_names_long(19) = {'Labor in consumption'};
M_.endo_names(20) = {'labor_share_obs'};
M_.endo_names_tex(20) = {'labor\_share\_obs'};
M_.endo_names_long(20) = {'labor_share_obs'};
M_.endo_names(21) = {'theta_z'};
M_.endo_names_tex(21) = {'{\theta_z}'};
M_.endo_names_long(21) = {'technology shock'};
M_.endo_names(22) = {'z_obs_m'};
M_.endo_names_tex(22) = {'z\_obs\_m'};
M_.endo_names_long(22) = {'z_obs_m'};
M_.endo_names(23) = {'z_obs'};
M_.endo_names_tex(23) = {'z\_obs'};
M_.endo_names_long(23) = {'z_obs'};
M_.endo_names(24) = {'u_obs'};
M_.endo_names_tex(24) = {'u\_obs'};
M_.endo_names_long(24) = {'u_obs'};
M_.endo_names(25) = {'v_obs'};
M_.endo_names_tex(25) = {'v\_obs'};
M_.endo_names_long(25) = {'v_obs'};
M_.endo_names(26) = {'u_imp'};
M_.endo_names_tex(26) = {'u\_imp'};
M_.endo_names_long(26) = {'u_imp'};
M_.endo_names(27) = {'v_imp'};
M_.endo_names_tex(27) = {'v\_imp'};
M_.endo_names_long(27) = {'v_imp'};
M_.endo_names(28) = {'theta_imp'};
M_.endo_names_tex(28) = {'theta\_imp'};
M_.endo_names_long(28) = {'theta_imp'};
M_.endo_names(29) = {'N_imp'};
M_.endo_names_tex(29) = {'N\_imp'};
M_.endo_names_long(29) = {'N_imp'};
M_.endo_names(30) = {'N_e_imp'};
M_.endo_names_tex(30) = {'N\_e\_imp'};
M_.endo_names_long(30) = {'N_e_imp'};
M_.endo_names(31) = {'C_imp'};
M_.endo_names_tex(31) = {'C\_imp'};
M_.endo_names_long(31) = {'C_imp'};
M_.endo_names(32) = {'AUX_ENDO_LAG_21_1'};
M_.endo_names_tex(32) = {'AUX\_ENDO\_LAG\_21\_1'};
M_.endo_names_long(32) = {'AUX_ENDO_LAG_21_1'};
M_.endo_partitions = struct();
M_.endo_partitions.lone_name = { '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' 'Labor share' '' '' '' '' '' '' '' '' '' '' '' '' };
M_.param_names = cell(13,1);
M_.param_names_tex = cell(13,1);
M_.param_names_long = cell(13,1);
M_.param_names(1) = {'w_ss'};
M_.param_names_tex(1) = {'{w_{ss}}'};
M_.param_names_long(1) = {'Steady-state wage'};
M_.param_names(2) = {'beta'};
M_.param_names_tex(2) = {'{\beta}'};
M_.param_names_long(2) = {'Discount factor'};
M_.param_names(3) = {'delta'};
M_.param_names_tex(3) = {'{\delta}'};
M_.param_names_long(3) = {'Product destruction rate'};
M_.param_names(4) = {'mu_net'};
M_.param_names_tex(4) = {'{\mu-1}'};
M_.param_names_long(4) = {'net markup'};
M_.param_names(5) = {'sigma'};
M_.param_names_tex(5) = {'{\sigma}'};
M_.param_names_long(5) = {'inverse of intertemporal elasticity of substitution'};
M_.param_names(6) = {'tau'};
M_.param_names_tex(6) = {'{\tau}'};
M_.param_names_long(6) = {'Aggregate separation rate'};
M_.param_names(7) = {'fbar'};
M_.param_names_tex(7) = {'fbar'};
M_.param_names_long(7) = {'fbar'};
M_.param_names(8) = {'qbar'};
M_.param_names_tex(8) = {'{\overline{q}}'};
M_.param_names_long(8) = {'Mean vacancy filling rate'};
M_.param_names(9) = {'b_ratio'};
M_.param_names_tex(9) = {'b'};
M_.param_names_long(9) = {'Replacement ratio: unemployment benefits to wages'};
M_.param_names(10) = {'eta_L'};
M_.param_names_tex(10) = {'\eta_L'};
M_.param_names_long(10) = {'Elasticity of matching function with respect to unemployment'};
M_.param_names(11) = {'labor_share'};
M_.param_names_tex(11) = {'labor\_share'};
M_.param_names_long(11) = {'labor_share'};
M_.param_names(12) = {'N_ss'};
M_.param_names_tex(12) = {'{N_{ss}}'};
M_.param_names_long(12) = {'N_ss'};
M_.param_names(13) = {'rho_z'};
M_.param_names_tex(13) = {'{\rho_z}'};
M_.param_names_long(13) = {'persistence technology shock'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 32;
M_.param_nbr = 13;
M_.orig_endo_nbr = 31;
M_.aux_vars(1).endo_index = 32;
M_.aux_vars(1).type = 1;
M_.aux_vars(1).orig_index = 22;
M_.aux_vars(1).orig_lead_lag = -1;
M_.aux_vars(1).orig_expr = 'z_obs_m(-1)';
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
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
M_.orig_eq_nbr = 31;
M_.eq_nbr = 32;
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
 0 8 40;
 0 9 0;
 1 10 0;
 0 11 0;
 2 12 0;
 0 13 41;
 3 14 0;
 4 15 0;
 0 16 0;
 0 17 42;
 0 18 43;
 0 19 44;
 0 20 0;
 0 21 0;
 0 22 45;
 0 23 0;
 0 24 0;
 0 25 0;
 0 26 0;
 0 27 0;
 5 28 0;
 6 29 0;
 0 30 0;
 0 31 0;
 0 32 0;
 0 33 0;
 0 34 0;
 0 35 0;
 0 36 0;
 0 37 0;
 0 38 0;
 7 39 0;]';
M_.nstatic = 19;
M_.nfwrd   = 6;
M_.npred   = 7;
M_.nboth   = 0;
M_.nsfwrd   = 6;
M_.nspred   = 7;
M_.ndynamic   = 13;
M_.dynamic_tmp_nbr = [12; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'Job creation condition' ;
  2 , 'name' , 'Marginal revenue product' ;
  3 , 'name' , 'Wage equation' ;
  4 , 'name' , 'Market tightness' ;
  5 , 'name' , 'Job finding probability' ;
  6 , 'name' , 'Vacancy filling probability' ;
  7 , 'name' , 'Aggregate labor' ;
  8 , 'name' , 'Composition of labor' ;
  9 , 'name' , 'Relative price' ;
  10 , 'name' , 'Lagrangian multiplier' ;
  11 , 'name' , 'Firm value' ;
  12 , 'name' , 'Retail output: resources' ;
  13 , 'name' , 'Retail output: expenditure' ;
  14 , 'Name' , 'Business entrants' ;
  14 , 'name' , 'N_e' ;
  15 , 'Name' , 'Firm value relative to price' ;
  15 , 'name' , 'nu_f' ;
  16 , 'Name' , 'Output: expenditure' ;
  16 , 'name' , 'Y' ;
  17 , 'Name' , 'Output: income' ;
  17 , 'name' , '17' ;
  18 , 'name' , 'LOM of unemployment' ;
  19 , 'name' , 'LOM of firms' ;
  20 , 'name' , 'labor share' ;
  21 , 'name' , 'Aggregate separation rate' ;
  22 , 'name' , 'z_obs' ;
  23 , 'name' , 'u_obs' ;
  24 , 'name' , 'v_obs' ;
  25 , 'name' , 'u_imp' ;
  26 , 'name' , 'v_imp' ;
  27 , 'name' , 'theta_imp' ;
  28 , 'name' , 'N_imp' ;
  29 , 'name' , 'N_e_imp' ;
  30 , 'name' , 'C_imp' ;
  31 , 'name' , 'Labor productivity process' ;
};
M_.mapping.w.eqidx = [1 3 20 ];
M_.mapping.v.eqidx = [4 13 24 26 ];
M_.mapping.u.eqidx = [4 7 18 23 25 ];
M_.mapping.theta.eqidx = [3 4 5 6 27 ];
M_.mapping.f.eqidx = [5 18 ];
M_.mapping.q.eqidx = [1 6 ];
M_.mapping.N.eqidx = [9 17 19 28 ];
M_.mapping.N_e.eqidx = [14 16 19 29 ];
M_.mapping.p.eqidx = [2 9 12 15 ];
M_.mapping.w_R.eqidx = [1 2 3 17 ];
M_.mapping.nu_f.eqidx = [11 15 16 ];
M_.mapping.d_f.eqidx = [11 17 ];
M_.mapping.Y_c.eqidx = [12 13 16 ];
M_.mapping.C.eqidx = [10 13 30 ];
M_.mapping.lam.eqidx = [1 10 11 ];
M_.mapping.Y.eqidx = [16 17 20 ];
M_.mapping.L.eqidx = [7 8 17 20 ];
M_.mapping.L_e.eqidx = [8 14 ];
M_.mapping.L_c.eqidx = [8 12 ];
M_.mapping.labor_share_obs.eqidx = [20 ];
M_.mapping.theta_z.eqidx = [2 12 14 21 31 ];
M_.mapping.z_obs_m.eqidx = [21 22 ];
M_.mapping.z_obs.eqidx = [22 ];
M_.mapping.u_obs.eqidx = [23 ];
M_.mapping.v_obs.eqidx = [24 ];
M_.mapping.u_imp.eqidx = [25 ];
M_.mapping.v_imp.eqidx = [26 ];
M_.mapping.theta_imp.eqidx = [27 ];
M_.mapping.N_imp.eqidx = [28 ];
M_.mapping.N_e_imp.eqidx = [29 ];
M_.mapping.C_imp.eqidx = [30 ];
M_.mapping.e_z.eqidx = [31 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [3 5 7 8 21 22 32 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(32, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(13, 1);
M_.endo_trends = struct('deflator', cell(32, 1), 'log_deflator', cell(32, 1), 'growth_factor', cell(32, 1), 'log_growth_factor', cell(32, 1));
M_.NNZDerivatives = [92; -1; -1; ];
M_.static_tmp_nbr = [10; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(1) = 1.0;
w_ss = M_.params(1);
M_.params(2) = 0.99673;
beta = M_.params(2);
xi_inv = 1; 
M_.params(3) = 0.00514;
delta = M_.params(3);
M_.params(6) = 0.031;
tau = M_.params(6);
M_.params(7) = 0.41;
fbar = M_.params(7);
M_.params(8) = 0.8;
qbar = M_.params(8);
M_.params(11) = 0.66;
labor_share = M_.params(11);
M_.params(12) = 1.0;
N_ss = M_.params(12);
M_.params(4) = 0.05;
mu_net = M_.params(4);
M_.params(5) = 1.5;
sigma = M_.params(5);
M_.params(9) = 0.71;
b_ratio = M_.params(9);
M_.params(10) = 0.6;
eta_L = M_.params(10);
M_.params(13) = 0.979;
rho_z = M_.params(13);
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
var_list_ = {'u_obs';'v_obs';'z_obs';'theta_z';'u_imp';'v_imp';'theta_imp';'N_imp';'N_e_imp';'C_imp';'labor_share_obs';'L_c';'L_e'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'flow_specification_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'flow_specification_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'flow_specification_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'flow_specification_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'flow_specification_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'flow_specification_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'flow_specification_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
