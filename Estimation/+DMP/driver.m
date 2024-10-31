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
M_.fname = 'DMP';
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
M_.exo_names(1) = {'e_x'};
M_.exo_names_tex(1) = {'{e_x}'};
M_.exo_names_long(1) = {'technology shock'};
M_.endo_names = cell(10,1);
M_.endo_names_tex = cell(10,1);
M_.endo_names_long = cell(10,1);
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
M_.endo_names(7) = {'theta_x'};
M_.endo_names_tex(7) = {'{\theta_x}'};
M_.endo_names_long(7) = {'technology shock'};
M_.endo_names(8) = {'u_imp'};
M_.endo_names_tex(8) = {'u\_imp'};
M_.endo_names_long(8) = {'u_imp'};
M_.endo_names(9) = {'v_imp'};
M_.endo_names_tex(9) = {'v\_imp'};
M_.endo_names_long(9) = {'v_imp'};
M_.endo_names(10) = {'theta_imp'};
M_.endo_names_tex(10) = {'theta\_imp'};
M_.endo_names_long(10) = {'theta_imp'};
M_.endo_partitions = struct();
M_.param_names = cell(12,1);
M_.param_names_tex = cell(12,1);
M_.param_names_long = cell(12,1);
M_.param_names(1) = {'x'};
M_.param_names_tex(1) = {'{x}'};
M_.param_names_long(1) = {'Steady-state productivity'};
M_.param_names(2) = {'beta'};
M_.param_names_tex(2) = {'{\beta}'};
M_.param_names_long(2) = {'Discount factor'};
M_.param_names(3) = {'delta'};
M_.param_names_tex(3) = {'{\delta}'};
M_.param_names_long(3) = {'Product destruction rate'};
M_.param_names(4) = {'s'};
M_.param_names_tex(4) = {'{s}'};
M_.param_names_long(4) = {'Worker separation rate'};
M_.param_names(5) = {'fbar'};
M_.param_names_tex(5) = {'fbar'};
M_.param_names_long(5) = {'fbar'};
M_.param_names(6) = {'qbar'};
M_.param_names_tex(6) = {'{\overline{q}}'};
M_.param_names_long(6) = {'Mean vacancy filling rate'};
M_.param_names(7) = {'alpha_L'};
M_.param_names_tex(7) = {'\alpha_L'};
M_.param_names_long(7) = {'Worker bargaining power'};
M_.param_names(8) = {'b'};
M_.param_names_tex(8) = {'b'};
M_.param_names_long(8) = {'Unemployment value'};
M_.param_names(9) = {'phi'};
M_.param_names_tex(9) = {'\phi'};
M_.param_names_long(9) = {'Elasticity of matching function with respect to unemployment'};
M_.param_names(10) = {'rho_x'};
M_.param_names_tex(10) = {'{\rho_x}'};
M_.param_names_long(10) = {'persistence technology shock'};
M_.param_names(11) = {'rho_b'};
M_.param_names_tex(11) = {'{\rho_Z}'};
M_.param_names_long(11) = {'persistence discount factor shock'};
M_.param_names(12) = {'rho_alphaL'};
M_.param_names_tex(12) = {'{\rho_{\alpha_L}}'};
M_.param_names_long(12) = {'persistence bargaining power shock'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 10;
M_.param_nbr = 12;
M_.orig_endo_nbr = 10;
M_.aux_vars = [];
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
M_.nonzero_hessian_eqs = [1 2 3 4 5 6 7 8 9];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 10;
M_.eq_nbr = 10;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 4 14;
 0 5 0;
 1 6 0;
 0 7 0;
 2 8 0;
 0 9 15;
 3 10 16;
 0 11 0;
 0 12 0;
 0 13 0;]';
M_.nstatic = 5;
M_.nfwrd   = 2;
M_.npred   = 2;
M_.nboth   = 1;
M_.nsfwrd   = 3;
M_.nspred   = 3;
M_.ndynamic   = 5;
M_.dynamic_tmp_nbr = [5; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'Job creation condition' ;
  2 , 'name' , 'Wage equation' ;
  3 , 'name' , 'LOM of unemployment' ;
  4 , 'name' , 'Market tightness' ;
  5 , 'name' , 'Job finding probability' ;
  6 , 'name' , 'Vacancy filling probability' ;
  7 , 'name' , 'u_imp' ;
  8 , 'name' , 'v_imp' ;
  9 , 'name' , 'theta_imp' ;
  10 , 'name' , 'Labor productivity process' ;
};
M_.mapping.w.eqidx = [1 2 ];
M_.mapping.v.eqidx = [4 8 ];
M_.mapping.u.eqidx = [3 4 7 ];
M_.mapping.theta.eqidx = [2 4 5 6 9 ];
M_.mapping.f.eqidx = [3 5 ];
M_.mapping.q.eqidx = [1 6 ];
M_.mapping.theta_x.eqidx = [1 2 10 ];
M_.mapping.u_imp.eqidx = [7 ];
M_.mapping.v_imp.eqidx = [8 ];
M_.mapping.theta_imp.eqidx = [9 ];
M_.mapping.e_x.eqidx = [10 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [3 5 7 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(10, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(12, 1);
M_.endo_trends = struct('deflator', cell(10, 1), 'log_deflator', cell(10, 1), 'growth_factor', cell(10, 1), 'log_growth_factor', cell(10, 1));
M_.NNZDerivatives = [26; 14; -1; ];
M_.static_tmp_nbr = [2; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(1) = 1.0;
x = M_.params(1);
M_.params(2) = 0.9966722160545233;
beta = M_.params(2);
M_.params(7) = 0.5;
alpha_L = M_.params(7);
M_.params(8) = 0.92;
b = M_.params(8);
M_.params(9) = 0.5;
phi = M_.params(9);
M_.params(4) = 0.035;
s = M_.params(4);
M_.params(5) = 0.4545454545454545;
fbar = M_.params(5);
M_.params(6) = 0.8024691358024691;
qbar = M_.params(6);
M_.params(10) = 0.979;
rho_x = M_.params(10);
stdx = 0.007;
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
options_.order = 2;
options_.periods = 0;
options_.pruning = true;
var_list_ = {'theta_x';'u_imp';'v_imp';'theta_imp'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'DMP_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'DMP_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'DMP_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'DMP_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'DMP_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'DMP_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'DMP_results.mat'], 'oo_recursive_', '-append');
end
disp('Note: 2 warning(s) encountered in the preprocessor')
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
