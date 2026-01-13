%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clearvars -global
clear_persistent_variables(fileparts(which('dynare')), false)
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info
options_ = [];
M_.fname = 'RBC_bas';
M_.dynare_version = '6.3';
oo_.dynare_version = '6.3';
options_.dynare_version = '6.3';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'eps_z'};
M_.exo_names_tex(1) = {'{\varepsilon_z}'};
M_.exo_names_long(1) = {'TFP shock'};
M_.endo_names = cell(14,1);
M_.endo_names_tex = cell(14,1);
M_.endo_names_long = cell(14,1);
M_.endo_names(1) = {'y'};
M_.endo_names_tex(1) = {'{y}'};
M_.endo_names_long(1) = {'output'};
M_.endo_names(2) = {'c'};
M_.endo_names_tex(2) = {'{c}'};
M_.endo_names_long(2) = {'consumption'};
M_.endo_names(3) = {'k'};
M_.endo_names_tex(3) = {'{k}'};
M_.endo_names_long(3) = {'capital'};
M_.endo_names(4) = {'l'};
M_.endo_names_tex(4) = {'{l}'};
M_.endo_names_long(4) = {'hours'};
M_.endo_names(5) = {'z'};
M_.endo_names_tex(5) = {'{z}'};
M_.endo_names_long(5) = {'TFP'};
M_.endo_names(6) = {'r'};
M_.endo_names_tex(6) = {'{r}'};
M_.endo_names_long(6) = {'annualized interest rate'};
M_.endo_names(7) = {'w'};
M_.endo_names_tex(7) = {'{w}'};
M_.endo_names_long(7) = {'real wage'};
M_.endo_names(8) = {'invest'};
M_.endo_names_tex(8) = {'{i}'};
M_.endo_names_long(8) = {'investment'};
M_.endo_names(9) = {'log_y'};
M_.endo_names_tex(9) = {'{\log(y)}'};
M_.endo_names_long(9) = {'log output'};
M_.endo_names(10) = {'log_k'};
M_.endo_names_tex(10) = {'{\log(k)}'};
M_.endo_names_long(10) = {'log capital stock'};
M_.endo_names(11) = {'log_c'};
M_.endo_names_tex(11) = {'{\log(c)}'};
M_.endo_names_long(11) = {'log consumption'};
M_.endo_names(12) = {'log_l'};
M_.endo_names_tex(12) = {'{\log(l)}'};
M_.endo_names_long(12) = {'log labor'};
M_.endo_names(13) = {'log_w'};
M_.endo_names_tex(13) = {'{\log(w)}'};
M_.endo_names_long(13) = {'log real wage'};
M_.endo_names(14) = {'log_invest'};
M_.endo_names_tex(14) = {'{\log(i)}'};
M_.endo_names_long(14) = {'log investment'};
M_.endo_partitions = struct();
M_.param_names = cell(7,1);
M_.param_names_tex = cell(7,1);
M_.param_names_long = cell(7,1);
M_.param_names(1) = {'beta'};
M_.param_names_tex(1) = {'{\beta}'};
M_.param_names_long(1) = {'discount factor'};
M_.param_names(2) = {'chi'};
M_.param_names_tex(2) = {'{\chi}'};
M_.param_names_long(2) = {'labor disutility parameter'};
M_.param_names(3) = {'psi'};
M_.param_names_tex(3) = {'{\psi}'};
M_.param_names_long(3) = {'Frisch elasticity of labor supply'};
M_.param_names(4) = {'sigma'};
M_.param_names_tex(4) = {'{\sigma}'};
M_.param_names_long(4) = {'risk aversion'};
M_.param_names(5) = {'delta'};
M_.param_names_tex(5) = {'{\delta}'};
M_.param_names_long(5) = {'depreciation rate'};
M_.param_names(6) = {'alpha'};
M_.param_names_tex(6) = {'{\alpha}'};
M_.param_names_long(6) = {'capital share'};
M_.param_names(7) = {'rhoz'};
M_.param_names_tex(7) = {'{\rho_z}'};
M_.param_names_long(7) = {'persistence TFP shock'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 14;
M_.param_nbr = 7;
M_.orig_endo_nbr = 14;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.learnt_shocks = [];
M_.learnt_endval = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
M_.matched_irfs = {};
M_.matched_irfs_weights = {};
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
options_.ramsey_policy = false;
options_.discretionary_policy = false;
M_.eq_nbr = 14;
M_.ramsey_orig_eq_nbr = 0;
M_.ramsey_orig_endo_nbr = 0;
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
 0 3 0;
 0 4 17;
 1 5 0;
 0 6 18;
 2 7 19;
 0 8 0;
 0 9 0;
 0 10 0;
 0 11 0;
 0 12 0;
 0 13 0;
 0 14 0;
 0 15 0;
 0 16 0;]';
M_.nstatic = 10;
M_.nfwrd   = 2;
M_.npred   = 1;
M_.nboth   = 1;
M_.nsfwrd   = 3;
M_.nspred   = 2;
M_.ndynamic   = 4;
M_.dynamic_tmp_nbr = [6; 2; 0; 0; ];
M_.equations_tags = {
  1 , 'name' , 'Euler equation' ;
  2 , 'name' , 'Labor FOC' ;
  3 , 'name' , 'Law of motion capital' ;
  4 , 'name' , 'resource constraint' ;
  5 , 'name' , 'production function' ;
  6 , 'name' , 'real wage/firm FOC labor' ;
  7 , 'name' , 'annualized real interest rate/firm FOC capital' ;
  8 , 'name' , 'exogenous TFP process' ;
  9 , 'name' , 'Definition log output' ;
  10 , 'name' , 'Definition log capital' ;
  11 , 'name' , 'Definition log consumption' ;
  12 , 'name' , 'Definition log hours' ;
  13 , 'name' , 'Definition log wage' ;
  14 , 'name' , 'Definition log investment' ;
};
M_.mapping.y.eqidx = [4 5 6 7 9 ];
M_.mapping.c.eqidx = [1 2 4 11 ];
M_.mapping.k.eqidx = [1 3 5 7 10 ];
M_.mapping.l.eqidx = [1 2 5 6 12 ];
M_.mapping.z.eqidx = [1 5 8 ];
M_.mapping.r.eqidx = [7 ];
M_.mapping.w.eqidx = [2 6 13 ];
M_.mapping.invest.eqidx = [3 4 14 ];
M_.mapping.log_y.eqidx = [9 ];
M_.mapping.log_k.eqidx = [10 ];
M_.mapping.log_c.eqidx = [11 ];
M_.mapping.log_l.eqidx = [12 ];
M_.mapping.log_w.eqidx = [13 ];
M_.mapping.log_invest.eqidx = [14 ];
M_.mapping.eps_z.eqidx = [8 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.block_structure.time_recursive = false;
M_.block_structure.block(1).Simulation_Type = 1;
M_.block_structure.block(1).endo_nbr = 1;
M_.block_structure.block(1).mfs = 1;
M_.block_structure.block(1).equation = [ 8];
M_.block_structure.block(1).variable = [ 5];
M_.block_structure.block(1).is_linear = true;
M_.block_structure.block(1).NNZDerivatives = 2;
M_.block_structure.block(1).bytecode_jacob_cols_to_sparse = [1 2 ];
M_.block_structure.block(2).Simulation_Type = 8;
M_.block_structure.block(2).endo_nbr = 6;
M_.block_structure.block(2).mfs = 6;
M_.block_structure.block(2).equation = [ 2 4 6 3 5 1];
M_.block_structure.block(2).variable = [ 7 8 1 3 4 2];
M_.block_structure.block(2).is_linear = false;
M_.block_structure.block(2).NNZDerivatives = 19;
M_.block_structure.block(2).bytecode_jacob_cols_to_sparse = [4 7 8 9 10 11 12 17 18 ];
M_.block_structure.block(3).Simulation_Type = 1;
M_.block_structure.block(3).endo_nbr = 7;
M_.block_structure.block(3).mfs = 7;
M_.block_structure.block(3).equation = [ 12 13 14 11 10 9 7];
M_.block_structure.block(3).variable = [ 12 13 14 11 10 9 6];
M_.block_structure.block(3).is_linear = true;
M_.block_structure.block(3).NNZDerivatives = 7;
M_.block_structure.block(3).bytecode_jacob_cols_to_sparse = [8 9 10 11 12 13 14 ];
M_.block_structure.block(1).g1_sparse_rowval = int32([]);
M_.block_structure.block(1).g1_sparse_colval = int32([]);
M_.block_structure.block(1).g1_sparse_colptr = int32([]);
M_.block_structure.block(2).g1_sparse_rowval = int32([4 5 1 3 2 4 2 3 5 4 6 1 3 5 1 2 6 6 6 ]);
M_.block_structure.block(2).g1_sparse_colval = int32([4 4 7 7 8 8 9 9 9 10 10 11 11 11 12 12 12 17 18 ]);
M_.block_structure.block(2).g1_sparse_colptr = int32([1 1 1 1 3 3 3 5 7 10 12 15 18 18 18 18 18 19 20 ]);
M_.block_structure.block(3).g1_sparse_rowval = int32([]);
M_.block_structure.block(3).g1_sparse_colval = int32([]);
M_.block_structure.block(3).g1_sparse_colptr = int32([]);
M_.block_structure.variable_reordered = [ 5 7 8 1 3 4 2 12 13 14 11 10 9 6];
M_.block_structure.equation_reordered = [ 8 2 4 6 3 5 1 12 13 14 11 10 9 7];
M_.block_structure.incidence(1).lead_lag = -1;
M_.block_structure.incidence(1).sparse_IM = [
 3 3;
 5 3;
 7 3;
 8 5;
];
M_.block_structure.incidence(2).lead_lag = 0;
M_.block_structure.incidence(2).sparse_IM = [
 1 2;
 1 3;
 2 2;
 2 4;
 2 7;
 3 3;
 3 8;
 4 1;
 4 2;
 4 8;
 5 1;
 5 4;
 5 5;
 6 1;
 6 4;
 6 7;
 7 1;
 7 6;
 8 5;
 9 1;
 9 9;
 10 3;
 10 10;
 11 2;
 11 11;
 12 4;
 12 12;
 13 7;
 13 13;
 14 8;
 14 14;
];
M_.block_structure.incidence(3).lead_lag = 1;
M_.block_structure.incidence(3).sparse_IM = [
 1 2;
 1 4;
 1 5;
];
M_.block_structure.dyn_tmp_nbr = 9;
M_.state_var = [5 3 ];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(14, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(7, 1);
M_.endo_trends = struct('deflator', cell(14, 1), 'log_deflator', cell(14, 1), 'growth_factor', cell(14, 1), 'log_growth_factor', cell(14, 1));
M_.NNZDerivatives = [39; -1; -1; ];
M_.dynamic_g1_sparse_rowval = int32([3 5 7 8 4 5 6 7 9 1 2 4 11 1 3 10 2 5 6 12 5 8 7 2 6 13 3 4 14 9 10 11 12 13 14 1 1 1 8 ]);
M_.dynamic_g1_sparse_colval = int32([3 3 3 5 15 15 15 15 15 16 16 16 16 17 17 17 18 18 18 18 19 19 20 21 21 21 22 22 22 23 24 25 26 27 28 30 32 33 43 ]);
M_.dynamic_g1_sparse_colptr = int32([1 1 1 4 4 5 5 5 5 5 5 5 5 5 5 10 14 17 21 23 24 27 30 31 32 33 34 35 36 36 37 37 38 39 39 39 39 39 39 39 39 39 39 40 ]);
M_.lhs = {
'c^(-sigma)'; 
'chi*l^(1/psi)'; 
'k'; 
'y'; 
'y'; 
'w'; 
'r'; 
'z'; 
'log_y'; 
'log_k'; 
'log_c'; 
'log_l'; 
'log_w'; 
'log_invest'; 
};
M_.static_tmp_nbr = [5; 2; 0; 0; ];
M_.block_structure_stat.block(1).Simulation_Type = 3;
M_.block_structure_stat.block(1).endo_nbr = 1;
M_.block_structure_stat.block(1).mfs = 1;
M_.block_structure_stat.block(1).equation = [ 8];
M_.block_structure_stat.block(1).variable = [ 5];
M_.block_structure_stat.block(2).Simulation_Type = 6;
M_.block_structure_stat.block(2).endo_nbr = 6;
M_.block_structure_stat.block(2).mfs = 6;
M_.block_structure_stat.block(2).equation = [ 2 3 4 5 6 1];
M_.block_structure_stat.block(2).variable = [ 7 8 2 1 4 3];
M_.block_structure_stat.block(3).Simulation_Type = 1;
M_.block_structure_stat.block(3).endo_nbr = 7;
M_.block_structure_stat.block(3).mfs = 7;
M_.block_structure_stat.block(3).equation = [ 12 13 14 11 10 9 7];
M_.block_structure_stat.block(3).variable = [ 12 13 14 11 10 9 6];
M_.block_structure_stat.variable_reordered = [ 5 7 8 2 1 4 3 12 13 14 11 10 9 6];
M_.block_structure_stat.equation_reordered = [ 8 2 3 4 5 6 1 12 13 14 11 10 9 7];
M_.block_structure_stat.incidence.sparse_IM = [
 1 2;
 1 3;
 1 4;
 1 5;
 2 2;
 2 4;
 2 7;
 3 3;
 3 8;
 4 1;
 4 2;
 4 8;
 5 1;
 5 3;
 5 4;
 5 5;
 6 1;
 6 4;
 6 7;
 7 1;
 7 3;
 7 6;
 8 5;
 9 1;
 9 9;
 10 3;
 10 10;
 11 2;
 11 11;
 12 4;
 12 12;
 13 7;
 13 13;
 14 8;
 14 14;
];
M_.block_structure_stat.tmp_nbr = 7;
M_.block_structure_stat.block(1).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(2).g1_sparse_rowval = int32([1 5 2 3 1 3 6 3 4 5 1 4 5 6 2 4 6 ]);
M_.block_structure_stat.block(2).g1_sparse_colval = int32([1 1 2 2 3 3 3 4 4 4 5 5 5 5 6 6 6 ]);
M_.block_structure_stat.block(2).g1_sparse_colptr = int32([1 3 5 8 11 15 18 ]);
M_.block_structure_stat.block(3).g1_sparse_rowval = int32([]);
M_.block_structure_stat.block(3).g1_sparse_colval = int32([]);
M_.block_structure_stat.block(3).g1_sparse_colptr = int32([]);
M_.static_g1_sparse_rowval = int32([4 5 6 7 9 1 2 4 11 1 3 5 7 10 1 2 5 6 12 1 5 8 7 2 6 13 3 4 14 9 10 11 12 13 14 ]);
M_.static_g1_sparse_colval = int32([1 1 1 1 1 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 6 7 7 7 8 8 8 9 10 11 12 13 14 ]);
M_.static_g1_sparse_colptr = int32([1 6 10 15 20 23 24 27 30 31 32 33 34 35 36 ]);
M_.params(1) = 0.99;
beta = M_.params(1);
M_.params(5) = 0.025;
delta = M_.params(5);
M_.params(2) = 0.9241;
chi = M_.params(2);
M_.params(3) = 4.0;
psi = M_.params(3);
M_.params(4) = 1;
sigma = M_.params(4);
M_.params(6) = 0.33;
alpha = M_.params(6);
M_.params(7) = 0.979;
rhoz = M_.params(7);
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(4) = 1.0;
oo_.steady_state(3) = ((1/M_.params(1)-(1-M_.params(5)))/M_.params(6))^(1/(M_.params(6)-1))*oo_.steady_state(4);
oo_.steady_state(8) = M_.params(5)*oo_.steady_state(3);
oo_.steady_state(1) = oo_.steady_state(3)^M_.params(6)*oo_.steady_state(4)^(1-M_.params(6));
oo_.steady_state(2) = oo_.steady_state(1)-oo_.steady_state(8);
oo_.steady_state(7) = (1-M_.params(6))*oo_.steady_state(1)/oo_.steady_state(4);
oo_.steady_state(6) = oo_.steady_state(1)*M_.params(6)*4/oo_.steady_state(3);
oo_.steady_state(9) = 100*log(oo_.steady_state(1));
oo_.steady_state(10) = 100*log(oo_.steady_state(3));
oo_.steady_state(11) = 100*log(oo_.steady_state(2));
oo_.steady_state(12) = 100*log(oo_.steady_state(4));
oo_.steady_state(13) = 100*log(oo_.steady_state(7));
oo_.steady_state(14) = log(oo_.steady_state(8));
oo_.steady_state(5) = 0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 0.0001;
options_resid_ = struct();
display_static_residuals(M_, options_, oo_, options_resid_);
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 40;
options_.order = 1;
var_list_ = {'log_y';'log_l';'z'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'RBC_bas_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'RBC_bas_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'RBC_bas_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'RBC_bas_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'RBC_bas_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'RBC_bas_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'RBC_bas_results.mat'], 'oo_recursive_', '-append');
end
if exist('options_mom_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'RBC_bas_results.mat'], 'options_mom_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
