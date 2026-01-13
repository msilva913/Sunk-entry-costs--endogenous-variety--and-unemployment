clc;
clear;
addpath c:\dynare\6.3\matlab

%% BGM model
dynare BGM.mod

%% BGM model without variety effects (set rho = 1)
dynare BGM_remove_variety_effects.mod 

%% RBC model
% comparable calibration
dynare RBC_bas.mod 

%% BGM (translog)
dynare BGM_translog.mod


%% Solution of model at first-order approximation for comparison
% dynare GS_HRW_linear.mod
% save("results_linear.mat", 'M_', 'oo_', 'options_', '-v7.3');
% 
% T = size(oo_.endo_simul, 2);
% 
% endo_names = cellstr(M_.endo_names);
% vars = {'u_l', 'v', 'theta_log', 'theta_x'}; 
% idx = cellfun(@(v) find(strcmp(endo_names, v), 1), vars);
% % Extract relevant simulations
% sim_matrix = oo_.endo_simul(idx, 1:T).'; % T x num_vars
% sim = cell2struct(num2cell(sim_matrix, 1), vars, 2);
% % Prepare variables for analysis (quarterly average, log or exp as needed)
% % Raw quarterly data 
% u_q      = log(month2quarter(exp(sim.u_l)));
% v_q      = log(month2quarter(sim.v));
% theta_q  = log(month2quarter(exp(sim.theta_log/100)));
% %
% save('simvars_linear.mat','u_q','v_q','theta_q','-v7.3');
% 
% %%
% %% Solution of model for different values of delta (third-order perturbation, pruning)
% fbar = 0.45;
% qbar = 0.80;
% theta = fbar/qbar;
% % Solve for destruction values in a loop
% delta_max = 1 - (1-0.034*0.99)^12;
% dest_ann_vals = [0.10, 0.14, 0.196, 0.236, delta_max];
% n = numel(dest_ann_vals);
% solve_models = true;
% if solve_models
%     for i = 1:n
%         % Create alpha string for safe filenames
%         dest_ann = dest_ann_vals(i);
% 
%         % Impute delta and adjust parameter nu_L
%         delta = 1 - (1 - dest_ann)^(1 / 12);
%         f = fbar/(1-delta);
% 
%         % Calibrate nu_L and include in dynare mod file 
%         jf = @(theta, nu_L) theta/(1+theta^( nu_L))^(1/ nu_L);
%         nu_L = fzero(@(nu_L) jf(theta, nu_L) - f, 1.46);
% 
%         % Obtain strings for dest_ann and nu_L in mod file
%         dest_ann_str = strrep(num2str(dest_ann), '.', '_');
% 
%         % Rename mod file
%         modfilename = sprintf('my_model_%s.mod', dest_ann_str);
% 
%         % Read template and substitute placeholder
%         fid = fopen('GS_HRW.mod','r');
%         if fid == -1
%             error('Could not open the template file.');
%         end
%         modtxt = fread(fid,'*char')';
%         fclose(fid);
% 
%         % Replace values of delta and nu_L;
%         modtxt2 = strrep(modtxt, '${dest_ann}', num2str(dest_ann));
%         modtxt2 = strrep(modtxt2, '${nu_L}', num2str(nu_L, '%.5f'));
% 
%         % Write temporary mod file
%         fid = fopen(modfilename, 'w');
%         fwrite(fid, modtxt2);
%         fclose(fid);
% 
%         % Run Dynare
%         dynare(modfilename, 'noclearall');
% 
%         % Save outputs, using same string
%         results_file = sprintf('results_dest_ann_%s.mat', dest_ann_str);
%         save(results_file, 'M_', 'oo_', 'options_', '-v7.3');
% 
%         % Clean up temp files
%         delete([modfilename(1:end-4) '.m']);
%         delete([modfilename(1:end-4) '.log']);
%         delete([modfilename(1:end-4) '_results.mat']);
%         delete(modfilename);
%     end
% end
% %% 
% % Computation of Quarterly Moments across delta specifications
% % Take Quarterly Average of Monthly Data, 
% % Then apply log, HP filter, and calculate moments
% % Verify accuracy through explicit quarterly average of monthly data
% lam = 1e5; % HP filter lambda for quarterly data
% 
% % Variables of interest; exponentiate log u to ensure positivity
% vars = {'u_l', 'v', 'theta_log', 'theta_x'}; 
% results_list = ["results_dest_ann_0_1.mat", "results_dest_ann_0_14.mat", "results_dest_ann_0_196.mat",...
%     "results_dest_ann_0_236.mat", "results_dest_ann_0_33693.mat"];
% num_vars = numel(vars);
% num_files = numel(results_list);
% nrows = 6;
% 
% %Store HP-filtered correlations
% corrs_mat = zeros(nrows, num_files); %
% % Store Hamilton-filtered correlations
% corrs_mat_ham = zeros(nrows, num_files);
% corrs_mat_growth = zeros(nrows, num_files);
% 
% % Dynamic correlations
% dyn_corrs_mat = zeros(4, num_files);
% 
% for i = 1:num_files
%     load(results_list(i), 'M_', 'oo_', 'options_');
% 
%     % Burn-in drop
%     %drop = getfield(options_, 'drop', 0);
%     %if isempty(drop), drop = 0; end
% 
%     % Index variables in oo_.endo_simul using names
%     endo_names = cellstr(M_.endo_names);
%     idx = cellfun(@(v) find(strcmp(endo_names, v), 1), vars);
%     %if any(cellfun(@isempty, idx))
%      %   error('One or more variables not found among endogenous names.');
%     %end
% 
%     T = size(oo_.endo_simul, 2);
%     %cols = (drop + 1):T;
% 
%     % Extract relevant simulations
%     sim_matrix = oo_.endo_simul(idx, 1:T).'; % T x num_vars
%     sim = cell2struct(num2cell(sim_matrix, 1), vars, 2);
% 
%     % Prepare variables for analysis (quarterly average, log or exp as needed)
%     % Raw quarterly data  in logs
%     u_m      = exp(sim.u_l);
%     u_q      = log(month2quarter(u_m));
%     v_m      = sim.v;
%     v_q      = log(month2quarter(v_m));
%     theta_q  = log(month2quarter(exp(sim.theta_log/100)));
%     x        = exp(sim.theta_x);
%     x_q      = log(month2quarter(x));
% 
%     % Save simulated data for baseline specification (BM) and (FR/MS)
%     if i == 1 % Baseline
%         save('simvars.mat','u_q','v_q','theta_q','x_q','-v7.3');
%         save('simvars_monthly.mat', 'u_m', 'v_m', 'x', '-v7.3'); % Note: monthly data in levels
%     end
%     if i == 4 % FR/MS
%        save('simvars_FR_monthly.mat', 'u_m', 'v_m', 'x', '-v7.3'); % Note: monthly data in levels
%     end 
%     % HP filtering
%     u_hp     = hp_filter(u_q, lam);
%     theta_hp = hp_filter(theta_q, lam);
%     v_hp     = hp_filter(v_q, lam);
%     x_hp     = hp_filter(x_q, lam);
% 
%     % Hamilton filter for robustness
%     u_ham = hamilton_filter(u_q, 8, 4);
%     theta_ham = hamilton_filter(theta_q, 8, 4);
%     v_ham = hamilton_filter(v_q, 8, 4);
%     x_ham = hamilton_filter(x_q, 8, 4);
% 
%     u_ham = u_ham(13:end);
%     theta_ham = theta_ham(13:end);
%     v_ham = v_ham(13:end);
%     x_ham = x_ham(13:end);
%     % Growth rates
%     u_growth = diff(u_q); u_growth = u_growth - mean(u_growth);
%     theta_growth = diff(theta_q); theta_growth = theta_growth - mean(theta_growth);
%     v_growth = diff(v_q); v_growth = v_growth - mean(v_growth);
%     x_growth = diff(x_q); x_growth = x_growth - mean(x_growth);
% 
%     % Collect correlation statistics and volatilities
%     corrs_mat(:, i) = [corr(u_hp, x_hp); corr(theta_hp, x_hp); corr(v_hp, x_hp); corr(u_hp, v_hp); std(u_hp); std(u_hp)/std(x_hp)];
%     corrs_mat_ham(:, i) = [corr(u_ham, x_ham); corr(theta_ham, x_ham); corr(v_ham, x_ham); corr(u_ham, v_ham); std(u_ham); std(u_ham)/std(x_ham)];
%     corrs_mat_growth(:, i) = [corr(u_growth, x_growth); corr(theta_growth, x_growth); corr(v_growth, x_growth); corr(u_growth, v_growth); std(u_growth); std(u_growth)/std(x_growth)];
%     % Dynamic correlations
%     maxLag = 4;
%     r = nan(maxLag,1);
%     for k = 1:maxLag
%         vk = v_hp((k+1):end);
%         pk = x_hp(1:end-k);
%         r(k) = corr(vk, pk);   % no NaNs because we trimmed
%     end
%     dyn_corrs_mat(:, i) = r;
% end
% 
% row_labels = {'u', '$\theta$', 'v', '$corr(u, v)$', '$std(u)$', '$std(u)/std(p)$'};
% corrs_tab = array2table(round(corrs_mat, 3, 'significant'), 'RowNames', row_labels);
% %corrs_tab_ham = array2table(round(corrs_mat_ham, 3, 'significant'), 'RowNames', row_labels);
% corrs_tab_growth = array2table(round(corrs_mat_growth, 3, 'significant'), 'RowNames', row_labels);
% disp(corrs_tab)
% %disp(corrs_tab_ham)
% disp(corrs_tab_growth)
% save('correlation_res.mat', 'corrs_tab', 'corrs_tab_growth', '-v7.3');
% 
% row_labels = {'$p_{t-1}$', '$p_{t-2}$', '$p_{t-3}$', '$p_{t-4}$'};
% % Dynamic correlations table 
% dyn_corrs_tab = array2table(round(dyn_corrs_mat, 2, 'significant'), 'RowNames', row_labels);
% disp(dyn_corrs_tab)
% %% Impulse responses 
% %load results_dest_ann_0_06
% results_list = ["results_dest_ann_0_1.mat", "results_dest_ann_0_196.mat",...
%     "results_dest_ann_0_236.mat", "results_dest_ann_0_33693.mat"];
% num_files = numel(results_list);
% 
% compute_irfs = true;
% if compute_irfs
%     for i = 1:num_files 
%          load(results_list(i), 'M_', 'oo_', 'options_');
%         % 
%         shk = 'e_x'; % name must match your shocks block
%         T = options_.irf; % IRF horizon
%         t = 1:T;
% 
%         % Groups of variables per figure
%         group1 = {'u_log','x_log','theta_log'};
%         group2 = {'e_log','v_log'};
% 
%         group1labels = {'u', 'p', 'theta'};
%         group2labels = {'e', 'v'};
% 
%         % Helper to fetch IRF series: oo_.irfs.
%         get_irf = @(v,s) oo_.irfs.(sprintf('%s_%s',v,s));
% 
%         % Group 1: one figure, all variables overlaid
%         figure('Name','IRFs - Group 1 & 2', 'Position', [100 100 900 400]);
%         for subplotIdx = 1:2
%             subplot(1,2,subplotIdx); hold on;
%             if subplotIdx==1
%                 vars = group1; labels = group1labels; ttl = 'IRFs - Group 1';
%             else
%                 vars = group2; labels = group2labels; ttl = 'IRFs - Group 2';
%             end
%             for i=1:numel(vars)
%                 plot(t, get_irf(vars{i},shk), 'LineWidth', 1.5);
%             end
%             if exist('yline','file'); yline(0,'k-'); else; plot([t(1) t(end)], [0 0], 'k-'); end
%             grid on;
%             legend(strrep(labels,'_','_'),'Location','best');
%             title(ttl);
%             xlabel('Periods');
%             hold off;
%         end
%         print('-clipboard','-dbitmap');
%     end
% end




