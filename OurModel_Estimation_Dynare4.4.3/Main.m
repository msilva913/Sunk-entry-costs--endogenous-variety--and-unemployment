clear; clc;
close all;

delete('*.asv');

addpath C:\dynare\4.4.3\matlab
firstrun = 1;

if firstrun == 1
    mh;
    data;
end

dynare model;

Data_vs_Model = [Values,model_moments];

%%% Save results %%%
save("model_moments", "model_moments") % for table
res = oo_;
M = M_;
save('res', 'res');
save('M', 'M');

posterior_density = res.posterior_density.parameters;
save('posterior_density', 'posterior_density')


posterior_mode = res.posterior_mode.parameters;
save('posterior_mode', 'posterior_mode');
posterior_mean = res.posterior_mean.parameters;
save('posterior_mean', 'posterior_mean');

irf = res.irfs;
save('irf.mat', 'irf');

% Load saveplots
for i_chain = 1:56
    mhname = sprintf('model_mh%d_blck1.mat', i_chain);
    addpath('model\metropolis');
    load(mhname);
    if i_chain == 1
        Traceplot_parameter = x2;
        Traceplot_likelihood = logpo2;
    else
        Traceplot_parameter = [Traceplot_parameter; x2];
        Traceplot_likelihood = [Traceplot_likelihood; logpo2];
    end
end
%mh_mode = mode(Traceplot_parameter)

subplot(3,4,1),
plot(Traceplot_likelihood(1:end));
title('log likelihood');

save("Traceplot_parameter", "Traceplot_parameter")

for i_plot = 1:11
    subplot(3,4,i_plot+1),
    plot(Traceplot_parameter(1:end,i_plot));
    title(parameter_names(i_plot),'Interpreter','none');
end

