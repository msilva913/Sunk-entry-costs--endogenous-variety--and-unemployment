clear,clc,

%%

for i_chain = 1:11
    mhname = sprintf('model_mh%d_blck1.mat',i_chain);
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

load('model_mh_mode');


%%

subplot(3,4,1),
plot(Traceplot_likelihood(1:end));
title('log likelihood');

for i_plot = 1:11
    subplot(3,4,i_plot+1),
    plot(Traceplot_parameter(1:end,i_plot));
    title(parameter_names(i_plot),'Interpreter','none');
end