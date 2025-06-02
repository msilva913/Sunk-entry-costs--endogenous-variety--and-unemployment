clear; clc;
close all;

delete('*.asv');

addpath C:\dynare\4.4.3\matlab

firstrun = 0;

if firstrun == 1
    mh;
    data;
end

dynare model;

Data_vs_Model = [Values,model_moments];

%%% Save results %%%
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
