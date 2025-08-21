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

% Check all of the 35 moments
Data_vs_Model = [Values, model_moments];

% Select the 27 moments that we target
Data_vs_Model_less = Data_vs_Model([1:5,7:11,13:16,18:20,22:23,25,27,29:33,35],:);

res = oo_;
M = M_;
%save('res', 'res');
%save('M', 'M');

save('model_moments', 'model_moments');

posterior_density = res.posterior_density.parameters;
save('posterior_density', 'posterior_density')


posterior_mode = res.posterior_mode.parameters;
save('posterior_mode', 'posterior_mode');
posterior_mean = res.posterior_mean.parameters;
save('posterior_mean', 'posterior_mean');

priors = res.prior;
save("priors", "priors")

posterior_std = res.posterior_std.parameters;
save("posterior_std", "posterior_std");

posterior_HPDlow = res.posterior_hpdinf.parameters;
save('posterior_HPDlow', 'posterior_HPDlow');

posterior_HPDhigh = res.posterior_hpdsup.parameters;
save('posterior_HPDhigh', 'posterior_HPDhigh');

irf = res.irfs;
save('irf.mat', 'irf');

clearvars -except Data_vs*;