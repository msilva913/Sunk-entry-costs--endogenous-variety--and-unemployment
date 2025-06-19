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

clearvars -except Data_vs*;