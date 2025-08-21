clear; clc;
%close all;

delete('*.asv');

addpath C:\dynare\4.4.3\matlab

firstrun = 0;

if firstrun == 1
    mh;
    data;
end

dynare model;

% Check all of the 35 moments
% Data_vs_Model = [Values', model_moments];