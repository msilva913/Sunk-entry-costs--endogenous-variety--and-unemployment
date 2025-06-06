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