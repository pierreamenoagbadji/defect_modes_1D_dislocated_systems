% TEST_STURM_LIOUVILLE_BVP
% The goal of this code is to test the
% SturmLiouvilleBVP.m code

% ************** %
% Initialization %
% ************** %
clear; clc;

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
coul = [035, 120, 057];

% ********** %
% Parameters %
% ********** %
% Coefficients
funP = @(x) ones(size(x, 1), 1);
funQ = @(x) 2 + cos(2*pi*x);

