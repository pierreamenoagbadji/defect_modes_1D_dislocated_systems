% TEST_STURM_LIOUVILLE_BVP
% The goal of this code is to test the
% SturmLiouvilleBVP.m code

% ************** %
% Initialization %
% ************** %
clear; clc;

% ********** %
% Parameters %
% ********** %
% Coefficients and source term
funP = @(x) ones(size(x, 1), 1);
funQ = @(x) 2 + cos(2*pi*x);
RHS  = @(x) cos(2*pi*x) + sin(2*pi*x);

% Boundary conditions
BC.A = eye(2);
BC.D = zeros(2);
BC.b = zeros(2, 1);

% Domain and mesh
xmin = 0;
xmax = 1;
N = 32;
msh = meshObject('uniform', xmin, xmax, N);

% *************** %
% Solution of BVP %
% *************** %
% Solve the boundary value problem
U = SturmLiouvilleBVP(msh, funP, funQ, RHS, BC);

% Compute the solution
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
vert = [035, 120, 057] / 255;
plot(msh.nodes, real(U), 'color', vert);
set(gca, 'FontSize', 16);

