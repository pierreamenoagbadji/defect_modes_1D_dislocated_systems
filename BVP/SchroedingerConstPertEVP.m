% function [eigval, eigfun] = SchroedingerConstPertBVP(msh, eigRan, funV, RtRcos)
% SCHROEDINGERCONSTPERT.m
% Computes defect modes for the Schroedinger operator
% 
% -u''(x) + V(x) u(x) acting on L²(R),
%
% where V(x) is constant outside a domain Ω_0
% and anything bounded inside.

% ************** %
% Initialization %
% ************** %
clear; clc;

% ********** %
% Parameters %
% ********** %
% Bounds of the interior domain
IDb = [-1 1];

% Potential inside the interior domain
Vint = @(x) -ones(size(x, 1), 1);
funP = @(x)  ones(size(x, 1), 1);

% Potential on each side of the interior domain
% Here the potential is constant oustide the
% interior domain
Vpos = 0;
Vneg = 0;

% Range for searching the eigenvalue
specRange = [-2 0];
numSpec   = 128;
specvec   = linspace(specRange(1), specRange(2), numSpec)';

% Number of eigenvalues
numEigs = 4;

% ******************** %
% Mesh interior domain %
% ******************** %
numNodes = 128;
msh = meshObject('uniform', IDb(1), IDb(2), numNodes);

% **************************** %
% Eigenvalue searching process %
% **************************** %
% We introduce the solutions of the half-line problems
% - (u_±)" + (V_± - E) u_± = 0 on Ω_±,
%       - (u_±)' + r_± u_± = 0 at x = a_±
% Since the potential is constant outside the
% interior domain, these solutions have
% an explicit expression
uNeg = @(x, E, RobinNeg) exp( sqrt(Vneg - E) * (x - IDb(1))) / (sqrt(Vneg - E) + RobinNeg); 
uPos = @(x, E, RobinPos) exp(-sqrt(Vpos - E) * (x - IDb(2))) / (sqrt(Vpos - E) + RobinPos); 

% Associated RtR coefficients
RtRneg = @(E, RobinNeg) (RobinNeg - sqrt(Vneg - E)) / (RobinNeg + sqrt(Vneg - E));
RtRpos = @(E, RobinPos) (RobinPos - sqrt(Vpos - E)) / (RobinPos + sqrt(Vpos - E));

% Initialize eigenvalues and eigenvectors
eigvals = zeros(numSpec, numEigs);
eigvecs = zeros(msh.numPoints, numEigs, numSpec);

for idI = 1:numSpec

  specVar = specvec(idI);

  % Compute Robin coefficient involved in interior EVP
  RobinNeg = 1i;
  RobinPos = 1i;
  BC.D = diag([-1 1]);
  BC.A = diag([RobinNeg * RtRneg(specVar, RobinNeg), RobinPos * RtRpos(specVar, RobinPos)]);

  % Solve eigenvalue problem
  [eigvals(idI, :), eigvecs(:, :, idI)] = SchroedingerEVP(msh, funP, Vint, BC, numEigs);
  
  % Plot progress
  fprintf('Progress: %3d%%\n', round(100*idI/numSpec));

end

plot(specvec, abs(eigvals(:, 1) - specvec));

%
f = @(k) k .* tan(k * IDb(2)) - sqrt(-Vint(0) - k.^2);
k0 = fzero(f, [0, sqrt(-Vint(0))]);
E = -k0^2;