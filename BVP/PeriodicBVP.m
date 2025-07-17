function [U, lambda, R, D, dU] = PeriodicBVP(msh, P, Q, F, BC, numCells, opts)
% function [U, lambda, dU] = PeriodicBVP(msh, P, Q, F, BC, numCells)
%
% This function computes a Lagrange P1 Finite Elements approximation
% of the solution of the regular Sturm Liouville equation
%
% - (P(x) * u')' + Q(x) * u = F(x),  ±x > 0,
%
% Here P, Q, and F are periodic coefficients.
% The ODE is combined with a boundary condition expressed in generic form
%
%   A*u(0) + D*u'(0) = B,
%
% where A, D, and B are coefficients.
%
% INPUT: msh (meshObject) describes the mesh,
%        P (function), the 2-order term coefficient,
%        Q (function), the 0-order term coefficient,
%        F (function), the source term,
%        BC is a structure which describes the boundary condition.
%           BC contains 3 fields:
%              + A, the trace coefficient
%              + D, the normal trace coefficient
%              + b, the boundary rhs
%        numCells, the number of cells on which the solution is computed
%        opts (struct) contains options:
%              + verbose (bool)
%              + compute_RtR (bool), to decide whether one should compute
%                the 'classical' RtR coefficient as lambda. Requires 
%                D to be ±1 and b to be 1.
%
% OUTPUT: U  (msh.numNodes x numCells) the solution of the Sturm-Liouville
%                                       boundary value problem
%         lambda (scalar), the transmission coefficient
%         dU (msh.numNodes x numCells) a weak approximation of the
%                                       derivative of U
%
% NOTE: For essential boundary conditions, we use elimination

if (nargin < 7)
  opts = struct('verbose', false, 'compute_RtR', false);
end

if ~isfield(opts, 'verbose'),     opts.verbose = false;     end
if ~isfield(opts, 'compute_RtR'), opts.compute_RtR = false; end

% Solve the local cell problems
% These are periodic cell problems
% The boundary conditions are similar to the one
% of the initial half-line problem

% Define the boundary conditions
cellsBC.A = diag([BC.A,  BC.A]);
cellsBC.D = diag([BC.D, -BC.D]);

% Compute the solutions of the local cell problems
cellsBC.b = [1; 0];
e0 = SturmLiouvilleBVP(msh, P, Q, F, cellsBC);

cellsBC.b = [0; 1];
e1 = SturmLiouvilleBVP(msh, P, Q, F, cellsBC);

% Compute the local transmission operators defined by
% t00 = conj(D)*e0(0) - conj(A)*e0'(0)
% t10 = conj(D)*e1(0) - conj(A)*e1'(0)
% t01 = conj(D)*e0(L) + conj(A)*e0'(L)
% t11 = conj(D)*e1(L) + conj(A)*e1'(L)
% We use strong or weak evaluation depending on whether
% the boundary condition involves the derivative (Neumann or Robin)
% or not (Dirichlet)
if (abs(BC.D) > eps)

  % Neumann or Robin boundary conditions
  alpha = (conj(BC.D) + abs(BC.A)^2 /BC.D);

  % Use strong evaluation
  t00 = alpha * e0(msh.boundsIds(1)) - (conj(BC.A)/BC.D);
  t10 = alpha * e1(msh.boundsIds(1));
  t01 = alpha * e0(msh.boundsIds(2));
  t11 = alpha * e1(msh.boundsIds(2)) - (conj(BC.A)/BC.D);

else

  % Dirichlet boundary conditions
  [~, ~, AA] = SturmLiouvilleBVP(msh, P, Q, F, cellsBC);

  % Use weak evaluation
  t00 = abs(BC.A)^2 * e0' * AA * e0;
  t10 = abs(BC.A)^2 * e0' * AA * e1;
  t01 = abs(BC.A)^2 * e1' * AA * e0;
  t11 = abs(BC.A)^2 * e1' * AA * e1;

end

if (opts.verbose)
  fprintf('t00: %d\t t10: %d\t t01: %d\t t11: %d\n', t00, t10, t01, t11);
end

% Compute the propagation coefficients

% Coefficients of the "Riccati" equation
% Ejk = [ej(k*L); ej'(k*L)]
E00 = [BC.A,  BC.D; conj(BC.D), -conj(BC.A)] \ [1; t00];
E10 = [BC.A,  BC.D; conj(BC.D), -conj(BC.A)] \ [0; t10];
E01 = [BC.A, -BC.D; conj(BC.D),  conj(BC.A)] \ [0; t01];
E11 = [BC.A, -BC.D; conj(BC.D),  conj(BC.A)] \ [1; t11];

% Solve the Riccati equation
[eigenVec, eigenVal] = eig([E01, E11], [E00, E10]);

solsP = diag(eigenVal);
solsD = eigenVec(2, :) ./ eigenVec(1, :);

% Choose the unique solution of the Riccati equation
tol = 1e-5;
tolrad = 1e-5;

if (length(find(abs(solsP) < 1-tolrad)) == 1) %#ok

  % Exactly one of the eigenvalues have a modulus strictly less than 1
  % This corresponds to the good solution
  IGood = find(abs(solsP) < 1-tolrad);

elseif (abs(abs(solsP) - 1) < tol)

  % Both the eigenvalues have a modulus equal to 1
  flux = imag((E00(2) + solsD * E10(2)) .* conj(E00(1) + solsD * E10(1)));

  % The good solution is the one with a negative flux
  IGood = find(sign(flux) == -1);

  % Verifications
  if (length(IGood) ~= 1)
    if (abs(solsP(1) - solsP(2)) < 1e-4)
      % Both solutions are equal: no ambiguity
      IGood = 1;
    else
      figure;
      set(groot,'defaultAxesTickLabelInterpreter','latex');
      set(groot,'defaulttextinterpreter','latex');
      set(groot,'defaultLegendInterpreter','latex');
      tvar = linspace(0, 2*pi, 100);
      plot(cos(tvar), sin(tvar), 'k'); hold on;
      plot(real(solsP), imag(solsP), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
      xlabel('$\Re (p)$');
      ylabel('$\Im (p)$');
      set(gca, 'DataAspectRatio', [1 1 1], 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'FontSize', 16);

      error('Solutions of the Riccati system are distinct despite having the same flux sign. The fluxes are %d and %d. The difference in absolute value between the solutions: %d', flux(1), flux(2), abs(solsP(1) - solsP(2)));
    end
  end

else

  error(['Probleme d''admissibilite : les solutions du systeme ', ...
         'de Riccati valent %0.2e et %0.2e. La difference en ',... 
         'valeur absolue: %d'], ...
         solsP(1), solsP(2), abs(solsP(1) - solsP(2)));

end

R = solsP(IGood);
D = solsD(IGood);

% Construct the solution cell by cell
U = BC.b * (e0 + D * e1) * R.^(0:numCells-1);

% Compute the transmission coefficient
% lambda = conj(D)*u(0) - conj(A)*u'(0)
lambda = t00 + D * t10;

if (abs(abs(BC.D) - 1.0) < eps && abs(BC.b - 1.0) < eps && opts.compute_RtR)
  % Change lambda to the usual expression of the RtR
  % coefficient if requested by the user
  % ONLY FOR BC.D = ±1 AND BC.b = 1.
  lambda = 2*BC.A * U(msh.boundsIds(1), 1) - BC.b;

end

% Construct a weak approximation of the derivative
[MM, ~, DD] = FEmatrices(msh, @(x) ones(size(x)), @(x) ones(size(x)));
dU = MM \ (DD * U);
