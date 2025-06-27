function [MM, KK, DD] = FEmatrices(msh, massCoeff, stiffnessCoeff, quadRules)
% [MM, KK, DD] = FEMATRICES(msh, massCoeff, stiffnessCoeff, quadRules)
%
% Computes Finite Elements matrices for P1 Lagrange elements.
%
% INPUT: msh (meshObject) describes the mesh
%        massCoeff (function object) is the mass matrix coefficient
%        stiffnessCoeff (function object) is the stiffness matrix coefficient
%        quadRules, a(n optional) structure with fields:
%         + mass.points, quadrature points for the mass matrix
%         + mass.weights, quadrature weights for the mass matrix
%         + stiffness.points, quadrature points for the stiffness matrix
%         + stiffness.weights, quadrature weights for the stiffness matrix
%
% OUTPUT: MM, the mass matrix
%         KK, the stiffness matrix
%         DD, matrix for the weak evaluation of the dervative
%
% NOTE: The current implementation of the FE matrices assumes that the
%       boundary nodes are the first and the last ones

if (nargin < 4)

  % Default quadrature points and weights
  [quadRules.mass.points, quadRules.mass.weights] = gaussLegendre(6, 0, 1);
  [quadRules.stiffness.points, quadRules.stiffness.weights] = ...
                                                    gaussLegendre(4, 0, 1);

end

if (nargin < 3)

  % Default stiffness coefficient
  stiffnessCoeff = @(x) ones(size(x));

end

if (nargin < 2)

  % Default mass coefficient
  massCoeff = @(x) ones(size(x));

end

% Extract some mesh information
N = msh.numPoints;    % Number of nodes
% Len contains the length of each element
Len = msh.points(msh.segments(:, 2))- msh.points(msh.segments(:, 1));

% =========== %
% Mass matrix %
% =========== %
x = quadRules.mass.points;       % Quadrature points
w = quadRules.mass.weights;      % Quadrature weights
Nquad = length(x);          % Number of quadrature points/weights

Mintegral = (Len * ones(1, Nquad)) .*...
            (massCoeff( msh.points(1:N-1)*ones(1, Nquad) + Len*x )) .*...
            (ones(N-1, 1) * w);

MM = sparse( (2:N-1), (2:N-1), Mintegral(1:N-2, :) * (x.^2).',     N, N ) + ...
     sparse( (2:N-1), (2:N-1), Mintegral(2:N-1, :) * ((1-x).^2).', N, N ) + ...
     sparse( (1:N-1), (2:N),   Mintegral(1:N-1, :) * (x.*(1-x)).', N, N ) + ...
     sparse( (2:N),   (1:N-1), Mintegral(1:N-1, :) * (x.*(1-x)).', N, N );

MM(1, 1) = Mintegral(1, :) * ((1-x).^2).';
MM(N, N) = Mintegral(N-1, :) * (x.^2).';

% ================ %
% Stiffness matrix %
% ================ %
x = quadRules.stiffness.points;       % Quadrature points
w = quadRules.stiffness.weights;      % Quadrature weights
Nquad = length(x);               % Number of quadrature points/weights

Kintegral = (((Len * ones(1, Nquad)) .*...
               stiffnessCoeff( msh.points(1:N-1)*ones(1,Nquad) + Len*x )) *...
               w.') ./ (Len.^2);

KK = sparse( (2:N-1), (2:N-1),  Kintegral(1:N-2) + Kintegral(2:N-1), N, N) + ...
     sparse( (2:N),   (1:N-1), -Kintegral, N, N) +...
     sparse( (1:N-1), (2:N),   -Kintegral, N, N);

KK(1, 1) = Kintegral(1);
KK(N, N) = Kintegral(N-1);

% ================= %
% Derivation matrix %
% ================= %
DD = sparse( (1:N-1), (2:N),  0.5, N, N ) + ...
     sparse( (2:N), (1:N-1), -0.5, N, N );

DD(1, 1) = -0.5;
DD(N, N) =  0.5;
