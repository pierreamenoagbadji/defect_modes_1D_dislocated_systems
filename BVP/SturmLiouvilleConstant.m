function U = SturmLiouvilleConstant(x, bounds, gamma, BC)
% function U = STURMLIOUVILLECONSTANT(x, bounds, lambda, BC, options)
%
% This function gives the analytical expression of the
% the exact solution of the Sturm Liouville equation with
% constant coefficients and without source term
%
% - u'' + (gamma*gamma) * u = 0,  a < x < b,
%
% where a and b can be infinite. If one of the bounds is
% finite, then the ODE is combined with boundary conditions
% expressed under the form
%
%   A*U + D*dU = B, with U = [u(a); u(b)] and dU = [u'(a); u'(b)]
%
% A and D are a NxN matrices and B a N-sized vector, where
% N is the number of finite bounds (N = 0, 1 or 2).
%
% INPUT: x (vector) contains the points in which the solution
%                   will be computed,
%        bounds (2x1 vector), contains the bounds of the domain;
%                             these bounds can be infinite,
%        gamma (complex), the square-root of 0-order term
%                         constant coefficient,
%        BC is an optional structure which describes the boundary
%           condition. It contains 3 fields:
%              + A, a NxN matrix corresponding to the trace
%              + D, a NxN matrix corresponding to the normal trace
%              + b, a N-sized vector corresponding to the rhs,
%           where N is the number of finite bounds (N = 0, 1, 2)
%
% OUTPUT: U the exact solution of the Sturm-Liouville boundary value
%           problem with constant coefficients and no source,

% Make sure first that the boundary conditions are coherent with
% the number of finite bounds
if (nargin < 4)

  % No boundary condition needs to be specified if the problem is
  % solved on the whole real axis
  if (~isinf(bounds(1)) || ~isinf(bounds(2)))

    error(['Aucune condition aux limites n''a ete ', ...
           'specifiee alors qu''il y a des bords.'])
  end

end

% The first bound should be smaller than the second one
bounds = sort(bounds);

% The solutions of the Sturm Liouville ODE are linear
% combinations of exp(gamma * x) and exp(-gamma * x)
if (~isinf(bounds(1)) && ~isinf(bounds(2)))

  % Case where both the bounds are finite
  % Make sure that there are two independent condtions
  if (rank([BC.A, BC.D]) < 2)

    % There are less than two valid boundary conditions
    error('Il n''y a pas assez de conditions aux limites.');

  end

  % Compute the coefficients of the linear combination
  tildeA = BC.A * [exp(gamma*bounds(:)),  exp(-gamma*bounds(:))];
  tildeD = BC.D * [exp(gamma*bounds(:)), -exp(-gamma*bounds(:))];

  C = (tildeA + gamma * tildeD) \ BC.b;

  % Deduce the expression of the solution
  U = C(1) * exp(gamma * x) + C(2) * exp(-gamma * x);

elseif (~isinf(bounds(1)) && isinf(bounds(2)))

  % Case where the upper bound is infinite
  warning('Attention a verifier ce resultat');
  U = BC.b * exp(gamma * (bounds(1) - x)) / (BC.A - gamma * BC.D);

elseif (~isinf(bounds(1)) && isinf(bounds(2)))

  % Case where the lower bound is infinite
  warning('Attention a verifier ce resultat');
  U = BC.b * exp(gamma * (x - bounds(2))) / (BC.A + gamma * BC.D);

else

  % Both bounds are infinite
  U = zeros(size(x));

end
