function [eigvals, eigvecs, AA, dEigvecs] = SchroedingerCellEVP(msh, funP, funQ, BC, numEigs, options)
% function [eigvals, eigvecs, AA] = SchroedingerCellEVP(msh, funP, funQ, BC, numEigs, options)
%
% This function computes a Lagrange P1 Finite Elements approximation
% of solutions of the eigenvalue equation
%
% - (funP(x) * u')' + funQ(x) * u = E * u,  a < x < b,
%
% combined with a boundary condition expressed in generic form
%
%   A*U + D*dU = 0, with U = [u(a); u(b)] and dU = [u'(a); u'(b)]
%
% where A and D are a 2x2 matrices.
%
% INPUT: msh (meshObject) describes the mesh,
%        funP (function), the 2-order term coefficient,
%        funQ (function), the 0-order term coefficient,
%        BC is a structure which describes the boundary condition.
%           BC contains 3 fields:
%              + A, a 2x2 matrix corresponding to the trace
%              + D, a 2x2 matrix corresponding to the normal trace
%        options is a structure with fields
%           + verbose, if true, provides additional information
%                      about the process (set to false by default)
%
% OUTPUT: U the solution of the Sturm-Liouville boundary value problem
%         dU a weak approximation of the derivative of U
%         AA, the finite elements matrix
%
% NOTE: For essential boundary conditions, we use elimination


if (nargin < 6)

  options.verbose = 0;

end

% Extract some mesh information
N = msh.numPoints;    % Number of nodes

% Compute the FE matrices
[MM, KK] = FEmatrices(msh, funQ, funP);
[MM_Id, ~, DD] = FEmatrices(msh);

% =================== %
% Boundary conditions %
% =================== %
if (rank([BC.A, BC.D]) < 2)

  % There are less than 2 valid boundary conditions due to redudancy
  error('Les conditions aux limites ne sont pas admissibles');

end

if (rank(BC.D) == 2)

  % Case where D is invertible
  % u'(a) and u'(b) therefore can be expressed as affine
  % combinations of u(a) and u(b).
  % This corresponds to generic Robin boundary conditions
  if (options.verbose)

    fprintf('Conditions de Robin reconnues\n');

  end

  % Compute boundary matrices
  S11 = sparse(N, N); S11(msh.boundsIds(1), msh.boundsIds(1)) = 1;
  S12 = sparse(N, N); S12(msh.boundsIds(1), msh.boundsIds(2)) = 1;
  S21 = sparse(N, N); S21(msh.boundsIds(2), msh.boundsIds(1)) = 1;
  S22 = sparse(N, N); S22(msh.boundsIds(2), msh.boundsIds(2)) = 1;

  tildeA = BC.D \ BC.A;

  % The FE matrix with the boundary contribution
  SS = - S11 * tildeA(1, 1) - S12 * tildeA(1, 2) +...
         S21 * tildeA(2, 1) + S22 * tildeA(2, 2);

  AA = MM + KK + SS;
  BB = MM_Id;

  % Eigenpairs
  [eigvecs, eigvals] = eigs(AA, BB, numEigs, 'smallestabs');
  eigvals = diag(eigvals).';

end

if (rank(BC.D) == 0)

  % Case where D = 0
  % Since [A D] is a 2-rank matrix, A is invertible
  % The boundary condition involves only u(a) and u(b)
  % This corresponds to Dirichlet boundary conditions
  % We introduce a lift and we use elimination
  if (options.verbose)

    fprintf('Conditions de Dirichlet reconnues\n');

  end

  % The Dirichlet projection matrix
  PP = sparse( (1:N-2)', (2:N-1)', ones(N-2, 1), N-2, N );

  % The FE matrices
  AA = MM + KK;
  BB = MM_Id;
  AA0 = PP * AA * PP';
  BB0 = PP * BB * PP';

  % Eigenpairs
  [eigvecs0, eigvals] = eigs(AA0, BB0, numEigs, 'smallestabs');
  eigvals = diag(eigvals).';
  eigvecs = PP' * eigvecs0;

end

if (rank(BC.D) == 1)

  % This a case encompasses variations of (quasi-)periodic conditions
  % mixed condition, or many exotic other boundary conditions which
  % cannot be treated using P1 Lagrange finite elements
  %
  % The idea is to write D = d1 * d2.', where d1 and d2 are 2x1 vectors.
  % From there, the boundary conditons can be reorganized as
  % d2.' * dU = b1 - a1.' * U (a natural condition)
  %  c.' *  U = b2            (an essential condition)
  %
  % A compatibility condition needs to be checked if we do not want
  % the normal traces of u to figure in the formulation.
  if (options.verbose)

    fprintf('Conditions (quasi-)periodiques, mixtes, ou exotiques\n');

  end

  % Write D as a product of 2x1 vectors d1 and d2
  % Find a non null element in D. Here we choose the maximum
  [~, maxId] = max(abs(BC.D(:)));

  % Find the indices of the max (Im, Km) and the other indices (I0, K0)
  Km = 1 + floor((maxId(1)-1)/2);  % d2(Km) will be non null
  Im = maxId(1) - 2*(Km - 1);      % d1(Im) will be non null

  I0 = mod(Im, 2) + 1; % If Im = 1 then I0 = 2, and vice versa
  K0 = mod(Km, 2) + 1; % If Km = 1 then K0 = 2, and vice versa

  % Compute d1 and d2 such that D = d1 * d2.'
  d1 = zeros(2, 1);
  d1(Im) = 1;
  d1(I0) = BC.D(I0, Km) / BC.D(Im, Km);

  d2 = zeros(2, 1);
  d2(Km) = BC.D(Im, Km);
  d2(K0) = BC.D(Im, K0);

  % One can deduce the essential condition
  % and express it under the form A0 * u = 0
  A0 = BC.A(I0, :) - BC.A(Im, :) * d1(I0) / d1(Im);

  % If we do not want the normal traces to appear
  % in the variational formulation, then the test
  % function also need to satisfy a given essential
  % condition.
  % Check if the essential conditions satisfied by
  % the solution and the test function are the same
  % or proportional
  if (rank([A0(1), A0(2); conj(d2(2)), conj(d2(1))]) ~= 1)

    % The case with non-equivalent essential conditions
    % is not handled
    error(sprintf(['Conditions de type Cauchy reconnues.\n',...
           'Si on veut faire disparaitre les traces ', ...
           'normales de la formulation, les fonctions ', ...
           'test doivent verifier une condition ', ...
           'essentielle differente de celle de ', ...
           'la solution.\nDe telles conditions ne sont ', ...
           'pas prises en charge.']));  % #ok

  else

    % Projection matrix
    PP = sparse((2:N-1)', (2:N-1)', ones(N-2, 1), N-1, N);
    PP(1, msh.boundsIds(Km)) =  d2(Km);
    PP(1, msh.boundsIds(K0)) = -d2(K0);

    % Boundary matrices
    Sm1 = sparse(N, N); Sm1(msh.boundsIds(Km), msh.boundsIds(1)) = 1;
    Sm2 = sparse(N, N); Sm2(msh.boundsIds(Km), msh.boundsIds(2)) = 1;

    SS = -((-1)^K0) * (BC.A(Im, 1) * Sm1 + BC.A(Im, 2) * Sm2) /...
                                            (d1(Im) * d2(Km));
    
    % FE matrices
    AA = MM + KK + SS;
    BB = MM_Id;
    AA0 = PP * AA * PP';
    BB0 = PP * BB * PP';

    % Eigenpairs
    [eigvecs0, eigvals] = eigs(AA0, BB0, numEigs, 'smallestabs');
    eigvals = diag(eigvals).';
    eigvecs = PP' * eigvecs0;

  end

end

dEigvecs = MM_Id \ (DD * eigvecs);