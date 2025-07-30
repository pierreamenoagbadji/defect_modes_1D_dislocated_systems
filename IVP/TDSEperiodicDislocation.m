% *************************************************** %
% TDSEPeriodicDislocation.m 
% Solve the time-dependent Schroedinger equation (TDSE)
%
% i(∂Ψ/∂t) (x, t) = -∂(μ(x)*(∂Ψ/∂x) (x, t))/∂x + V(x) Ψ(x, t) + f(x, t) on R,
%
%       Ψ(x, t=0) = Ψ0(x),
%
% where there exist a_+ > a_- such that for
% Ω_- := (-∞, a_-), 
% Ω_0 := (a_-, a_+), and
% Ω_+ := (a_+, +∞), 
% the coefficient μ(x) (resp. V(x)) coincides in Ω_± 
% with a T_±-periodic function μ_±(x) (resp. V_±(x))
% and with anything bounded in Ω_0.
% *************************************************** %

%  ************** %
%% Initialization %
%  ************** %
clear; clc;
close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
vert = [035, 120, 057] / 255;

%  ********** %
%% Parameters %
%  ********** %
IDb = [-1.0 1.0]; % Bounds of the interior domain
numCellsPos = 3; % Number of cells on which to compute solution on Ω_+ 
numCellsNeg = 3; % Number of cells on which to compute solution on Ω_-

% Coefficients outside the interior domain
% The coefficients are periodic outside the
% interior domain
kappa = @(X) -1 + 2*domainwall(X + 0.5);
Veven = @(x) 10*cos(4*pi*x);
Vodd  = @(x) 30*cos(2*pi*x);

perPos = 1;
muPos = @(x) ones(size(x, 1), 1);
Vpos  = @(x) zeros(size(x, 1), 1); % cos(2*pi*x); % Veven(x) + kappa(x) .* Vodd(x);


perNeg = 1;
muNeg  = @(x) ones(size(x, 1), 1);
Vneg   = @(x) zeros(size(x, 1), 1); % cos(2*pi*x); % Veven(x) + kappa(x) .* Vodd(x);

% Make sure the periods are positive
if ~(perPos > 0 && perNeg > 0)
  error(['The period of coefficients on Ω_+ (', num2str(perPos, '%0.3e'),...
          ') and on Ω_- (', num2str(perNeg, '%0.3e'), ') should be positive.']);
end

% Coefficients inside the interior domain
muInt = @(x) ones(size(x, 1), 1);
Vint  = @(x) zeros(size(x, 1), 1); % cos(2*pi*x); % Veven(x) + kappa(x) .* Vodd(x);

% Source term and initial data
% are compactly supported in interior domain
rhsInt = @(x, t) zeros(size(x, 1), 1);
asol0 = 30;
sol0   = @(x) exp(-asol0 * x.^2);% cutoff(x, -0.5, 0.5);

%  ************************* %
%% Discretization parameters %
%  ************************* %
numNodes = 256;
delta_t = 1e-3;
final_t = 1;
theta = 0.5;

%  **** %
%% Mesh %
%  **** %
numNodesPos = ceil(numNodes * max(perPos, 1));
numNodesNeg = ceil(numNodes * max(perNeg, 1));
numNodesInt = ceil(numNodes * max(abs(IDb(2) - IDb(1)), 1));

% Meshes
mshPos = meshObject('uniform', 0,  perPos, numNodesPos);
mshNeg = meshObject('uniform', 0, -perNeg, numNodesNeg);
mshInt = meshObject('uniform', IDb(1), IDb(2), numNodesInt);

Npos = mshPos.numPoints;
Nneg = mshNeg.numPoints;
Nint = mshInt.numPoints;

%  ***************** %
%% Plot coefficients %
%  ***************** %
figCoeffs = figure;

for idI = 1:numCellsPos
  X = IDb(2) + mshPos.points + (idI - 1.0) * perPos;
  subplot(3, 1, 1); plot(X, muPos(X), 'r'); hold on;
  subplot(3, 1, 2); plot(X,  Vpos(X), 'b'); hold on;
end

for idI = 1:numCellsNeg
  X = IDb(1) + mshNeg.points - (idI - 1.0) * perNeg;
  subplot(3, 1, 1); plot(X, muNeg(X), 'r');
  subplot(3, 1, 2); plot(X,  Vneg(X), 'b');
end

X = mshInt.points;
xmin = IDb(1) + mshNeg.bounds(2) - (numCellsNeg - 1.0) * perNeg;
xmax = IDb(2) + mshPos.bounds(2) + (numCellsPos - 1.0) * perPos;

subplot(3, 1, 1); 
plot(X, muInt(X), 'r');
title('$x \mapsto \mu (x)$');
set(gca, 'FontSize', 16);
xlim([xmin, xmax]);

subplot(3, 1, 2); 
plot(X, Vint(X), 'b');
title('$x \mapsto V (x)$');
set(gca, 'FontSize', 16);
xlim([xmin, xmax]);

subplot(3, 1, 3); 
plot(X, rhsInt(X, 0), 'Color', vert);
title('$x \mapsto f (x)$');
set(gca, 'FontSize', 16);
xlim([xmin, xmax]);

%  ****************** %
%% Initialize outputs %
%  ****************** %
numTsteps = ceil(final_t / delta_t);
out.pos.prop = zeros(numTsteps, 1);
out.pos.DtN  = zeros(numTsteps, 1);
auxSolPos    = zeros(Npos, numCellsPos, numTsteps);
out.pos.sol  = zeros(Npos, numCellsPos, numTsteps);

out.neg.prop = zeros(numTsteps, 1);
out.neg.DtN  = zeros(numTsteps, 1);
auxSolNeg    = zeros(Nneg, numCellsNeg, numTsteps);
out.neg.sol  = zeros(Nneg, numCellsNeg, numTsteps);

out.int.sol = zeros(Nint, numTsteps);

%  ****************************** %
%% First time step: positive side %
%  ****************************** %
muPosTrans = @(x) muPos(x + IDb(2));
VposTrans  = @(x)  Vpos(x + IDb(2));

% Local cell solutions %
% ******************** %
% FE matrices
[MM, KK] = FEmatrices(mshPos, VposTrans, muPosTrans);
MM0 = FEmatrices(mshPos);
out.pos.PP = sparse((1:Npos-2)', (2:Npos-1)', ones(Npos-2, 1), Npos-2, Npos);
% ********************************************************* %
out.pos.AA =       theta  * delta_t * (MM + KK) - 1i * MM0; %
out.pos.BB = (-1 + theta) * delta_t * (MM + KK) - 1i * MM0; %
% ********************************************************* %
AA0 = out.pos.PP * out.pos.AA * out.pos.PP';

% Right-hand side
Ub = sparse(Npos, 2);
Ub(mshPos.boundsIds, :) = eye(2);
LL = -out.pos.AA * Ub;
LL0 = out.pos.PP * LL;

% Solve the linear system and deduce the solution
eEliminated = AA0 \ LL0;
e = Ub + out.pos.PP' * eEliminated;

out.pos.e0 = e(:, 1);
out.pos.e1 = e(:, 2);

% DtN coefficients %
% **************** %
out.pos.t00 = full(out.pos.e0' * out.pos.AA * out.pos.e0);
out.pos.t10 = full(out.pos.e0' * out.pos.AA * out.pos.e1);
out.pos.t01 = full(out.pos.e1' * out.pos.AA * out.pos.e0);
out.pos.t11 = full(out.pos.e1' * out.pos.AA * out.pos.e1);

% Propagation coefficient                       %
% The unique solution of the quadratic equation %
% t10 * p^2 + (t00 + t11) * p + t01 = 0         %
% with a strictly less than 1 absolute value.   %
% ********************************************* %
disc = sqrt((out.pos.t00 + out.pos.t11)^2 - 4*out.pos.t10 * out.pos.t01);
R = min([(-(out.pos.t00 + out.pos.t11) - disc) / (2 * out.pos.t10),...
         (-(out.pos.t00 + out.pos.t11) + disc) / (2 * out.pos.t10)],...
        [], 'ComparisonMethod', 'abs');

% R should be strictly less than one, as we are
% in the resolvent set.
if (abs(R) >= 1.0 - eps)
  error('The propagation coefficient should be less than 1 in absolute value: %0.5e', abs(R));
end

% Compute auxiliary solution at first time step %
% ********************************************* %
auxSolPos(:, :, 1) = (out.pos.e0 + out.pos.e1 * R) * R.^(0:numCellsPos-1);

% Update solution %
% *************** %
out.pos.prop(1) = R;
out.pos.DtN(1)  = out.pos.t10 * R + out.pos.t00;

%  ****************************** %
%% First time step: negative side %
%  ****************************** %
muNegTrans = @(x) muNeg(x + IDb(1));
VnegTrans  = @(x)  Vneg(x + IDb(1));

% Local cell solutions %
% ******************** %
% FE matrices
[MM, KK] = FEmatrices(mshNeg, VnegTrans, muNegTrans);
MM0 = FEmatrices(mshNeg);
out.neg.PP = sparse((1:Nneg-2)', (2:Nneg-1)', ones(Nneg-2, 1), Nneg-2, Nneg);
% ********************************************************* %
out.neg.AA =       theta  * delta_t * (MM + KK) - 1i * MM0; %
out.neg.BB = (-1 + theta) * delta_t * (MM + KK) - 1i * MM0; %
% ********************************************************* %
AA0 = out.neg.PP * out.neg.AA * out.neg.PP';

% Right-hand side
Ub = sparse(Nneg, 2);
Ub(mshNeg.boundsIds, :) = eye(2);
LL = -out.neg.AA * Ub;
LL0 = out.neg.PP * LL;

% Solve the linear system and deduce the solution
eEliminated = AA0 \ LL0;
e = Ub + out.neg.PP' * eEliminated;

out.neg.e0 = e(:, 1);
out.neg.e1 = e(:, 2);

% DtN coefficients                                 %
% Le signe '-' est nécessaire, contrairement à la  % 
% théorie, car le produit (ej' * AA * el) est égal %
% à une intégrale qui ***part de 0 à -perNeg***, à %
% cause de la façon dont mon maillage est défini.  %
% ************************************************ %
out.neg.t00 = -full(out.neg.e0' * out.neg.AA * out.neg.e0);
out.neg.t10 = -full(out.neg.e0' * out.neg.AA * out.neg.e1);
out.neg.t01 = -full(out.neg.e1' * out.neg.AA * out.neg.e0);
out.neg.t11 = -full(out.neg.e1' * out.neg.AA * out.neg.e1);

% Propagation coefficient                       %
% The unique solution of the quadratic equation %
% t10 * p^2 + (t00 + t11) * p + t01 = 0         %
% with a strictly less than 1 absolute value.   %
% ********************************************* %
disc = sqrt((out.neg.t00 + out.neg.t11)^2 - 4*out.neg.t10 * out.neg.t01);
R = min([(-(out.neg.t00 + out.neg.t11) - disc) / (2 * out.neg.t10),...
         (-(out.neg.t00 + out.neg.t11) + disc) / (2 * out.neg.t10)],...
        [], 'ComparisonMethod', 'abs');

% R should be strictly less than one, as we are
% in the resolvent set.
if (abs(R) >= 1.0 - eps)
  error('The propagation coefficient should be less than 1 in absolute value: %0.5e', abs(R));
end

% Compute auxiliary solution at first time step %
% ********************************************* %
auxSolNeg(:, :, 1) = (out.neg.e0 + out.neg.e1 * R) * R.^(0:numCellsNeg-1);

% Update solution %
% *************** %
out.neg.prop(1) = R;
out.neg.DtN(1)  = out.neg.t10 * R + out.neg.t00;

% solrefNeg = exp(sqrt(-1i / (delta_t * theta))*mshNeg.points);
% plot(mshNeg.points, imag(solrefNeg)); hold on;
% plot(mshNeg.points, imag(auxSolNeg{1}(:, 1)));
% fprintf('%0.5d\n', sqrt(((auxSolNeg{1}(:, 1) - solrefNeg)' * MM0 * (auxSolNeg{1}(:, 1) - solrefNeg)) / (solrefNeg' * MM0 * solrefNeg)));

%  ********************************* %
%% First time step: interior problem %
%  and complete reconstruction.      %
%  ********************************* %
[MM, KK] = FEmatrices(mshInt, Vint, muInt);
out.int.MM0 = FEmatrices(mshInt);
% ***************************************************************** %
out.int.AA =       theta  * delta_t * (MM + KK) - 1i * out.int.MM0; %
out.int.BB = (-1 + theta) * delta_t * (MM + KK) - 1i * out.int.MM0; %
% ***************************************************************** %

% Surface contributions
S11 = sparse(Nint, Nint); 
S22 = sparse(Nint, Nint); 
S11(mshInt.boundsIds(1), mshInt.boundsIds(1)) = 1;
S22(mshInt.boundsIds(2), mshInt.boundsIds(2)) = 1;

out.int.AA = out.int.AA + out.neg.DtN(1) * S11 + out.pos.DtN(1) * S22;

% Right-hand side
rhsTheta = @(idT)       theta  * rhsInt(mshInt.points, (idT  ) * delta_t) +...
                  (-1 + theta) * rhsInt(mshInt.points, (idT-1) * delta_t);
LLint = out.int.BB * sol0(mshInt.points) - out.int.MM0 * rhsTheta(1);

% Interior solution % 
% ***************** %
out.int.sol(:, 1) = out.int.AA \ LLint;

% Construct entire solution
out.pos.sol(:, :, 1) = out.int.sol(mshInt.boundsIds(2), 1) * auxSolPos(:, :, 1);
out.neg.sol(:, :, 1) = out.int.sol(mshInt.boundsIds(1), 1) * auxSolNeg(:, :, 1);

figSol = figure;

for idI = 1:numCellsPos
  X = IDb(2) + mshPos.points + (idI - 1.0) * perPos;
  plot(X, real(out.pos.sol(:, idI, 1)), 'b'); hold on;
end

for idI = 1:numCellsNeg
  X = IDb(1) + mshNeg.points - (idI - 1.0) * perNeg;
  plot(X, real(out.neg.sol(:, idI, 1)), 'b'); hold on;
end

plot(mshInt.points, real(out.int.sol(:, 1)), 'b');
xmin = IDb(1) + mshNeg.bounds(2) - (numCellsNeg - 1.0) * perNeg;
xmax = IDb(2) + mshPos.bounds(2) + (numCellsPos - 1.0) * perPos;
title(['$t_1 = ', num2str(delta_t, '%0.5e'), '$'])
axis([xmin, xmax, -0.5, 1.0]);

% Reference solution
% if initial data is a gaussian
Uref = @(x, t) (1 / sqrt(1 + 4i * asol0*t)) * exp(-asol0 * x.^2 / (1 + 4i * asol0*t));

x = linspace(xmin, xmax, 2^10);
plot(x, real(Uref(x, delta_t)), 'r*', 'MarkerSize', 1);
% error;
pause;

%  *************** %
%% Next time steps %
%  *************** %
for idT = 2:numTsteps

  %  ********************************* %
  %% Auxiliary solution: positive side %
  %  ********************************* %
  % Local cell solution with right-hand side
  AA0 = out.pos.PP * out.pos.AA * out.pos.PP';
  LL  = out.pos.BB * auxSolPos(:, 1, idT-1);
  LL0 = out.pos.PP * LL;

  eF  = out.pos.PP' * (AA0 \ LL0);

  % Source-to-Neumann coefficients
  G0pos = full(out.pos.e0' * (out.pos.AA * eF - LL));
  G1pos = full(out.pos.e1' * (out.pos.AA * eF - LL));

  % Compute the current propagation coefficient
  out.pos.prop(idT) = -(G1pos + out.pos.prop(1) * G0pos +...
    out.pos.prop(idT-1:-1:2).' * out.pos.DtN(2:idT-1)) /...
    (out.pos.t11 + out.pos.t10 * out.pos.prop(1) + out.pos.DtN(1));

  % Deduce DtN coefficient
  out.pos.DtN(idT) = G0pos + out.pos.t10 * out.pos.prop(idT);

  % Compute auxiliary solution cell by cell
  auxSolPos(:, 1, idT) = eF + out.pos.prop(idT) * out.pos.e1;
  for idCell = 2:numCellsPos
    auxSolPos(:, idCell, idT) = reshape(auxSolPos(:, idCell-1, 1:idT), Npos, idT)...
                              * out.pos.prop(idT:-1:1);
  end

  %  ********************************* %
  %% Auxiliary solution: negative side %
  %  ********************************* %
  % Local cell solution with right-hand side
  AA0 = out.neg.PP * out.neg.AA * out.neg.PP';
  LL  = out.neg.BB * auxSolNeg(:, 1, idT-1);
  LL0 = out.neg.PP * LL;

  eF  = out.neg.PP' * (AA0 \ LL0);

  % source-to-Neumann coefficients
  % Again, the '-' sign is important.
  G0neg = -full(out.neg.e0' * (out.neg.AA * eF - LL));
  G1neg = -full(out.neg.e1' * (out.neg.AA * eF - LL));

  % Compute the current propagation coefficient
  out.neg.prop(idT) = -(G1neg + out.neg.prop(1) * G0neg + ...
    out.neg.prop(idT-1:-1:2).' * out.neg.DtN(2:idT-1)) /...
    (out.neg.t11 + out.neg.t10 * out.neg.prop(1) + out.neg.DtN(1));

  % Deduce DtN coefficient
  out.neg.DtN(idT) = G0neg + out.neg.t10 * out.neg.prop(idT);

  % Compute auxiliary solution cell by cell
  auxSolNeg(:, 1, idT) = eF + out.neg.prop(idT) * out.neg.e1;
  for idCell = 2:numCellsNeg
    auxSolNeg(:, idCell, idT) = reshape(auxSolNeg(:, idCell-1, 1:idT), Nneg, idT)...
                              * out.neg.prop(idT:-1:1);
  end

  %  **************** %
  %% Interior problem %
  %  **************** %
  LLint = out.int.BB * out.int.sol(:, idT-1) - out.int.MM0 * rhsTheta(idT);

  % Surface contributions
  LLint(mshInt.boundsIds(1)) = LLint(mshInt.boundsIds(1)) -...
    out.int.sol(mshInt.boundsIds(1), idT-1:-1:1) * out.neg.DtN(2:idT);
  
  LLint(mshInt.boundsIds(2)) = LLint(mshInt.boundsIds(2)) -...
    out.int.sol(mshInt.boundsIds(2), idT-1:-1:1) * out.pos.DtN(2:idT);

  % Interior solution
  out.int.sol(:, idT) = out.int.AA \ LLint;

  %  ************************* %
  %% Construct entire solution %
  %  ************************* %
  for idCell = 1:numCellsPos
    out.pos.sol(:, idCell, idT) = reshape(auxSolPos(:, idCell, 1:idT), Npos, idT) *...
      out.int.sol(mshInt.boundsIds(2), idT:-1:1).';
  end
  
  for idCell = 1:numCellsNeg
    out.neg.sol(:, idCell, idT) = reshape(auxSolNeg(:, idCell, 1:idT), Nneg, idT) *...
      out.int.sol(mshInt.boundsIds(1), idT:-1:1).';
  end

  figure(figSol);
  hold off;

  for idI = 1:numCellsPos
    X = IDb(2) + mshPos.points + (idI - 1.0) * perPos;
    plot(X, real(out.pos.sol(:, idI, idT)), 'b'); hold on;
  end

  for idI = 1:numCellsNeg
    X = IDb(1) + mshNeg.points - (idI - 1.0) * perNeg;
    plot(X, real(out.neg.sol(:, idI, idT)), 'b'); hold on;
  end

  plot(mshInt.points, real(out.int.sol(:, idT)), 'b');

  % Reference solution
  plot(x, real(Uref(x, idT * delta_t)), 'r*', 'MarkerSize', 1);
  
  axis([xmin, xmax, -0.5, 1.0]);
  title(['$t_{', int2str(idT), '} = ', num2str(idT*delta_t, '%0.5e'), '$']);
  pause(0.05);

end