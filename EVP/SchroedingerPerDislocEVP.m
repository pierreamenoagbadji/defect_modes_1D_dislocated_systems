% *************************************************** %
% SCHROEDINGERPERDISLOC.m 
% Computes defect modes for the Schroedinger operator
%
% -(μ(x)*u')'(x) + V(x) u(x) acting on L²(R),
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

%  ********** %
%% Parameters %
%  ********** %
pbInputs.IDb = [-0.5 0.5]; % Bounds of the interior domain
pbInputs.pos.numCells = 8; % Number of cells on which to compute solution on Ω_+ 
pbInputs.neg.numCells = 8; % Number of cells on which to compute solution on Ω_-

% Coefficients outside the interior domain
% The coefficients are periodic outside the
% interior domain
kappa = @(X) domainwall(X+0.5);
Veven = @(x) 10*cos(4*pi*x);
Vodd  = @(x) 30*cos(2*pi*x);

pbInputs.pos.per = 1;
pbInputs.pos.mu  = @(x) ones(size(x, 1), 1);
pbInputs.pos.V   = @(x) Veven(x) + kappa(x) .* Vodd(x);

pbInputs.neg.per = 1;
pbInputs.neg.mu  = @(x) ones(size(x, 1), 1);
pbInputs.neg.V   = @(x) Veven(x) + kappa(x) .* Vodd(x);

% Coefficients inside the interior domain
pbInputs.int.mu = @(x) ones(size(x, 1), 1);
pbInputs.int.V  = @(x) Veven(x) + kappa(x) .* Vodd(x);

% Gap range
pbInputs.specRange = [pi^2 + 5, pi^2 + 10];
pbInputs.numSpec = 1024;
pbInputs.numEigs = 4; % Number of eigenvalues for interior problem

% Mesh nodes
pbInputs.numNodes = 256;

%  ************ %
%% Run the code %
%  ************ %
EVPlocPert(pbInputs);

function EVPlocPert(pbInputs)

  %  ********** %
  %% Parameters %
  %  ********** %
  IDb = pbInputs.IDb;

  numCellsPos = pbInputs.pos.numCells;
  numCellsNeg = pbInputs.neg.numCells;

  perPos = pbInputs.pos.per;
  muPos  = pbInputs.pos.mu;
  Vpos   = pbInputs.pos.V; 

  perNeg = pbInputs.neg.per;
  muNeg  = pbInputs.neg.mu;
  Vneg   = pbInputs.neg.V; 

  % Make sure the periods are positive
  if ~(perPos > 0 & perNeg > 0)
    error(['The period of coefficients on Ω_+ (', num2str(perPos, '%0.3e'), ') and on Ω_- (', num2str(perNeg, '%0.3e'), ') should be positive.'])
  end

  % Coefficients inside the interior domain
  muInt = pbInputs.int.mu;
  Vint  = pbInputs.int.V;

  % Range for searching the eigenvalue
  specRange = pbInputs.specRange;
  numSpec   = pbInputs.numSpec;
  specvec   = linspace(specRange(1), specRange(2), numSpec)';

  % Number of eigenvalues to compute
  numEigs = pbInputs.numEigs;

  %% **** %
  %  Mesh %
  %% **** %
  numNodesPos = ceil(pbInputs.numNodes * max(perPos, 1));
  numNodesNeg = ceil(pbInputs.numNodes * max(perNeg, 1));
  numNodesInt = ceil(pbInputs.numNodes * max(abs(IDb(2) - IDb(1)), 1));

  % Meshes
  mshPos = meshObject('uniform', 0,  perPos, numNodesPos);
  mshNeg = meshObject('uniform', 0, -perNeg, numNodesNeg);
  mshInt = meshObject('uniform', IDb(1), IDb(2), numNodesInt);

  %% ***************** %
  %  Plot coefficients %
  %% ***************** %
  figCoeffs = figure;

  for idI = 1:numCellsPos
    X = IDb(2) + mshPos.points + (idI - 1.0) * perPos;
    subplot(2, 1, 1); plot(X, muPos(X), 'r'); hold on;
    subplot(2, 1, 2); plot(X,  Vpos(X), 'b'); hold on;
  end

  for idI = 1:numCellsNeg
    X = IDb(1) + mshNeg.points - (idI - 1.0) * perNeg;
    subplot(2, 1, 1); plot(X, muNeg(X), 'r');
    subplot(2, 1, 2); plot(X,  Vneg(X), 'b');
  end

  X = mshInt.points;
  xmin = IDb(1) + mshNeg.bounds(2) - (numCellsNeg - 1.0) * perNeg;
  xmax = IDb(2) + mshPos.bounds(2) + (numCellsPos - 1.0) * perPos;

  subplot(2, 1, 1); 
  plot(X, muInt(X), 'r');
  title('$x \mapsto \mu (x)$');
  set(gca, 'FontSize', 16);
  xlim([xmin, xmax]);

  subplot(2, 1, 2); 
  plot(X, Vint(X), 'b');
  title('$x \mapsto V (x)$');
  set(gca, 'FontSize', 16);
  xlim([xmin, xmax]);

  %  **************************** %
  %% Eigenvalue searching process %
  %  **************************** %
  % Initialize eigenvalues and eigenvectors
  eigvals = zeros(numSpec, numEigs);
  eigvecs = zeros(mshInt.numPoints, numEigs, numSpec);

  RobinNegVec = 1i * sqrt(abs(specvec));
  RobinPosVec = 1i * sqrt(abs(specvec));

  BCpos = struct('A', [], 'D', -1, 'b', 1);
  BCneg = struct('A', [], 'D',  1, 'b', 1);
  BCint = struct('A', zeros(2), 'D', diag([-1 1]));

  muPosTrans = @(x) muPos(x + IDb(2)); VposTrans  = @(x)  Vpos(x + IDb(2));
  muNegTrans = @(x) muNeg(x + IDb(1)); VnegTrans  = @(x)  Vneg(x + IDb(1));
  F = @(x) zeros(size(x, 1), 1);

  for idI = 1:numSpec

    specVar  = specvec(idI);
    RobinNeg = RobinNegVec(idI);
    RobinPos = RobinPosVec(idI);
    
    % Step 1
    % ****** %
    % Compute RtR coefficients obtained by solving the
    % auxiliary half-line problems
    %
    % - (μ_± (x) * u'_±)' + (V_± (x) - E) u_± = 0 on Ω_±,
    %                      - (u_±)' + r_± u_± = 1 at x = a_±
    % Plus side
    BCpos.A = RobinPos;
    [~, RtRpos] = PeriodicBVP(mshPos, muPosTrans, VposTrans, F,...
      BCpos, numCellsPos, struct('compute_RtR', true));

    % Minus side
    BCneg.A = RobinNeg;
    [~, RtRneg] = PeriodicBVP(mshNeg, muNegTrans, VnegTrans, F,...
      BCneg, numCellsNeg, struct('compute_RtR', true));

    % Step 2
    % ****** %
    % Construct boundary condition for eigenvalue
    % problem in interior domain
    BCint.A(1, 1) = RobinNeg * (RtRneg - 1) / (RtRneg + 1);
    BCint.A(2, 2) = RobinPos * (RtRpos - 1) / (RtRpos + 1);

    % Solve eigenvalue problem
    [eigvals(idI, :), eigvecs(:, :, idI)] = SchroedingerCellEVP(mshInt, muInt, Vint, BCint, numEigs);

    % Print progress
    fprintf('%3d%%\n', round(100*idI/numSpec));

  end

  % ********************************************* %
  % Plot the difference in absolute value between %
  % the eigenvalues and the spectral parameter    %
  % ********************************************* %
  figEigvals = figure;
  hold off;
  idEig  = 3;
  subplot(1, 2, 1);
  specIndic = abs(eigvals(:, idEig) - specvec);
  plot(specvec, specIndic, 'b');
  xlabel('$E$')
  title(['$E \mapsto |\mu_', int2str(idEig), ' (E) - E|$']);
  set(gca, 'FontSize', 16);
  hold on;

  subplot(1, 2, 2);
  vert = [035, 120, 057];
  cols = [255, 0, 0; vert; 0, 0, 255; 213, 52, 235] / 255;
  for idE = 1:numEigs
    X = specvec;
    Y = real(eigvals(:, idE));
    
    % Remove sharp jumps coming from
    % how eigenvalues are labelled
    idX = [false; abs(diff(Y)) > 0.9];
    X(idX) = NaN;
    Y(idX) = NaN;
    
    plot(X, Y, 'b', 'Color', cols(idE, :)); hold on;
  end
  plot(specvec, specvec, 'r--');
  xlabel('$E$');
  title('$E \mapsto \Re \mu (E)$');
  set(gca, 'FontSize', 16);

  error('test');
  %% Uncomment for validation for constant interior potential
  %
  % % ************************************************** %
  % % The exact eigenvalue satisfies                     %
  % %     sqrt(-E) = sqrt(C + E) * tan(A * sqrt(C + E))  %
  % % provided that Vint is constant, equal to -C,       %
  % % and IDb = [-A, A].                                 %
  % % ************************************************** %
  % C = -Vint(0);
  % A = IDb(2);
  % dispfun = @(E) sqrt(-E) - sqrt(C + E) * tan(A * sqrt(C + E));
  % E0 = fzero(dispfun, [-C, 0]);

  % figure;
  % subplot(1, 2, 1);
  % plot([E0 E0], [min(specIndic) max(specIndic)], '--r');

  % subplot(1, 2, 2);
  % plot(xneg, cos(A * sqrt(C + E0)) * exp(-sqrt(-E0) * (abs(xneg) - A)), 'b--'); hold on;
  % plot(xpos, cos(A * sqrt(C + E0)) * exp(-sqrt(-E0) * (abs(xpos) - A)), 'b--');
  % plot(msh.points, cos(sqrt(C + E0) * msh.points), 'b--');

  %%
  % ********************* %
  % Numerical defect mode %
  % ********************* %
  tol = 1e-3;
  defectId = find(abs(eigvals(:, idEig) - specvec) < tol);

  if ~isempty(defectId)
    
    numModes = length(defectId);

    % Print the eigenvalues
    fprintf('%d defect mode(s) found: ', numModes);
    defect.val = specvec(defectId);
    fprintf(repmat('%5e\t', [1 3]), specvec(defectId));
    fprintf('\n');

    % Compute eigenfunctions (defect modes)
    for idJ = 1:numModes

      d = defectId(idJ);
      specVar  = specvec(d);

      defect.modeInt(:, idJ) = eigvecs(:, idEig, d);
      
      RobinNeg = RobinNegVec(d);
      RobinPos = RobinPosVec(d);
      
      coeffNeg = -2*defect.modeInt(msh.boundsIds(1)) * RobinNegVec(d) / (RtRneg(specVar, RobinNeg) + 1);
      coeffPos = -2*defect.modeInt(msh.boundsIds(2)) * RobinPosVec(d) / (RtRpos(specVar, RobinPos) + 1);

      defect.modeNeg(:, idJ) = coeffNeg * uNeg(xneg, specvec(d), RobinNeg);
      defect.modePos(:, idJ) = coeffPos * uPos(xpos, specvec(d), RobinPos);

    end

    %% Plot defect modes
    figure(2);
    subplot(2, 1, 1);
    plot(xneg, Vneg, 'color', vert); hold on;
    plot(msh.points, Vint(msh.points), 'color', vert);
    plot(xpos, Vpos, 'color', vert);
    xline(IDb(1));
    xline(IDb(2));
    xlabel('$x$');
    title('Potential');
    set(gca, 'FontSize', 16);

    subplot(2, 1, 2);
    cst = max(abs(defect.modeInt(:, 1)));
    plot(xneg,       real(defect.modeNeg(:, 1)/cst), 'r'); hold on;
    plot(msh.points, real(defect.modeInt(:, 1)/cst), 'r');
    plot(xpos,       real(defect.modePos(:, 1)/cst), 'r');
    xline(IDb(1));
    xline(IDb(2));
    xlabel('$x$');
    title('Defect mode');
    set(gca, 'FontSize', 16);

  else

    warning('No defect mode found. Use a finer mesh, or change the tolerance.')

  end

end
