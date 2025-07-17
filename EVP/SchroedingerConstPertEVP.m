% *************************************************** %
% SCHROEDINGERCONSTPERT.m 
% Computes defect modes for the Schroedinger operator
%
% -u''(x) + V(x) u(x) acting on L²(R),
%
% where V(x) is constant outside a domain Ω_0
% and anything bounded inside.
% *************************************************** %

% ************** %
% Initialization %
% ************** %
clear; clc;
close all;
vert = [035, 120, 057] / 255;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% ********** %
% Parameters %
% ********** %
% Bounds of the interior domain
IDb = [-1.5 1.5];
ODb = [-10 10];

% Potential inside the interior domain
% Vint = @(x) -10.0*ones(size(x, 1), 1);
Vint = @(x) -8.0*cutoff(x, -1-1/4, -1+1/4, 0.1)...
            -8.0*cutoff(x, +1-1/4, +1+1/4, 0.1);
funP = @(x)      ones(size(x, 1), 1);

% Potential on each side of the interior domain
% Here the potential is constant outside the
% interior domain
Vpos = 0;
Vneg = 0;

% Range for searching the eigenvalue
specRange = [-10 -1e-2];
% specRange = [-1 -1e-2];
numSpec   = 2048;
specvec   = linspace(specRange(1), specRange(2), numSpec)';

% Number of eigenvalues
numEigs = 4;

% ******************** %
% Mesh interior domain %
% ******************** %
numNodes = 512;
msh = meshObject('uniform', IDb(1), IDb(2), numNodes);
xneg = linspace(ODb(1), IDb(1), ceil(numNodes * (IDb(1) - ODb(1))));
xpos = linspace(IDb(2), ODb(2), ceil(numNodes * (ODb(2) - IDb(2))));

% **************************** %
% Eigenvalue searching process %
% **************************** %
% We introduce the solutions of the half-line problems
% - (u_±)" + (V_± - E) u_± = 0 on Ω_±,
%       - (u_±)' + r_± u_± = 1 at x = a_±
% Since the potential is constant outside the
% interior domain, these solutions have
% an explicit expression
uNeg = @(x, E, RobinNeg) exp( sqrt(Vneg - E) * (x - IDb(1))) / (sqrt(Vneg - E) - RobinNeg); 
uPos = @(x, E, RobinPos) exp(-sqrt(Vpos - E) * (x - IDb(2))) / (sqrt(Vpos - E) - RobinPos); 

% Associated RtR coefficients
RtRneg = @(E, RobinNeg) (RobinNeg + sqrt(Vneg - E)) / (RobinNeg - sqrt(Vneg - E));
RtRpos = @(E, RobinPos) (RobinPos + sqrt(Vpos - E)) / (RobinPos - sqrt(Vpos - E));

% Initialize eigenvalues and eigenvectors
eigvals = zeros(numSpec, numEigs);
eigvecs = zeros(msh.numPoints, numEigs, numSpec);
dEigvecs = zeros(msh.numPoints, numEigs, numSpec);

RobinNegVec = 1i * sqrt(abs(specvec));
RobinPosVec = 1i * sqrt(abs(specvec));

for idI = 1:numSpec

  specVar  = specvec(idI);
  RobinNeg = RobinNegVec(idI);
  RobinPos = RobinPosVec(idI);

  % Compute Robin coefficient involved in interior EVP
  BC.D = diag([-1 +1]);
  BC.A = zeros(2);
  BC.A(1, 1) = RobinNeg * (RtRneg(specVar, RobinNeg) - 1) / (RtRneg(specVar, RobinNeg) + 1);
  BC.A(2, 2) = RobinPos * (RtRpos(specVar, RobinPos) - 1) / (RtRpos(specVar, RobinPos) + 1);

  % Solve eigenvalue problem
  [eigvals(idI, :), eigvecs(:, :, idI), ~, dEigvecs(:, :, idI)] = SchroedingerCellEVP(msh, funP, Vint, BC, numEigs);
  
  % fprintf('%5e\n', abs(BC.D(1, 1)*dEigvecs(1, 1, idI) + BC.A(1, 1)*eigvecs(1, 1, idI)));

  % Plot progress
  fprintf('%3d%%\n', round(100*idI/numSpec));

end

% ********************************************* %
% Plot the difference in absolute value between %
% the eigenvalues and the spectral parameter    %
% ********************************************* %
figure(1);
hold off;
idEig  = 1;
subplot(1, 2, 1);
specIndic = abs(eigvals(:, idEig) - specvec);
plot(specvec, specIndic, 'b');
xlabel('$E$')
title('$E \mapsto |\mu_1 (E) - E|$');
set(gca, 'FontSize', 16);
hold on;

subplot(1, 2, 2);
for idE = 1:numEigs
  plot(specvec, real(eigvals(:, idE)), 'b'); hold on;
end
plot(specvec, specvec, 'r--');
xlabel('$E$');
title('$E \mapsto \Re \mu_1 (E)$');
set(gca, 'FontSize', 16);

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

