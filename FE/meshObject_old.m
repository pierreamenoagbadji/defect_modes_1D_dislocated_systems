% meshObject.m
% + ========================================================================= +
% Class for one-dimensional meshes                                            %
%                                                                             %
% A mesh here consists of nodes and elements (segments connecting the nodes)  %
% of a bounded interval. A finer mesh may be generated from a coarser mesh by %
% some refinement procedure.                                                  %
% + ========================================================================= +
classdef meshObject < handle

   properties (SetAccess = protected)

       % Name of the mesh
      name = '';

      % Number of nodes in the mesh
      numNodes = 0;

      % Coordinates of the nodes
      nodes = [];

      % Number of elements in the mesh
      numElements = 0;

      % Matrix which contains the definition of the elements. The i-th row
      % contains the indices of the nodes connected by the i-th element.
      elements = [];

      % 2-sized vector which contains the segment bounds
      bounds = [];

      % 2-sized vector which contains the node indices of the segment bounds
      boundsIds = [];

      % Whenever the mesh has been generated from refining another mesh, this
      % object describes the refinement
      refinementRelation = [];

   end

   methods

      % ======================================================= %
      % Generate a uniform mesh given the bounds of the segment %
      % and the number of nodes                                 %
      % ======================================================= %
      function uniformMesh(mesh, a, b, numNodes)

        % Make sure the mesh contains more than two nodes
        if (numNodes < 2)
          error('Le maillage doit contenir au moins deux elements.');
        end

        % The mesh
        mesh.nodes = linspace(a, b, numNodes).';
        mesh.elements = [(1:numNodes-1).', (2:numNodes).'];
        mesh.numNodes = numNodes;
        mesh.numElements = numNodes - 1;
        mesh.bounds = [a; b];
        mesh.boundsIds = [1; numNodes];
        mesh.refinementRelation = [];

      end

      % ======================================== %
      % Create a mesh from a given set of points %
      % ======================================== %
      function meshFromPoints(mesh, points)

        % The vector containing the nodes has to be
        % a sorted column vector
        points = sort(points(:));


        % Make sure the mesh contains more than two nodes
        if (length(points) < 2)
          error('Le maillage doit contenir au moins deux elements.');
        end

        % The mesh
        mesh.nodes = points;
        mesh.numNodes = length(mesh.nodes);
        mesh.elements = [(1:mesh.numNodes-1).', (2:mesh.numNodes).'];
        mesh.numElements = mesh.numNodes - 1;
        mesh.bounds = [min(mesh.nodes); max(mesh.nodes)];
        mesh.boundsIds = [1; mesh.numNodes];
        mesh.refinementRelation = [];

      end

      % ==================================================== %
      % Construct an interpolated function from nodal values %
      % INPUT: x (Npts x 1) vector containing the points     %
      %        mesh, the mesh object                         %
      %        U (m.numNodes x Nu) the vectors               %
      %                                                      %
      % OUTPUT: val (Npts x Nu)                              %
      % ==================================================== %
      function val = interpFun(x, mesh, U)

        % Make sure U has the adequate number of rows
        if (size(U, 1) ~= mesh.numNodes)

          error(['U a %d lignes, alors qu''il lui en faut %d ',...
                 'pour l''interpolation.'], size(U, 1), mesh.numNodes);

        end

        % Make sure the points belong to the interval
        if (max(max((x - mesh.bounds(1)) .* (mesh.bounds(2) - x) < 0)))

          error('Des elements de x ne sont pas dans le segment de mesh.');

        end

        % x is treated as a column vector
        x0 = x(:);
        Npt = length(x0);
        N = mesh.numNodes;

        % Compute the barycentric coordinate
        % of each point in each segment
        Len = mesh.nodes(2:end) - mesh.nodes(1:end-1); % Length of elements
        tau = (x0 * ones(1, N-1) - ones(Npt, 1) * mesh.nodes(1:end-1).') ./...
                                                  (ones(Npt, 1) * Len.');

        % Find the indices such that tau * (1 - tau) is non-negative
        [~, idElt] = max(tau .* (1 - tau), [], 2);

        % Extract the elements to which each point belongs
        % as well as the corresponding barycentric coordinates
        tau = tau(sub2ind([Npt, N-1], (1:Npt).', idElt)) * ones(1, size(U, 2));

        % Deduce the interpolated value
        val = (1 - tau) .* U(idElt, :) + tau  .* U(idElt + 1, :);

        % Reshape the interpolated value if needed
        if (size(U, 2) == 1)

          val = reshape(val, size(x));

        end

      end

      % ========================================================= %
      % Construct an interpolated function from nodal values      %
      % with a cartesian product of two segments                  %
      % INPUT: x1 (Npts x 1) vector containing the y1 coordinates %
      %        x2 (Npts x 1) vector containing the y2 coordinates %
      %        mesh1, the mesh object in the y1 direction         %
      %        mesh2, the mesh object in the y2 direction         %
      %        U (mesh1.numNodes x mesh2.numNodes) the vectors    %
      %                                                           %
      % OUTPUT: val (Npts x 1)                                    %
      % ========================================================= %
      function val = interpFun2D(x1, x2, mesh1, mesh2, U)

        % Check that U has the adequate number of rows
        if (size(U, 1) ~= mesh1.numNodes)

          error(['U a %d lignes, alors qu''il lui en faut %d ',...
                 'pour l''interpolation.'], size(U, 1), mesh1.numNodes);

        end

        % Check that U has the adequate number of columns
        if (size(U, 2) ~= mesh2.numNodes)

          error(['U a %d colonnes, alors qu''il lui en faut %d ',...
                 'pour l''interpolation.'], size(U, 2), mesh2.numNodes);

        end

        % Make sure the points in x1 belong to the interval meshed by mesh1
        if (max(max((x1 - mesh1.bounds(1)) .* (mesh1.bounds(2) - x1) < -eps)))

          error('Des elements de x1 ne sont pas dans le segment de mesh1.');

        end

        % Make sure the points in x2 belong to the interval meshed by mesh2
        if (max(max((x2 - mesh2.bounds(1)) .* (mesh2.bounds(2) - x2) < -eps)))

          error('Des elements de x2 ne sont pas dans le segment de mesh2.');

        end

        % x1 and x2 are treated as column vectors
        x1 = x1(:);
        x2 = x2(:);

        Npt1 = length(x1); N1 = mesh1.numNodes;
        Npt2 = length(x2); N2 = mesh2.numNodes;

        % Compute the barycentric coordinate
        % of each point in each segment
        Len1 = mesh1.nodes(2:end) - mesh1.nodes(1:end-1);
        Len2 = mesh2.nodes(2:end) - mesh2.nodes(1:end-1);

        tau1 = (x1*ones(1, N1-1) - ones(Npt1, 1)*mesh1.nodes(1:end-1).') ./...
                                                (ones(Npt1, 1) * Len1.');
        tau2 = (x2*ones(1, N2-1) - ones(Npt2, 1)*mesh2.nodes(1:end-1).') ./...
                                                (ones(Npt2, 1) * Len2.');

        % Find the indices such that tau * (1 - tau) is non-negative
        [~, idElt1] = max(tau1 .* (1 - tau1), [], 2);
        [~, idElt2] = max(tau2 .* (1 - tau2), [], 2);

        % Extract the elements to which each point belongs
        % as well as the corresponding barycentric coordinates
        tau1 = tau1(sub2ind([Npt1, N1-1], (1:Npt1).', idElt1));
        tau2 = tau2(sub2ind([Npt2, N2-1], (1:Npt2).', idElt2));

        % Deduce the interpolated value
        val = (1 - tau1) .* (1 - tau2) .* U(idElt1 + (idElt2-1) * N1)   + ...
              (1 - tau1) .* tau2 .* U(idElt1 + idElt2 * N1) + ...
                   tau1 .* (1 - tau2) .* U(idElt1 + 1 + (idElt2-1) * N1) + ...
                   tau1 .* tau2  .* U(idElt1 + 1 + idElt2 * N1);

      end

   end


end
