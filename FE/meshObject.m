%> @file meshObject.m
%> @brief Contains the meshObject class.
% =========================================================================== %
%> @brief Class for one-dimensional meshes
%>
%> A mesh here consists of nodes and elements (segments connecting the nodes)
%> of a bounded interval.
% =========================================================================== %
classdef meshObject < handle
  % meshObject < handle

  properties (SetAccess = protected)

    %> @brief Dimension of the geometry
    dimension = 1;

    %> @brief Total number of vertices in the mesh
    numPoints = 0;

    %> @brief Matrix that contains the coordinates of the points. The k-th
    %> point is represented by the k-th line
    points = [];

    %> @brief Number of segments
    %> NOTE: a segment might be a volumic or surfacic entity depending on
    %> the dimension
    numSegments = 0;

    %> @brief Matrix that contains the definition of the segments. The k-th
    %> segment is represented by the k-th line. Each line contains the indices
    %> of the nodes that form the segment
    segments = [];

    %> @brief 2-sized vector which contains the segment bounds
    bounds = [];

    %> @brief 2-sized vector which contains the node indices of the segment bounds
    boundsIds = [];

  end

  methods

    % ============= %
    % Create a mesh %
    % ============= %
    function mesh = meshObject(varargin)
      % meshObject constructor for segment mesh
      % The first argument (string) indicates the type of mesh:
      %     - for a uniform mesh, the syntax is
      %         mesh = meshes.meshObject('uniform', a, b, N)
      %       a, b are the bounds of the segment and N the number of nodes
      %
      %     - for a mesh generated from points, the syntax is
      %         mesh = meshes.meshObject('custom', points)
      %       points is a list of the nodes
      if (nargin < 1)

        % Default argmuments
        numPoints = 8;
        points = linspace(0, 1, numPoints).';

      elseif strcmpi(varargin{1}, 'uniform')

        % Generate a mesh of equispaced points
        xmin = varargin{2};
        xmax = varargin{3};
        numPoints = varargin{4};
        points = linspace(xmin, xmax, numPoints).';

      elseif strcmpi(varargin{1}, 'custom')

        % Generate a mesh from a list of vertices
        points = unique(varargin{2});
        numPoints = length(points);
        xmin = points(1);
        xmax = points(end);

      else

        error(['Variable ''', varargin{1}, ''' unknown. ',...
               'Options are ''uniform'' and ''custom''.']);

      end

      % Properties of the mesh
      mesh.dimension = 1;
      mesh.numPoints = numPoints;
      mesh.points = points(:);
      mesh.numSegments = numPoints - 1;
      mesh.segments = [(1:numPoints-1).', (2:numPoints).'];
      mesh.bounds = [xmin; xmax];
      mesh.boundsIds = [1; numPoints];

    end

  end

end
