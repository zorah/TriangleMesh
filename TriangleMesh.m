classdef TriangleMesh < handle
    %TRIANGLEMESH Class for Triangle Meshes. 
    %   Loads triangle mesh from a given file, can perform basic
    %   operations and automatically stores computed data for later 
    %   usage (if wanted). 
    % 
    %   The fields VERT, TRIV, n and m are always filled, other fields
    %   are empty by default. If fields are accessed via the corresponding
    %   function, it will check if the data already exists and read it
    %   from there or compute it if necessary. 
    %   Precomputed data can be safely removed by setting the field to [].
    % 
    %   copyright (c) 2016 Zorah Lähner (laehner@in.tum.de)
    
    properties
        VERT; % 3D vertex coordinates, n x 3
        TRIV; % vertex indices of triangles, m x 3
        n; % #vertices
        m; % #faces
        EVALS; % ascending eigenvalues of the LBO
        EVECS; % eigenvectors of the LBO, sorted like EVALS
        X; % first column of VERT, empty by default
        Y; % second column of VERT, empty by default
        Z; % third column of VERT, empty by default
        MASS; % mass matrix, n x n, sparse
        MASS_MODE = 'barycentric'; % type of mass matrix, 'barycentric', 'voroni' or 'full'
        STIFFNESS; % stiffness matrix
        BOUNDARY; % indices of boundary vertices
        TRI_AREAS; % areas of triangles
        TRI_CENTER; % 3D coordinates of the barycenter of each triangle
        DIFFERENTIAL; % differential of each triangle, 3 x 2 x m
        METRIC_TENSOR; % metric tensor of each triangle, 2 x 2 x m
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CONSTRUCTOR 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Reads a triangle mesh out of a file. Possible types are
        % .off, .ply and tosca .mat format (reading into a surface
        % struct with TRIV, X, Y and Z). Additionally, the first input
        % can be a struct with fields VERT and TRIV that will be converted
        % into a TriangleMesh type.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = TriangleMesh(filename)
            
            if isstruct(filename) % construct TriangleMesh from struct
                obj.VERT = filename.VERT;
                obj.TRIV = filename.TRIV;
            else % read mesh from file
                % load file
                [~, ~, ext] = fileparts(filename);
                if strcmp('.off', ext)
                    [vert, tri] = read_off(filename);
                    obj.VERT = vert';
                    obj.TRIV = tri';
                elseif strcmp('.ply', ext)
                    [vert, tri] = read_ply(filename);
                    obj.VERT = vert;
                    obj.TRIV = tri;
                elseif strcmp('.mat', ext)
                    load(filename);
                    obj.TRIV = surface.TRIV;
                    obj.VERT = [surface.X, surface.Y, surface.Z];
                else
                    error('TriangleMesh: The file extension %s is not known.', ext);
                end
                
            end
            
            obj.n = size(obj.VERT, 1);
            obj.m = size(obj.TRIV, 1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % READ VERTICES, X, Y, Z
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns 1st, 2nd or 3rd columns of obj.VERT, respectively.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = x(obj)
           res = obj.VERT(:,1);
        end
        
        function res = y(obj)
           res = obj.VERT(:,2);
        end
        
        function res = z(obj)
           res = obj.VERT(:,3);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % POPULATE_XYZ
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fills the properties X, Y, Z with the columns of 
        % VERT. The entries will not automatically change
        % if VERT is changed
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = populate_xyz(obj)
            obj.X = obj.VERT(:,1);
            obj.Y = obj.VERT(:,2);
            obj.Z = obj.VERT(:,3);
        end
     
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MASS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns the mass matrix.
        % If the input 'mode' is not set, the property 
        % MASS_MODE is used. Possible are 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = mass(obj, mode)
            if ~exist('mode', 'var')
                mode = obj.MASS_MODE;
            end
            
            if (isempty(obj.MASS) || ~strcmp(mode, obj.MASS_MODE))
                obj.MASS_MODE = mode;
                
                if (strcmp(mode, 'barycentric'))
                    obj.MASS= sparse(1:obj.n, 1:obj.n, vertex_voronoi(obj), obj.n, obj.n);
                elseif strcmp(mode, 'voronoi')
                    obj.MASS = sparse(1:obj.n, 1:obj.n, zeros(1,obj.n), obj.n, obj.n);
                    for i=1:obj.m
                        ids = obj.TRIV(i,:);
                        [~, areas] = triangle_hybridvoronoi( obj.VERT(ids, :)' );
                        for j=1:3
                            obj.MASS(ids(j), ids(j)) = obj.MASS(ids(j), ids(j)) + areas(j);
                        end
                    end
                elseif strcmp(mode, 'full')
                    obj.MASS = massmatrix(obj);
                else
                    error('The mass_mode is set to incompatible type %s. Possible are barycentric, voronoi and full. \n', obj.mass_mode);
                end
            end
            res = obj.MASS;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STIFFNESS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns the stiffness/cotan matrix.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = stiffness(obj)
            if isempty(obj.STIFFNESS)
                obj.STIFFNESS = stiffnessmatrix(obj);
            end
            
            res = obj.STIFFNESS;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EIGEN
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns n eigenvectors and corresponding eigen-
        % values of the laplacian of this mesh. 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [evecs, evals] = eigen(obj, n)
           
           if (n > length(obj.EVALS))
               options.disp = 0;
               M = obj.mass();
               S = obj.stiffness();
               [evecs, obj.EVALS] = eigs(S, M, n, 1e-10, options);
               obj.EVALS = -diag(obj.EVALS);
               [obj.EVALS, sortIDs] = sort(obj.EVALS);
               obj.EVECS = evecs(:, sortIDs);
           end
           
           evecs = obj.EVECS(:,1:n);
           evals = obj.EVALS(1:n);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EVALS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns n eigenvalues of the laplacian of this mesh. 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = evals(obj, n)
           if n > length(obj.EVALS)
              eigen(obj, n); 
           end
           
           res = obj.EVALS(1:n);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EVECS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns n eigenvectors of the laplacian of this mesh. 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = evecs(obj, n)
           if n > length(obj.EVALS)
               eigen(obj, n);
           end
           
           res = obj.EVECS(:,1:n);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LAPLACIAN
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns the laplacian. The matrix will not be stored
        % for later use due to memory reasons.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = laplacian(obj)
            res = inv(obj.mass()) * obj.stiffness();
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GRADIENT
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns the gradient of a function f defined on the vertices.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = gradient(obj, f)
            if isempty(obj.DIFFERENTIAL)
                obj.DIFFERENTIAL = trimesh_differential(obj);
            end
            if isempty(obj.METRIC_TENSOR)
                obj.METRIC_TENSOR = trimesh_fff(obj.DIFFERENTIAL);
            end
            res = trimesh_gradient(obj, f);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DIVERGENCE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns the divergence of a vector field V defined on the faces.
        % V must be m x 3.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = divergence(obj, V)
            res = trimesh_divergence(obj, V);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FFF
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns the metric tensor / first fundamental form on each
        % triangle. Output is 2 x 2 x m.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function g = fff(obj)
           if isempty(obj.DIFFERENTIAL)
                obj.DIFFERENTIAL = trimesh_differential(obj);
           end
           if isempty(obj.METRIC_TENSOR)
                obj.METRIC_TENSOR = trimesh_fff(obj.DIFFERENTIAL);
           end
           g = obj.METRIC_TENSOR;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DIFFERENTIAL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns the differential on each
        % triangle. Output is 3 x 2 x m.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Dx = differential(obj)
           if isempty(obj.DIFFERENTIAL)
                obj.DIFFERENTIAL = trimesh_differential(obj);
           end
           Dx = obj.DIFFERENTIAL;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TRI_AREAS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns the areas of all triangles.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = tri_areas(obj)
            if isempty(obj.TRI_AREAS)
                obj.TRI_AREAS = calc_tri_areas(obj);
            end
            
            res = obj.TRI_AREAS;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TRI_CENTER
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns the barycenters of all triangles.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = tri_center(obj)
           if isempty(obj.TRI_CENTER)
              obj.TRI_CENTER = 1/3 * (obj.VERT(obj.TRIV(:,1),:) + obj.VERT(obj.TRIV(:,2),:) + obj.VERT(obj.TRIV(:,3),:));
           end
           
           res = obj.TRI_CENTER;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BOUNDARY_VERTICES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns a logical array marking the boundary vertices.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = boundary_vertices(obj)
            if isempty(obj.BOUNDARY)
                [~, ~, ~, ind] = get_boundary(obj.TRIV); 
                obj.BOUNDARY = false(obj.n, 1);
                obj.BOUNDARY(ind) = true;
            end
           res = obj.BOUNDARY;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NONBOUNDARY_VERTICES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns a logical array marking the non-boundary vertices.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = nonboundary_vertices(obj)
           res = ~obj.boundary_vertices();
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PLOT_SCALAR_MAP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots a scalar function on the mesh. Default colormap is jet.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = plot_scalar_function(obj, f )
            trisurf(obj.TRIV, obj.x(), obj.y(), obj.z(), f, 'EdgeColor', 'none'),
            axis equal, caxis auto,
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PLOT_VECTOR_FIELD
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots a vector field on the mesh. V must be m x 3. f will be
        % plotted on the underlying mesh, if f does not exist the coloring
        % will be constant.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = plot_vector_field(obj, V, f)
            if ~exist('f', 'var')
                f = ones(obj.n, 1);
            end
            trisurf(obj.TRIV, obj.x(), obj.y(), obj.z(), f, 'EdgeColor', 'none'),
            axis equal, caxis auto, hold on,
            C = obj.tri_center();
            quiver3(C(:,1), C(:,2), C(:,3), V(:,1), V(:,2), V(:,3)),
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % VISUALIZE_POINTS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots a vector field on the mesh. V must be m x 3.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = visualize_points(obj, ind)
            trisurf(obj.TRIV, obj.x(), obj.y(), obj.z(), ones(obj.n, 1), 'EdgeColor', 'none'),
            axis equal, caxis auto, hold on,
            scatter3(obj.VERT(ind,1), obj.VERT(ind,2), obj.VERT(ind,3)),
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % VISUALIZE_ISOLINES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots isolines of the function f on the mesh. n equidistant
        % isolines will be plotted.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = visualize_isolines(obj, f, n)
            isolines_trianglemesh(obj, f, n);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PICK_POINTS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Opens a figure where the vertex ids will be shown when selected
        % with the cursor.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = pick_points(obj)
            function [txt] = myupdatefcn(~,event_obj)
                pos = get(event_obj,'Position');
                [~,idx] = ismember(pos,obj.VERT,'rows');
                txt = {['idx: ',num2str(idx)]};
            end
            
            figure, plot_mesh(obj)
            
            h = datacursormode;
            set(h,'SnapToDataVertex','on');
            set(h,'UpdateFcn',@myupdatefcn);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DIJKSTRA
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Computes the distance between the vertices with index a and b
        % using the Dijkstra algorithm. If inputs are
        % - a and b - distance between them
        % - only - distance map of a
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dist = dijkstra(obj, a, b)
            
            if exist('b', 'var')
                dist = dijkstra_geodesic(obj, a, b);
            else
                dist = dijkstra_geodesic(obj, a);
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EXACT_GEODESIC
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Computes the distance between vertices using the exact polyhedral 
        % method. If inputs are
        % - a and b - distance between them, can be arrays
        % - only a - distance map of a
        % - neither - complete distance matrix (use with care)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dist = exact_geodesic(obj, a, b)
            
            if exist('b', 'var')
                dist = exact_geodesics(obj, a, b);
            elseif exist('a', 'var')
                dist = exact_geodesics_distmap(obj, a);
            else
                dist = exact_geodesics_distmatrix(obj);
            end
            
        end
        
    end
    
end

