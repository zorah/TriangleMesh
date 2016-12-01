function [ distances, sources ] = exact_geodesics( M, sources, targets )
%EXACT_GEODESICS Returns the distances of all vertices in targets to the
%closest point in sources. 
%
% M - triangle mesh, needs to contain M.VERT(nx3) and M.TRIV(mx3) (1-based)
% sources - 1-based ids of the source points
% targets - 1-based ids of the target points
%
% distances[i] - distance of targets[i] to the closest point in sources
% sources[i] - 1-based id of the closest point in sources to targets[i]

[distances, sources] = geodesic_exact(M.VERT', (M.TRIV-1)', sources-1, targets-1);
sources = sources + 1;

end

