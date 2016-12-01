function [ dist ] = exact_geodesics_distmap( M, source )
%EXACT_GEODESICS_DISTMAP Returns the distance of each point on M to the source
%point source. 
% 
% M - triangle mesh, needs to contain M.VERT (nx3) and M.TRIV (mx3)
% source - 1-based id of the source point
%
% dist - distances of each point in M to source

dist = geodesic_exact_distmap(M.VERT', (M.TRIV-1)', source-1);

end

