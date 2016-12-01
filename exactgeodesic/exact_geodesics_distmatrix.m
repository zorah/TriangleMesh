function [ dist ] = exact_geodesics_distmatrix( M )
%EXACT_GEODESICS_DISTMATRIX Returns all pointwise distances between
%vertices on M.
% 
% M - triangle mesh, needs to contain M.VERT (nx3) and M.TRIV (mx3)
% source - 1-based id of the source point
%
% dist - n x n, distances of each point to each other pount

dist = geodesic_exact_distmatrix(M.VERT', (M.TRIV-1)');

end

