function [ Dx ] = trimesh_differential( M )
%TRIMESH_DIFFERENTIAL Calculates the differential of each triangle
% on a triangle mesh consisting of M.VERT for the vertices and
% M.TRIV for the triangles.
% The output is a 3x2xm matrix.
%
%   copyright (c) 2016 Zorah Lähner (laehner@in.tum.de)

Dx = zeros(3,2,size(M.TRIV,1));

for k=1:size(M.TRIV,1)
   Dx_1 = M.VERT(M.TRIV(k,2),:)' - M.VERT(M.TRIV(k,1),:)';
   Dx_2 = M.VERT(M.TRIV(k,3),:)' - M.VERT(M.TRIV(k,1),:)';
   Dx(:,:,k) = [Dx_1 Dx_2];
end

end

