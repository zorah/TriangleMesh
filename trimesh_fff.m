function [ fff ] = trimesh_fff( Dx )
%TRIMESH_FFF Calculates the metric tensor out of the differential of 
% a triangle mesh. The differential must be given triangle-wise
% by a 3x2xm matrix. 
% The output is a 2x2xm matrix containing one 2x2 metric tensor for
% each triangle.
% 
%   copyright (c) 2016 Zorah Lähner (laehner@in.tum.de)

fff = zeros(2,2,size(Dx,3));

for k=1:size(Dx,3)
   fff_k = squeeze(Dx(:,:,k))' * squeeze(Dx(:,:,k));
   fff(:,:,k) = fff_k;
end

end

