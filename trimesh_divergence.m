function [ res ] = trimesh_divergence( M, V )
%TRIMESH_DIVERGENCE Calculates the divergence of a vector field V defined
%on M. 
%   V - 3 x m
%   M needs to contain the following fields:
%   - n - number of vertices
%   - TRIV - indices of vertices of each triangle
%   res - n x 1
%
%   copyright (c) 2016 Zorah Lähner (laehner@in.tum.de) 

res = zeros(M.n, 1);

% formula for divergence as found in geodesics in heat, crane et al, 2014
for i=1:M.n
    triangles = find(M.TRIV(:,1) == i | M.TRIV(:,2) == i | M.TRIV(:,3) == i);
    sum = 0;
    for j=1:length(triangles)
       point_ids = M.TRIV(triangles(j), M.TRIV(triangles(j),:) ~= i); % points in triangle not i
       e1 = M.VERT(point_ids(1),:) - M.VERT(i, :);
       e2 = M.VERT(point_ids(2),:) - M.VERT(i, :);
       e3 = M.VERT(point_ids(1),:) - M.VERT(point_ids(2),:);
       theta2 = cotd(acosd((e1*e3')./(norm(e1)*norm(e3))));
       theta1 = cotd(acosd((-e2*e3')./(norm(e2)*norm(e3))));
       sum = sum + theta1 * (e1*(V(triangles(j),:)')) + theta2 * (e2 * (V(triangles(j),:)'));
    end
    res(i) = 0.5*sum;
end

end

