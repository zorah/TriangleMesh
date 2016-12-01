function [ grad_f] = trimesh_gradient( M, f )
%TRIMESH_GRADIENT Calculates the gradient of a function f on the triangle mesh M. 
%   f - n x 1
%   M needs to contain the follownig fields:
%   - n - number of vertices
%   - m - number of faces
%   - METRIC_TENSOR - 2 x 2 x m containing the first fundamental form for
%                     each triangle
%   - TRIV - m x 3, indices of the vertices in each triangle
%   - DIFFERENTIAL - 3 x 2 x m containing the differential of each triangle
%   grad_f - m x 3
% 
%   copyright (c) 2016 Zorah Lähner (laehner@in.tum.de)

n = M.n;
m = M.m;

if length(f) ~= n
   error('trimesh_gradient: Number of elements of f is unequal the number of vertices of M.'); 
end

grad_f = zeros(m,3);

% calculate g^-1

g_inv = zeros(size(M.METRIC_TENSOR));

for k = 1:m
   g_inv(:,:,k) = 1/det(squeeze(M.METRIC_TENSOR(:,:,k))) * ...
       [ M.METRIC_TENSOR(2,2,k), -M.METRIC_TENSOR(1,2,k); -M.METRIC_TENSOR(2,1,k), M.METRIC_TENSOR(1,1,k) ];
end

% calculate gradient on the reference triangle (this part has to be
% consistent with the construction of the differential!)

grad_ftilde = zeros(m,2);

for k=1:m
   grad_ftilde(k,1) = f(M.TRIV(k,2)) - f(M.TRIV(k,1));
   grad_ftilde(k,2) = f(M.TRIV(k,3)) - f(M.TRIV(k,1));
end

% gradient = Dx*g^-1*grad_ftilde

for k=1:m
   grad_f(k,:) = squeeze(M.DIFFERENTIAL(:,:,k))*squeeze(g_inv(:,:,k))*grad_ftilde(k,:)';
end

end

