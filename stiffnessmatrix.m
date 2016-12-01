function [ S ] = stiffnessmatrix( M )
%STIFFNESSMATRIX Returns the (sparse) stiffness matrix created by using the
%cotan-scheme. 
%   M needs to contain the following fields: 
%   - VERT - n x 3, 3D coordinates of vertices
%   - TRIV - m x 3, indices of vertices in each triangle
%   S - n x n, sparse,
% 
%   copyright (c) 2016 Matthias Vestner

V = M.VERT;
F = M.TRIV;

if  (size(V,2)~=3)
    V = V';
end
assert(size(V,2)==3);
n = length(V);
if  (size(F,2)~=3)
    F = F';
end
assert(size(F,2)==3);

% Compute matrix containing the cotangens-values
vm = V(F(:,1),:);
vn = V(F(:,2),:);
vk = V(F(:,3),:);

size(vm);
% calculate edge lengths
vmn = vm-vn;
vnk = vn-vk;
vkm = vk-vm;

% inner products
ipk = sum(-vkm.*vnk,2);
ipn = sum(vmn.*(-vnk),2);
ipm = sum(-vmn.*vkm,2);

cotam = ipm ./ sqrt(sum(vmn.^2,2) .* sum(vkm.^2,2) - ipm.^2);
cotan = ipn ./ sqrt(sum(vmn.^2,2) .* sum(vnk.^2,2) - ipn.^2);
cotak = ipk ./ sqrt(sum(vkm.^2,2) .* sum(vnk.^2,2) - ipk.^2);

A = 0.5*[cotam cotan cotak];

i = [    F(:,1);    F(:,1);  F(:,1);  F(:,2);      F(:,2);    F(:,2);  F(:,3);  F(:,3);      F(:,3)];
j = [    F(:,1);    F(:,2);  F(:,3);  F(:,1);      F(:,2);    F(:,3);  F(:,1);  F(:,2);      F(:,3)];
s = [A(:,3)+A(:,2); -A(:,3); -A(:,2); -A(:,3); A(:,3)+A(:,1); -A(:,1); -A(:,2); -A(:,1); A(:,1)+A(:,2)];

S = -sparse(i,j,s,n,n);

end

