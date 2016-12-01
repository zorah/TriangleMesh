function [ M ] = massmatrix( N )
%MASSMATRIX Generates a (sparse) n-by-n mass matrix.
%  Input:
%   - V : n-by-3 matrix contains the 3D coordinates of the n vertices in
%         each of it's rows.
%   - F : m-by-3 matrix contains the indeces of the triangles' vertices.
%
%  Output:
%   - M : n-by-n mass matrix contains the values as described in the
%         lecture.

V = N.VERT;
F = N.TRIV;

if  (size(V,2)~=3)
    V = V';
end
assert(size(V,2)==3);
if  (size(F,2)~=3)
    F = F';
end
assert(size(F,2)==3);
m = size(F,1);
n = size(V,1);

%triangle areas
A = N.tri_areas();

indecesI = [F(:,1);F(:,2);F(:,3);F(:,3);F(:,2);F(:,1)];
indecesJ = [F(:,2);F(:,3);F(:,1);F(:,2);F(:,1);F(:,3)];
values   = [A(:)  ;A(:)  ;A(:)  ;A(:)  ;A(:)  ;A(:)  ]*(1/12);
M = sparse(indecesI, indecesJ, values,n,n);
M = M+sparse(1:n,1:n,sum(M));


end

