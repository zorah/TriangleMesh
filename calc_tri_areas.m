function S_tri = calc_tri_areas(M)
%CALC_TRI_AREAS Calculates the area of each triangle in M.
%   M needs to have fields VERT, TRIV, n and m.
% 
%   copyright (c) 2016 Emanuele Rodolà 

S_tri = zeros(M.n,1);

for k=1:M.m
    e1 = M.VERT(M.TRIV(k,3),:) - M.VERT(M.TRIV(k,1),:);
    e2 = M.VERT(M.TRIV(k,2),:) - M.VERT(M.TRIV(k,1),:);
    S_tri(k) = 0.5*norm(cross(e1,e2));
end

end
