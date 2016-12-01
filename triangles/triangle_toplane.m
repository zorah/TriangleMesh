function [ S, tr, rot ] = triangle_toplane( T )
%TRIANGLE_TOPLANE Rigidly moves the triangle T of form [A, B, C] 3x3 onto
%the x-y plane such that A is in the origin and the last coordinate of B 
%and C is zero.
% S = rot * (T - tr);
%
% copyright (c) 2016 Zorah Lähner (laehner@in.tum.de)

% move A to origin
tr = T(:,1);
T = T - repmat(tr, 1, 3);

% rotate on x-y plane
u = T(:,2) ./ norm(T(:,2), 2);
w = cross(u, T(:,3));
w = w ./ norm(w, 2);
v = cross(u, w);

rot = [u, v, w]';

S = rot*T;

end

