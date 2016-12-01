function [ centers ] = triangle_edgecenters( T )
%TRIANGLE_EDGECENTERS Returns the edge centers of triangle T (3x3, [A, B,
%C]) is form [a, b, c] such that a is the center opposing A.
%
% copyright (c) 2016 Zorah Lähner (laehner@in.tum.de)

centers = [ (T(:,2) + T(:,3))/2, (T(:,1) + T(:,3))/2, (T(:,1) + T(:,2))/2 ];

end

