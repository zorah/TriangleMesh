function [ edges ] = triangle_edges( T )
%TRIANGLE_EDGES Calculates the vectores a, b, c describing the edges of T.
%edges = [a, b, c] 3x3 where a is the vector BC opposing A and so on.
% T is a 3x3 matrix of form [A, B, C].
%
% copyright (c) 2016 Zorah Lähner (laehner@in.tum.de)

edges = [T(:,3) - T(:,2), T(:,1) - T(:,3), T(:,2) - T(:,1)];

end

