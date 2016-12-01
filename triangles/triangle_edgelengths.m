function [ length ] = triangle_edgelengths( T )
%TRIANGLE_EDGELENGTHS Calculates the edge lengths of T. T is 3x3 of form
%[A, B, C] and length 1x3 of form [length(a), length(b), length(c)] such
%that a is opposing the vertex A and so on.
%
% copyright (c) 2016 Zorah Lähner (laehner@in.tum.de)

length = sqrt(sum(triangle_edges(T).^2, 1));

end

