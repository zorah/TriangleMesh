function [ angles ] = triangle_angles( T )
%TRIANGLE_ANGLES Calculates the angles of T in rad. T is 3x3 with each
%column a coordinate. 
% angles - 1x3, ordered corresponding to the vertices
%
%   copyright (c) 2016 Zorah Lähner (laehner@in.tum.de)

edges = triangle_edges(T);

angles = zeros(1,3);

angles(1) = acos( -edges(:,3)'*edges(:,2) / (norm(edges(:,3),2) * norm(edges(:,2),2)) );
angles(2) = acos( -edges(:,1)'*edges(:,3) / (norm(edges(:,1),2) * norm(edges(:,3),2)) );
angles(3) = acos( -edges(:,1)'*edges(:,2) / (norm(edges(:,1),2) * norm(edges(:,2),2)) );

end

