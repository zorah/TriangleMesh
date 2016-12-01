function [ center, areas ] = triangle_hybridvoronoi( T )
%TRIANGLE_HYBRIDVORONOI Calculates the center and areas of the hybrid
%voronoi cells in the triangle T (3x3, coordiantes in columns). Hybrid
%means the center is set to the middle of the cut edge if the center would
%lie outside of the triangle. 
% center - 3x1, inside or on the edge of T
% areas - 1x3, area of the voronoi cell belonging to each vertex
%
% edges can be easily computed by connecting the center to the midpoint of
% each edge of T
%
% copyright (c) 2016 Zorah Lähner (laehner@in.tum.de)
% using instructions and visualizations from alecjacobson.com

% calculate center of circumcircle
[center, radius] = triangle_circumcircle( T );

edge_centers = triangle_edgecenters( T );
areas = zeros(1,3);

% check if center is inside of T 
[maxi, obtuse_id] = max(triangle_edgelengths(T));
if maxi/2 <= radius
    % hybrid
    center = sum(T(:, 1:3 ~= obtuse_id), 2)/2;
    for i=1:3
       opp = find(1:3 ~= i);
       if i ~= obtuse_id
           areas(i) = triangle_area( [T(:,i), edge_centers(:,opp(2)), edge_centers(:,opp(1)) ] );
       else
          areas(i) = triangle_area( [T(:,i), center, edge_centers(:,opp(1))] ) ...
              + triangle_area( [T(:,i), center, edge_centers(:,opp(2))] );
       end
    end
else
    % normal voronoi
    for i=1:3
        opp = find(1:3 ~= i);
        areas(i) = triangle_area( [T(:,i), center, edge_centers(:,opp(1))] ) ...
              + triangle_area( [T(:,i), center, edge_centers(:,opp(2))] );
    end
end

end

