function [ area ] = triangle_area( T )
%TRIANGLE_AREA Calculates the area of T given by a 3x3 with a 3D coordinate
%in each column.
%
% copyright (c) 2016 Zorah Lähner (laehner@in.tum.de)

e1 = T(:,3) - T(:,1);
e2 = T(:,2) - T(:,1);
area = 0.5*norm(cross(e1,e2));

end

