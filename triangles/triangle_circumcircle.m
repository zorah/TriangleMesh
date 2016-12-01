function [ center, radius ] = triangle_circumcircle( T )
%TRIANGLE_CIRCUMCIRCLE Calculates the center and radius of the circumcircle
%of T (3x3, [A, B, C]).
%
% copyright (c) 2016 Zorah Lähner (laehner@in.tum.de)

[T, tr, rot] = triangle_toplane(T);
lengths = sum(T(1:2, :).^2, 1);

Sx = 0.5 * det([ lengths', T(2,:)', ones(3,1) ]);
Sy = 0.5 * det([ T(1,:)', lengths', ones(3,1) ]);

a = det([ T(1:2, :)', ones(3,1) ]);
b = det([ T(1:2, :)', lengths' ]); 

c = [Sx; Sy; 0] ./ a;
radius = ( b/a + (Sx^2 + Sy^2)/(a^2));
center = (rot' * c) + tr;

end

