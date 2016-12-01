function [ ] = isolines_trianglemesh( M, f, n )
%ISOLINES_TRIANGLEMESH Plots n equidistant isolines of the functions f on M. 
%   M needs to contain the fields VERT and TRIV.
%   f - function defined on vertices of M
%   n - number of isolines
%
%   code (slightly modificated) by Alec Jacobson (alecjacobson.com)

F = M.TRIV;
V = M.VERT;
S = f;

colormap(jet(n));
trisurf(F,V(:,1),V(:,2),V(:,3),'CData',S,'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
axis equal;
[LS,LD,I] = isolines(V,F,S,linspace(min(S),max(S),n+1));
hold on;
plot3([LS(:,1) LD(:,1)]',[LS(:,2) LD(:,2)]',[LS(:,3) LD(:,3)]','k','LineWidth',3);
hold off
set(gcf, 'Color', [1,1,1]);
set(gca, 'visible', 'off');


end

