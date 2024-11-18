function P = polytope_approx(x,y)
centroid_x = mean(x);
centroid_y = mean(y);
[theta,r] = cart2pol(x-centroid_x,y-centroid_y);
% sort theta in ascending order
[theta_sorted,ind] = sort(theta);
r_sorted = r(ind);
% remove duplicates before interpolation
[theta_unique,IA,IC] = unique(theta_sorted);
r_unique = r_sorted(IA);


% load('bound.mat')
% 
% the=linspace(min(theta_unique), max(theta_unique),100);
% r_=bound(the)'+0.01;

[u,v] = pol2cart(theta_unique,r_unique);
u = u + centroid_x;
v = v + centroid_y;

Verts=[u;v]';

P=Polyhedron('V',Verts);

end

