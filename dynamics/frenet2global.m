function xg = frenet2global(x, xc, xcg)
N=size(x,2);
xg=[];
for t=1:N
    [ind,dist]=dsearchn(xc(1,:)', x(1,t)');
    xt=xcg(1:2,ind)+x(2,t)*[cos(pi/2+xc(3,ind));sin(pi/2+xc(3,ind))];
    xg= [xg [xt;xc(3,ind)+x(3,t)]];
end
end