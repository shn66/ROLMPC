function XS = smoother(X)
nx=size(X,1);
XS=[];
for i=1:nx
    XS=[XS;smoothdata(X(i,:),'movmean',2)];
end
end