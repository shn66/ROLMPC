function [Xb, Ub] = tightened_constraints_0(X,U, UBx, LBx, UBu, LBu)
nx=3;
nu=2;
Xh=X;
Uh=U;

r=3;

R=20;
L=2;
dt=0.2;

Xb_vert=Xh.computeVRep.V';
Ub_vert=Uh.computeVRep.V';
nVx=size(Xb_vert,2);
nVu=size(Ub_vert,2);

cost = 0;
constr = [];

a_x=sdpvar(1); v_x=sdpvar(nx,1);

constr=[a_x>=0.9];

cost = cost -exp(a_x) + v_x'*v_x;
u_var=sdpvar(nu,r-1,nVx);
for n=1:nVx
    x=[a_x*Xb_vert(:,n)+v_x];
    y=[x(1:nu,1)];
    constr=[constr;
            Xh.A*x(:,1)<=Xh.b];
    for i=1:r-1
        K=2/R/pi*atan(100-.5*x(1,i)^2);
        x=[x x(:,i)+dt*[u_var(1,i,n)*cos(x(3,i))/(1-x(2,i)*K);u_var(1,i,n)*sin(x(3,i));u_var(1,i,n)*(tan(u_var(2,i,n))/L-K*cos(x(3,i))/(1-x(2,i)*K))]];
        y=[y x(1:nu,i+1)];
%         constr=[constr;
%                 Xh.A*x(:,i+1)<=Xh.b;
%                 Uh.A*u_var(:,i,n)<=Uh.b;
%                 ];
    end

    constr=[constr; 
            
            y(2,2)-y(2,1)<=tan(UBx(3))*(1-1/R)*(y(1,2)-y(1,1));
            y(2,2)-y(2,1)>=tan(LBx(3))*(1+1/R)*(y(1,2)-y(1,1));
%             (1+1/R)^2*(y(1,2)-y(1,1))^2+(y(2,2)-y(2,1))^2<=UBu(1)^2
            ];
end


a_u=sdpvar(nu,1); v_u=sdpvar(nu,1);
x_var=sdpvar(nx,r, nVu);
u_var2=sdpvar(nu,r-1,nVu);
constr = [constr; a_u>=0.9];

cost = cost - 1*a_u(1,1)-1*a_u(2,1) + 10*v_u'*v_u;
for n=1:nVu
    y=x_var(1:nu,1,n);
    constr=[constr;
        u_var2(:,1,n)==diag(a_u)*Ub_vert(:,n)+v_u]; 
%         Xh.A*x_var(:,1,n)<=Xh.b];
%         Uh.A*u_var2(:,1,n)<=Uh.b;

    for i=1:r-1
        
        K=2/R/pi*atan(100-.5*x_var(1,i,n)^2);
        y=[y x_var(1:nu,i+1,n)];
        x_tp1 = x_var(:,i,n)+dt*[u_var2(1,i,n)*cos(x_var(3,i,n))/(1-x_var(2,i,n)*K);u_var2(1,i,n)*sin(x_var(3,i,n));u_var2(1,i,n)*(tan(u_var2(2,i,n))/L-K*cos(x_var(3,i,n))/(1-x_var(2,i,n)*K))];
        constr=[constr; x_var(:,i+1,n)==x_tp1];
%                 Xh.A*x_var(:,i+1,n)<=Xh.b];
    end

    constr=[constr; 
            sqrt((1+1/R)^2*(y(1,2)-y(1,1))^2+(y(2,2)-y(2,1))^2)<=UBu(1)];
end

ops=sdpsettings('solver','ipopt');
sol=optimize(constr, cost, ops);
  
Xb=value(a_x)*Xh+value(v_x); Xb = Xb.intersect(X);
Ub=value(diag(a_u))*Uh+value(v_u); 
Ub = Ub.intersect(U);

end