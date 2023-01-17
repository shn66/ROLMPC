function [Xb, Ub] = tightened_constraints_j(UBx, LBx, UBu, LBu, E, K, Xp, Up)
nx=size(UBx,1);
nu=size(UBu,1);
X=Polyhedron('A',[eye(nx);-eye(nx)], 'B', [UBx;-LBx] );
Xh=X-E.outerApprox();
U=Polyhedron('A',[eye(nu);-eye(nu)], 'B', [UBu;-LBu] );
Uh=K*E;
Uh=U-Uh.outerApprox();

a_x=sdpvar(1); v_x=sdpvar(nx,1);
a_u=sdpvar(1); v_u=sdpvar(nu,1);

cost=-a_x-a_u;

r=3;


R=20;
L=6;
dt=0.2;

Xb_vert=Xb.computeVRep.V';
Ub_vert=Ub.computeVRep.V';
nVx=size(Xb_vert,2);
nVu=size(Ub_vert,2);
constr=[a_x>=1; a_u>=1];


x_var=sdpvar(nx,r, nVu);
u_var=sdpvar(nu,r-1,nVx);
for n=1:nVx

    x=[a_x*Xb_vert(:,n)+v_x];
    y=[x(1:nu,1)];
    constr=[constr;
            Xh.A*x(:,1)<=Xh.B];
    for i=1:r-1
        K=2/R*(-0.5+1/(1+exp(5*x(1,i)-100)));
        x=[x x(:,i)+dt*[u_var(1,i,n)*cos(x(3,i))/(1-x(2,i)*K);u_var(1,i,n)*sin(x(3,i));u_var(1,i,n)*(tan(u_var(2,i,n))/L-K*cos(x(3,i))/(1-x(2,i)*K))]];
        y=[y x(1:nu,i+1)];
        constr=[constr;
                Xh.A*x(:,i+1)<=Xh.B;
                Uh.A*u_var(:,i,n)<=Uh.B;
                ];
    end

    constr=[constr; 
            y(2,2)-y(2,1)<=tan(UBx(3))*(1-3/R)*(y(1,2)-y(1,1));
            y(2,2)-y(2,1)<=tan(UBx(3))*(1+1/R)*(y(1,2)-y(1,1));
            y(2,2)-y(2,1)>=tan(LBx(3))*(1-3/R)*(y(1,2)-y(1,1));
            y(2,2)-y(2,1)>=tan(LBx(3))*(1+1/R)*(y(1,2)-y(1,1));
            (1+1/R)^2*(y(1,2)-y(1,1))^2+(y(2,2)-y(2,1))^2<=UBu(1)^2];
end

u_var2=sdpvar(nu,r-1,nVu);
for n=1:nVu
    y=[x_var(1:nu,1,n)];
    constr=[constr;u_var2(:,1,n)==a_u*Ub_vert(:,n)+v_u;
            Xh.A*x_var(:,1,n)<=Xh.B];
    for i=1:r-1
        K=2/R*(-0.5+1/(1+exp(5*x_var(1,i,n)-100)));
        y=[y x_var(1:nu,i+1,n)];
        constr=[constr;
                x_var(:,i+1,n)==x_var(:,i,n)+dt*[u_var2(1,i,n)*cos(x_var(3,i,n))/(1-x_var(2,i,n)*K);u_var2(1,i,n)*sin(x_var(3,i,n));u_var2(1,i,n)*(tan(u_var2(2,i,n))/L-K*cos(x_var(3,i,n))/(1-x_var(2,i,n)*K))]; 
                Xh.A*x_var(:,i+1,n)<=Xh.B;
                Uh.A*u_var2(:,i,n)<=Uh.B;
                ];
    end

    constr=[constr; 
            y(2,2)-y(2,1)<=tan(UBx(3))*(1-3/R)*(y(1,2)-y(1,1));
            y(2,2)-y(2,1)<=tan(UBx(3))*(1+1/R)*(y(1,2)-y(1,1));
            y(2,2)-y(2,1)>=tan(LBx(3))*(1-3/R)*(y(1,2)-y(1,1));
            y(2,2)-y(2,1)>=tan(LBx(3))*(1+1/R)*(y(1,2)-y(1,1));
            (1+1/R)^2*(y(1,2)-y(1,1))^2+(y(2,2)-y(2,1))^2<=UBu(1)^2];
end

ops=sdpsettings('solver','ipopt');
sol=optimize(constr, cost, ops);
    
Xb=value(a_x)*Xh+value(v_x);
Ub=value(a_u)*Uh+value(v_u);

end