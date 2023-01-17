function [xn, un] = make_nominal(X,U, E, x, x_guess, u_guess)
nx=3;
nu=2;
xn=sdpvar(nx, size(x,2));
un=sdpvar(nu, size(x,2)-1);
tn=sdpvar(1, size(x,2));
assign(xn, x_guess);
assign(un,u_guess);
cost=0;

R=20;
L=2;
dt=0.2;
s_f=40;
constr=[];
for t=1:size(x,2)-1
    cost=cost+10*tn(t);
    constr=[constr;s_f-xn(1,t+1)<=tn(t);tn(t)>=0];
    constr=[constr;
            E.A*(x(:,t)-xn(:,t))<=E.b;
            X.A*[2;xn(2:end,t+1)]<=X.b;
            U.A*un(:,t)<=U.b];
    K=2/R/pi*atan(100-.5*xn(1,t)^2);
    constr=[constr;
                xn(:,t+1)==xn(:,t)+dt*[un(1,t)*cos(xn(3,t))/(1-xn(2,t)*K);un(1,t)*sin(xn(3,t));un(1,t)*(sin(un(2,t))/L-K*cos(xn(3,t))/(1-xn(2,t)*K))]; 
                ];

end

    
ops=sdpsettings('solver','ipopt', 'verbose',1);
ops.ipopt.max_iter=3000;
options.ipopt.tol=1e-2;
sol=optimize(constr, cost , ops);
    
xn=double(xn);
un=double(un);

end