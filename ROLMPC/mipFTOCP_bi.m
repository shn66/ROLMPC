 function [ x, u, x_guess, u_guess, lmbd_guess, t_guess, sol_time ] = mipFTOCP_bi( x_t , xn_tm1, time, N, Ny,C, P, Qfun, SS, s_f, X, U, E, x_guess, u_guess, lmbd_guess, t_guess, solver)
% Define Yalmip Variables

x=sdpvar(3,N+1);
t=sdpvar(1,N);
u=sdpvar(2,N+1);
ny=size(C,1);

eps=0.0001;
lambda = binvar(length(Qfun), 1); % Number of multipliers used for the convex hull
R=20;
L=2;
dt=0.2;

% Select state and input constraints matrices from polyhedrons
Hx  = X.A;  bx  = X.b;
Hu  = U.A;  bu  = U.b;
He  = E.A;  be  = E.b;
assign(x, x_guess);
assign(u, u_guess);
assign(lambda, lmbd_guess);
assign(t, t_guess);
% ======= Constraints Definition ======

% Initial Condition
if time==1
    Constraints = [He*(x_t-x(:,1))<=be; u(:,1)==[0;0]; t>=0];
else
    K=2/R/pi*atan(100-.5*xn_tm1(1,1)^2);
    Constraints = [He*(x_t-x(:,1))<=be;
                   Hu*u(:,1)<=bu; 
                   x(:,1) == xn_tm1+dt*[u(1,1)*cos(xn_tm1(3,1))/(1-xn_tm1(2,1)*K);u(1,1)*sin(xn_tm1(3,1));u(1,1)*(sin(u(2,1))/L-K*cos(xn_tm1(3,1))/(1-xn_tm1(2,1)*K)) ];
                   t>=0];
end
% System Dynamics
for i = 1:N
    K=2/R/pi*atan(100-.5*x(1,i)^2);
    Constraints=[Constraints;
                 x(:,i+1) == x(:,i)+dt*[u(1,i+1)*cos(x(3,i))/(1-x(2,i)*K);u(1,i+1)*sin(x(3,i));u(1,i+1)*(sin(u(2,i+1))/L-K*cos(x(3,i))/(1-x(2,i)*K)) ];
                 Hx * x(:,i+1) <= bx;
                 Hu * u(:,i+1) <= bu];
end

% Terminal Constraint: enforce predicted state in the convex safe set
K=2/R/pi*atan(100-.5*x(1,N+1)^2);
Constraints=[Constraints;                       % Constraint the multipliers to be positive
             x(:,N+1) == SS*lambda;% Terminal point in the convex hull
             ones(1,length(Qfun))*lambda == 1];  % Must be convex combination --> sum to 1


% ======= Cost Definition ======
% Running Cost
Cost=0;
for i=1:N        
    Cost = Cost + P*t(1,i);
    Constraints=[Constraints; s_f-x(1,i)<=t(1,i)];
end 

% Terminal cost as convex combination of stored cost-to-go
Cost = Cost + Qfun*lambda;

% Solve the FTOCP
options = sdpsettings('verbose',0, 'solver', 'bnb', 'debug', 1);

options.ipopt.tol=1e-3;

Problem = optimize(Constraints,Cost,options);
if Problem.problem==1
    display("sub-optimal")
    x=x_guess;
    u=u_guess;
    K=2/R/pi*atan(100-.5*x_guess(1,end)^2);
    u_f=u_guess(:,end);
    u_guess=[u_guess(:,2:end) u_f];
    x_guess=[x_guess(:,2:end) x_guess(:,end)+dt*[u_f(1)*cos(x_guess(3,end))/(1-x_guess(2,end)*K);u_f(1)*sin(x_guess(3,end));u_f(1)*(sin(u_f(2))/L-K*cos(x_guess(3,end))/(1-x_guess(2,end)*K))]];
    lmbd_guess=zeros(length(Qfun),1);
    t_guess=max(s_f-x_guess(1,2:end),0);
    sol_time = nan;
else
y_f=reshape(double(x(1:2,end-1:end)), [4,1]);
K=2/R/pi*atan(100-.5*y_f(1)^2);
u_f=[1/dt*norm([(1-y_f(2)*K)*(y_f(1)-y_f(3));y_f(2)-y_f(4)],2);0];
u_guess=[double(u(:,2:end)) u_f];
x_guess=[double(x(:,2:end)) double(x(:,end))+dt*[u_f(1)*cos(double(x(3,end)))/(1-y_f(2)*K);u_f(1)*sin(double(x(3,end)));u_f(1)*(sin(u_f(2))/L-K*cos(double(x(3,end)))/(1-double(x(2,end))*K))]];
lmbd_guess=[0;double(lambda(2:end))];
t_guess=max(s_f-x_guess(1,2:end),0);

Objective = double(Cost);
sol_time = Problem.solvertime;
end


end

