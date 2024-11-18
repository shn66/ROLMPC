dt=0.1;
f=@(x,u) [x(1)+dt*(x(1)-x(2)^3);x(2)+dt*(x(1)^2-u)];
t_max=80;
x_traj_1=[0.4;0.4];
u_traj_1 = [];
t=size(x_traj_1,1);

while t<=t_max
    u=x_traj_1(1,end)^2-6*(x_traj_1(1,end)-1)+4*(x_traj_1(2,end)-1);
    x_traj_1=[x_traj_1 f(x_traj_1(:,end),u)];
    u_traj_1=[u_traj_1 u];
    t=t+1;
end
y_1=[x_traj_1(1,1:end-1);x_traj_1(1,2:end)];
Y1=Polyhedron('V', y_1');


x_traj_0=[0.4;0.4];
u_traj_0=[];
t=size(x_traj_0,1);

while t<=t_max
    u=x_traj_0(1,end)^2-6*(x_traj_0(1,end)-0)+4*(x_traj_0(2,end)-0);
    x_traj_0=[x_traj_0 f(x_traj_0(:,end),u)];
    u_traj_0=[u_traj_0 u];
    t=t+1;
end

y_0=[x_traj_0(1,1:end-1);x_traj_0(1,2:end)];
Y0=Polyhedron('V', y_0');

Y01=Polyhedron('V', [Y0.V;Y1.V]);
subplot(2,1,1);
plot(Y01,'Color','b', 'alpha', 0.5, 'DisplayName', 'CS_y');
Y_v= Y01.V';
hold on
plot(Y_v(1,:),Y_v(2,:),'ob', 'DisplayName', 'Y_{data}');

xlabel("y_{t}")
ylabel("y_{t+1}")
legend("Location","southeast");

X_ss=[];
Fx= @(y)[y(1);((y(1)*(1+dt)-y(2))/dt)^(1/3)];

for i=1:size(Y_v,2)
    X_ss=[X_ss Fx(Y_v(:,i))];
end

y_e1=[0.5673;0.5584];
y_e2=[1.0642;1.0202];

lmbd=rand(1,10);

Y_extra=y_e1*lmbd+y_e2*(1-lmbd);

for i=1:10
    X_ss=[X_ss Fx(Y_extra(:,i))];
end

y_e1=[0.4;0.4336];
y_e2=[.5435;.5383];

lmbd=rand(1,10);

Y_extra=y_e1*lmbd+y_e2*(1-lmbd);

for i=1:10
    X_ss=[X_ss Fx(Y_extra(:,i))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2);
Inv = polytope_approx(X_ss(1,:), X_ss(2,:));
plot(Inv, 'alpha',0.5)
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlabel("x[0]")
ylabel("x[1]")
legend("Location","southeast");

plot(x_traj_0(1,:), x_traj_0(2,:),'ok','HandleVisibility','off');
plot(x_traj_1(1,:), x_traj_1(2,:), 'ok', 'DisplayName', 'x_{data}');
x_traj_2=[0.4335;0.5729];
t=size(x_traj_2,1);

while t<=t_max  
    u = data_based_control(x_traj_2(:,end), Fx, [y_0 y_1], [u_traj_0 u_traj_1]);
    x_traj_2=[x_traj_2 f(x_traj_2(:,end),u)];
    t=t+1;
end

plot(x_traj_2(1,:), x_traj_2(2,:),'+g','DisplayName', 'x_{test}');


x_traj_2=[0.5027;0.7631];
t=size(x_traj_2,1);

while t<=t_max
    u = data_based_control(x_traj_2(:,end), Fx, [y_0 y_1], [u_traj_0 u_traj_1]);
    x_traj_2=[x_traj_2 f(x_traj_2(:,end),u)];
    t=t+1;
end

plot(x_traj_2(1,:), x_traj_2(2,:),'+g','HandleVisibility', 'off');




