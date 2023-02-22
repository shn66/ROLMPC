clc
clear all
close all


%% Options
LMPC_options.solver = 'ipopt'; % Options are 'gurobi' or 'quadprog'. IMPORTANT: Use gurobi for better precision;

%%
% Define your system constraints and cost 

[P, C, U, X, Q, R, N, Ny] = DefineSystem_bicycle();


%% Load the first feasible solution

load('./data/bicycle_iter0.mat');
load('./data/lip.mat');
y_feasible=C*xn_feasible;
un_feasible=u_feasible;
s_f=40;
%% Plot the first feasible solution

ind_it=dsearchn(x_feasible(1,:)',s_f);
ind_g=dsearchn(xc(1,:)',s_f);
plot( xcg(1,1:ind_g), xcg(2,1:ind_g),'k-')
hold on
plot( xcgu(1,1:ind_g), xcgu(2,1:ind_g),'r-')
plot( xcgd(1,1:ind_g), xcgd(2,1:ind_g),'r-')
plot( xg(1,1:ind_it), xg(2,1:ind_it),'b-', 'LineWidth', 2.0)


%% ========================================================================
%  ======================= Now Run Learning MPC ===========================
%  ========================================================================
%% Initialize Safe Set and Q-funtion
x_cl{1} = x_feasible; % Vector collecting the state of the performed iterations
u_cl{1} = u_feasible; % Vector collecting the input of the performed iterations
x_clg{1} =xg;

xn_cl{1} = xn_feasible;  % Vector collecting the nominal state of the performed iterations
un_cl{1} = un_feasible;  % Vector collecting the nominal input of the performed iterations

% Modify nominal state trajectory for safe set construction
[xn,un]=make_nominal(Xt,Ut,E,x_cl{1}, xn_cl{1}, un_cl{1});

xn_cl{1}=xn;
un_cl{1}=un;

Y_flat=[];
y_feasible=C*xn;
for i=1:Ny
Y_flat=[Y_flat;y_feasible(:,i:end-Ny+i)];
end
y_cl{1}=Y_flat;
Qfun = ComputeCost_bi(xn, un, P);
IterationCost{1} = Qfun;
d_cl{1}=d;

%% Convert constructed nominal trajectory into global coordinates
xng=frenet2global(xn, xc, xcg);

%% Get Error Invariant
Lip=[l1, l2, l3];
gamm= [w1, w2, w3];
UBx=[X.b(1)-10;X.b(3);X.b(5)];
LBx=[X.b(2);-X.b(4);-X.b(6)];
d_data_1=[d d_data];
x_data_1=[x_feasible x_data];
[S1,c1]=D_set1(d_data_1(1,:),x_data_1(1,:),Lip(1), gamm(1), UBx(1), LBx(1));
[S2,c2]=D_set1(d_data_1(2,:),x_data_1(2,:),Lip(2), gamm(2), UBx(2), LBx(2));
[S3,c3]=D_set1(d_data_1(3,:),x_data_1(3,:),Lip(3), gamm(3), UBx(3), LBx(3));
S=diag([S1, S2, S3]);
c=[c1;c2;c3];
d_b=[sqrt(S(1,1))+c(1);sqrt(S(1,1))-c(1);sqrt(S(2,2))+c(2);sqrt(S(2,2))-c(2);sqrt(S(3,3))+c(3);sqrt(S(3,3))-c(3)];
[E,W,F]=ErrInv(d_b);
%% Define tightened constraints
Xt=X-E.outerApprox();
Ut=F*E;
Ut=U-Ut.outerApprox();

%%
% %% Now run the LMPC
% Pick number of iterations to run
Iterations = 30;

% Run the LMPC
x0 = x_cl{1}(:,1);

[x_LMPC, u_LMPC, y_cl, x_cl, u_cl, d_cl, xn_cl, un_cl, x_clg, Xt_it, Ut_it, W_it, x_data, d_data, IterationCost] = LMPC_bi(x0, y_cl, x_cl, u_cl, d_cl, xn_cl, un_cl, x_clg, xcg, xc, ...
                                    IterationCost, Lip, gamm, C, P, R, N, Ny, Iterations, X, U, x_data, d_data, LMPC_options);

%% Display Cost and Plot the Results
format long
clc
for i = 1:Iterations+1
    TrajCosts(i)=IterationCost{i}(1);
end

%%
SS = [];
SSn= [];
SSn_goal=[];
SSng= [];
SSng_goal=[];
SSu = [];
SSy = [];
SSg = [];
for i=1:Iterations+1
    SS=[SS, smoother(x_cl{Order(i)})];
    SSu=[SSu, smoother(u_cl{Order(i)})];
    SSg=[SSg, smoother(x_clg{Order(i)})];

    ind_ngl=dsearchn(xn_cl{Order(i)}(1,:)',s_f);
    SSn_goal=[SSn_goal, smoother(xn_cl{Order(i)}(:,ind_ngl:end)) ];
    SSn=[SSn, smoother(xn_cl{Order(i)}(:,1:ind_ngl))];
    SSng=[SSng, smoother(frenet2global(xn_cl{Order(i)}(:,1:ind_ngl), xc, xcg))];
    SSng_goal=[SSng_goal, smoother(frenet2global(xn_cl{Order(i)}(:,ind_ngl:end), xc, xcg)) ];
end

%%
figure
ind_g=dsearchn(xc(1,:)',s_f);
plot( xcg(1,1:ind_g), xcg(2,1:ind_g),'k--','DisplayName',"Centerline")
hold on
plot( xcgu(1,1:ind_g), xcgu(2,1:ind_g),'-','LineWidth', 2.0, 'color', [0.7,0.7,0.7],'DisplayName',"Boundary")
plot( xcgd(1,1:ind_g), xcgd(2,1:ind_g),'-','LineWidth', 2.0, 'color', [0.7,0.7,0.7],'HandleVisibility','off')
plot( [xcgu(1,ind_g+1),xcgd(1,ind_g+1)] , [xcgu(2,ind_g+1), xcgd(2,ind_g+1)],'-','LineWidth', 2.0, 'color', [0,0,.7],'DisplayName',"Goal")
plot( xcgu(1,ind_g+1:ind_g+30), xcgu(2,ind_g+1:ind_g+30),'-','LineWidth', 2.0, 'color', [0.,0.,0.7],'HandleVisibility','off')
plot( xcgd(1,ind_g+1:ind_g+30), xcgd(2,ind_g+1:ind_g+30),'-','LineWidth', 2.0, 'color', [0,0,0.7],'HandleVisibility','off')
plot( [xcgu(1,ind_g+30),xcgd(1,ind_g+30)] , [xcgu(2,ind_g+30), xcgd(2,ind_g+30)],':','LineWidth', 2.0, 'color', [0,0,.7],'HandleVisibility','off')

clr=[interp1([1, round(0.5*(Iterations+1)), Iterations+1], [0,1,1], [1:Iterations+1], 'cubic');
     interp1([1, round(0.5*(Iterations+1)), Iterations+1], [1,1,0], [1:Iterations+1], 'cubic');
     interp1([1, round(0.5*(Iterations+1)), Iterations+1], [0,0,0], [1:Iterations+1],'cubic')];
for i=1:2:Iterations+1
    it=i;
    ind_it=dsearchn(x_cl{it}(1,:)',s_f);
    plot( x_clg{it}(1,1:ind_it), x_clg{it}(2,1:ind_it),'-', 'LineWidth', 1.+4/3*(Iterations+1-i)/(Iterations+1), 'color',clr(:,i)', 'DisplayName',"Iter: "+num2str(Iterations+1-i));
end
axis equal
legend("Location","southeast","NumColumns",3)
xlabel("X")
ylabel("Y")
%%
% figure
% plot( xcgu(1,1:ind_g), xcgu(2,1:ind_g),'-','LineWidth', 2.0, 'color', [0.7,0.7,0.7],'DisplayName',"Boundary")
% hold on
% plot( xcgd(1,1:ind_g), xcgd(2,1:ind_g),'-','LineWidth', 2.0, 'color', [0.7,0.7,0.7],'HandleVisibility','off')
clr=[interp1([1, round(0.5*(Iterations+1)), Iterations+1], [0,1,1], [1:Iterations+1], 'cubic');
     interp1([1, round(0.5*(Iterations+1)), Iterations+1], [1,1,0], [1:Iterations+1], 'cubic');
     interp1([1, round(0.5*(Iterations+1)), Iterations+1], [0,0,0], [1:Iterations+1],'cubic')];
plot(0:25, 18*ones(1,26),'-','LineWidth', 2.0, 'color', [0.7,0.7,0.7],'DisplayName',"Constraints");
hold on
plot(0:25, 0*ones(1,26),'-','LineWidth', 2.0, 'color', [0.7,0.7,0.7],'HandleVisibility','off');

for i=1:2:Iterations+1
    it=i;
    ind_it=dsearchn(x_cl{it}(1,:)',s_f-0.5)-2;
    plot( [0:ind_it-1], smoother(u_cl{it}(1,1:ind_it)),'-', 'LineWidth', 1.+4/3*(Iterations+1-i)/(Iterations+1), 'color',clr(:,i)', 'DisplayName',"Iter: "+num2str(Iterations+1-i));
end
plot(0:25, 8*ones(1,26),'-','LineWidth', 1.0, 'color', [1.,0,0], 'DisplayName',"Iter: 0");
% axis equal
legend("Location","southeast","NumColumns",2)
xlabel("time")
ylabel("speed")

%%
clr=[interp1([1, round(0.5*(Iterations+1)), Iterations+1], [0,1,1], [1:Iterations+1], 'cubic');
     interp1([1, round(0.5*(Iterations+1)), Iterations+1], [1,1,0], [1:Iterations+1], 'cubic');
     interp1([1, round(0.5*(Iterations+1)), Iterations+1], [0,0,0], [1:Iterations+1],'cubic')];
plot(0:25, (-pi/2+pi/20)*ones(1,26),'-','LineWidth', 2.0, 'color', [0.7,0.7,0.7],'DisplayName',"Constraints");
hold on
plot(0:25, (pi/2-pi/20)*ones(1,26),'-','LineWidth', 2.0, 'color', [0.7,0.7,0.7],'HandleVisibility','off');
for i=1:2:Iterations+1
    it=i;
    ind_it=dsearchn(x_cl{it}(1,:)',s_f-1);
    plot( [0:ind_it-1], smoother(u_cl{it}(2,1:ind_it)),'-', 'LineWidth', 1.+4/3*(Iterations+1-i)/(Iterations+1), 'color',clr(:,i)', 'DisplayName',"Iter: "+num2str(Iterations+1-i));
end
% axis equal
legend("Location","southeast","NumColumns",3)
xlabel("time")
ylabel("steering")

