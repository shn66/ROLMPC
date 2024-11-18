function [ x_LMPC, u_LMPC, y_cl_out, x_cl_out, u_cl_out, d_cl_out, xn_cl_out, un_cl_out, x_clg_out, Xt_out, Ut_out, W_out, x_data, u_data, d_data, IterationCost_out, sol_times_out]...
 = RLMPC_bi(x0, y_cl, x_cl, u_cl,d_cl, xn_cl, un_cl, x_clg, xcg, xc, IterationCost, Lip, gamm, C,P, R, N, Ny, Iterations, X, U, x_data,u_data, d_data, LMPC_options)
%% Learning Model Predictive Control
% This function runs the LMPC for a user-defined number of iterations.
% The input to this functions are:
% - x0: initiol condition used at each iteration
% - SS: the initial safe set constructed using a feasible trajectory
% - (A,B): matrices defining the system dynamics
% - N: LMPC horizon length
% - Iterations: max number of iterations
% - X: polyhedron representng the state constraints
% - U: polyhedron representng the input constraints
% - IterationCost: vector collecting the iteration cost
% - Solver: solver to use for the FTOCP
ny=size(C,1);

% Now start iteration loop
j = 1;
SS = xn_cl{1}(:,1:end-1);
% SS=y_cl{1};
ss{1}=SS;
Qfun= IterationCost{1};
qfun{1}=Qfun;
IterationCost_out{1} = IterationCost{1};
s_f=40;
eps=0.01;
x_cl_out{1}=x_cl{1};
u_cl_out{1}=u_cl{1};
d_cl_out{1}=d_cl{1};
y_cl_out{1}=y_cl{1};
xn_cl_out{1}=xn_cl{1};
un_cl_out{1}=un_cl{1};
x_clg_out{1}=x_clg{1};
TrajCost=[];

while (j <= Iterations)
    t = 1;        % Initialize time
    x_LMPC = x0;  % Initial condition
    x_LMPCg = x0;
    xn_LMPC = x0;
    d=[];
    sol_times = [];
    UBx=[X.b-10;X.b(3);X.b(5)];
    LBx=[-X.b(2);-X.b(4);-X.b(6)];
    UBu=[max(U.V(:,1)); max(U.V(:,2))];
    LBu=[min(U.V(:,1)); min(U.V(:,2))];
    if j==1
        [S1,c1]=D_set1(d_data(1,:),[x_data(1,:);u_data(1,:)],Lip(1), gamm(1), [UBx(1);UBu(1)], [LBx(1);LBu(1)]);
        [S2,c2]=D_set1(d_data(2,:),[x_data(2,:);u_data(1,:)],Lip(2), gamm(2), [UBx(2);UBu(1)], [LBx(2);LBu(1)]);
        [S3,c3]=D_set1(d_data(3,:),[x_data(3,:);u_data(:,:)],Lip(3), gamm(3), [UBx(3);UBu], [LBx(3);LBu]);

        S=diag([S1, S2, S3]);
        c=[c1;c2;c3];
    else
        d_data=[d_data d_cl_out{j}];
        x_data=[x_data x_cl_out{j}(:,1:end-1)];
        u_data=[u_data u_cl_out{j}];
        
        [S1,c1]=D_setj(d_data(1,:),[x_data(1,:);u_data(1,:)],Lip(1), gamm(1), [UBx(1);UBu(1)], [LBx(1);LBu(1)], S1,c1);
        [S2,c2]=D_setj(d_data(2,:),[x_data(2,:);u_data(1,:)],Lip(2), gamm(2), [UBx(2);UBu(1)], [LBx(2);LBu(1)], S2, c2);
        [S3,c3]=D_setj(d_data(3,:),[x_data(3,:);u_data(:,:)],Lip(3), gamm(3), [UBx(3);UBu], [LBx(3);LBu], S3, c3);

        S=diag([S1, S2, S3]);
        c=[c1;c2;c3];
    end
    d_b=[sqrt(S(1,1))+c(1);sqrt(S(1,1))-c(1);sqrt(S(2,2))+c(2);sqrt(S(2,2))-c(2);sqrt(S(3,3))+c(3);sqrt(S(3,3))-c(3)];
    [E, W, F]=ErrInv(d_b);
    
    Xt=X-E.outerApprox(); Xt.minHRep();
    Ut=F*E;
    Ut=U-Ut.outerApprox(); Ut.minHRep();
    
    Xt_out{j}=Xt; Ut_out{j}=Ut; W_out{j}=W;    
    flag=0;
    SSQfun = Polyhedron([SS', Qfun']);
%     SSQfun.computeHRep()
    SSQfun.computeVRep();
%     SS = SSQfun.V(:,1:2)';
    SS=SSQfun.V(:,1:end-1)';
    Qfun = SSQfun.V(:, end)';
%     [Qfun,I]=sort(Qfun, 'Descend');
%     SS=SS(:,I);
    
    % This is the time loop for the j-th iteration. This loop terminates when the
    % closed-loop trajectory reaches the terminal point x_F=0.
    TrajCost=[TrajCost qfun{j}(1)];
    [~,I]=sort(TrajCost);
    SS_used=[ss{1}];
    Qfun_used=[qfun{1}];
    for k=1:j
        if I(k)~=1
            SS_used=[SS_used ss{I(k)}];
            Qfun_used=[Qfun_used qfun{I(k)}];
        end
    end
    
    u_guess=[[0;0] un_cl_out{1}(:,1:N)];
    x_guess=xn_cl_out{1}(:,1:N+1);
    lmbd_guess=zeros(length(Qfun_used),1);
    t_guess=zeros(1,N);
    while flag==0       
        x_frenet = global2frenet(x_LMPCg(:,t), xc, xcg);
        s = x_frenet(1); ey = x_frenet(2); ephi = x_frenet(3);
        d=[d [s;ey;ephi]-x_LMPC(:,t)];
        x_LMPC(:,t)=[s;ey;ephi];
        disp(['Time step: ', num2str(t), ' s: ', num2str(x_LMPC(1,t)),' ey: ', num2str(x_LMPC(2,t)),' epsi: ', num2str(x_LMPC(3,t)),' Iteration: ', num2str(j), ' Best Cost: ',num2str(IterationCost_out{j}(1))])
        % Solve the LMPC at time t of the j-th iteration
        if norm(x_LMPC(1,t)-s_f-29)<=1 || x_LMPC(1,t)>=s_f+28 || t>=4 
                flag=1;
        end
        if t==1
            xntm1=xn_LMPC(:,t);
        else
            xntm1=xn_LMPC(:,t-1);
        end
        
        if x_LMPC(1,t)>=s_f+1
            R=20;L=2;
            K=2/R/pi*atan(100-.5*x_LMPC(1,t)^2);
            un_LMPC(:,t) = [12;atan(L*K-.09*x_LMPC(2,t)-.1*x_LMPC(3,t))];
            xPred=[xn_LMPC(:,t) xn_LMPC(:,t-1)];
            
        else

                                             
            [xPred, uPred, x_guess, u_guess, lmbd_guess, t_guess, sol_time ] = mipFTOCP_bi(x_LMPC(:,t),xntm1, t, N, Ny, C, P, Qfun_used, SS_used,...
                                                 s_f, Xt, Ut, E, x_guess, u_guess, lmbd_guess, t_guess, LMPC_options.solver);
            un_LMPC(:,t) = double(uPred(:,2));
            
            sol_times = [sol_times sol_time];
            
        end
        u_LMPC(:,t) = un_LMPC(:,t)+F*(x_LMPC(:,t)-xn_LMPC(:,t));
        xn_LMPC(:,t)=double(xPred(:,1));
        % Update system position
        
        
        x_LMPC(:,t+1)=state_update_bicycle(x_LMPC(:,t),u_LMPC(:,t));
        xn_LMPC(:,t+1)=state_update_bicycle(xn_LMPC(:,t),un_LMPC(:,t));
        x_LMPCg(:,t+1)=state_update_bicycle_global(x_LMPCg(:,t),u_LMPC(:,t));
        t = t + 1;
       
     end

    % Now save the data, update cost and safe set.
    x_cl_out{j+1} = smoother(x_LMPC);
    u_cl_out{j+1} = smoother(u_LMPC);
    d_cl_out{j+1}= d;
    xn_cl_out{j+1}=xn_LMPC;
    un_cl_out{j+1}=un_LMPC;
    x_clg_out{j+1}=smoother(x_LMPCg);
    
    
    IterationCost_out{j+1} = ComputeCost_bi(xn_LMPC, un_LMPC,P);
    SS   = [SS, xn_LMPC(:,1:end-1)];
    ss{j+1}=xn_LMPC(:,1:end-1);
    Qfun = [Qfun, IterationCost_out{j+1}];
    qfun{j+1}=IterationCost_out{j+1};
    
    sol_times_out{j} = [nanmin(sol_times) nanmax(sol_times)];
    
    % increase Iteration index and restart
    j = j + 1;
    if j <= Iterations
        clear x_LMPC
        clear x_LMPCg
        clear u_LMPC
        clear xn_LMPC
        clear un_LMPC
        clear sol_times
    end
    
end
end