function [E, W, F]=ErrInv(d)
%% Computing positive invariant set 
% In this example, we compute the positive invariant set for the inverted 
% pendulum system which is formulated into an autonomous system 
% x(t+1) = (Ad + Bd*F)*x(t). We use MPT command invariantSet() to compute 
% the positive invariant set directly, also we compute it manually with the 
% algorithm which was introduced in the lecture.

%% Fix obselete use of NARGCHK
Wstate = warning('off');


%% Transform the continuous-time state-space model into discretized state-space model
Ts = 0.1; 
L=2.0;

v_bar=15;
v_max=18;
epsi_max=pi/8;
ey_max=5;
del_max=pi/8;
curv=0;
curv_bnd=0;



Ad=[1 Ts*curv 0;
    0 1 Ts*v_bar;
    0 -Ts*v_bar*curv 1];

Bd=Ts*[1 0;
       0 0;
      -curv v_bar/L];
  
%% Define LQR cost matrices

% Q = diag([100 10 200]);
% R = diag([1 0.1]);

% Q = diag([10 200 100]);
% R = diag([.1 10]);

% Q = diag([100 500 500]);
% R = diag([100 200]);
Q = diag([100 500 500]);
R = diag([100 200]);

%% Compute the state feedback gain F for LQR controller
[K, P] = dlqr(Ad,Bd,Q,R);
F = -K;
% [K1,P1]=dlqr(Ad(1,1), Bd(1,1), Q(1,1), R(1,1));
% [K2, P2]=dlqr(Ad(2:end,2:end), Bd(2:end,2), Q(2:end,2:end), R(2,2));
% F=[-K1 0 0;0 -K2];

%% Disturbance set

Hd = [1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1];
ds_max= max(-2.,0.5*(2*curv_bnd^2*ey_max^2+max(-2*curv_bnd*(v_max-v_bar),2*curv_bnd*(v_max-v_bar)))*Ts);
ds_min= min(2.,0.5*(-epsi_max^2*v_bar+min(-2*curv_bnd*(v_max-v_bar),2*curv_bnd*(v_max-v_bar)))*Ts);

dey_max=0.5*(v_max-v_bar)*epsi_max*Ts;
dey_min=-0.5*(v_max-v_bar)*epsi_max*Ts;

depsi_max= 0.5*(2*abs(curv_bnd)^3*ey_max^2*v_bar+curv_bnd*epsi_max^2*v_bar+2/L*(v_max-v_bar)*del_max...
            +2*curv_bnd^2*(v_max-v_bar)*ey_max)*Ts;
depsi_min= 0.5*(-2*abs(curv_bnd)^3*ey_max^2*v_bar-curv_bnd*epsi_max^2*v_bar-2/L*(v_max-v_bar)*del_max-2*curv_bnd^2*(v_max-v_bar)*ey_max)*Ts;

Kd = [ds_max -ds_min dey_max -dey_min  depsi_max -depsi_min]';
% Kd = [Ts/20 Ts/20 Ts^2*v_bar/20*3 Ts^2*v_bar/20*3  Ts*v_bar/20 Ts*v_bar/20]';
Kd =Kd+d;
% Kd= max(Kd, 0.001);
%% Constraint set of S, W
fA = Ad+Bd*F;
% Hs = [Hx;Hu*F];
% Ks = [Kx;Ku];
% S = Polyhedron('A',Hs,'B',Ks);
W = Polyhedron('A',1*Hd,'B',real(1*Kd));

%% Compute mRPI


OmegaK =W; % Set the initial Omega
maxiter = 100; % Set the maximum iterations  

% OmegaKs =Ws;
% OmegaKl=Wl;
for i = 1:maxiter 

    Reach = fA*OmegaK;
    OmegaKplusOne = Reach+W;
     if OmegaK==OmegaKplusOne 
            fprintf(['Completed at iteration ' num2str(i) ' \n'])
            break 
     end
%     
    % the OmegaKplusOne we've got in this iteration will be the OmegaK for 
    % the next iteration
    OmegaK = OmegaKplusOne;

%     fprintf([ num2str(i) '\n'])



end
mRPI = OmegaKplusOne;

%% Plot Oinf
% figure
% plot(mRPI, 'color', 'r', 'alpha', 0.2, W, 'color', 'b', 'alpha', 0.8, mRPI.outerApprox(), 'color', 'g', 'alpha', 0.01)
% title('mRPI')
% axis equal
E=mRPI.minHRep();
%% Restore warning state
warning(Wstate);

% end