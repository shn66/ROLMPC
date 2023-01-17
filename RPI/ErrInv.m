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
Ts = 0.2; % Sample time
% [Ad, Bd] = c2d(Ac,Bc,Ts);
v_bar=15;
L=2;
Ad=[1 0 0;
    0 1 Ts*v_bar;
    0 0 1];

Bd=Ts*[1 0;
       0 0;
       0 v_bar/L];
  
%% Define LQR cost matrices

Q = diag([100 100 500]);
R = diag([10 10]);

%% Compute the state feedback gain F for LQR controller
[K, P] = dlqr(Ad,Bd,Q,R);
F = -K;
%% Constraints on input
u1max  = 1;
u2max  = pi/10;
u1min  = -1;
u2min  = -pi/10;

%% Constraints on states
x1max = 1;
x2max = 1;
x3max = pi/8;
x1min = 0;
x2min = -1;
x3min = -pi/8;

%% State constraint set
Hx = [1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1];
Kx = [x1max -x1min x2max -x2min x3max -x3min]';

%% Input constraint set
Hu = [1 0;-1 0;0 1;0 -1];
Ku = [u1max -u1min u2max -u2min]';

%% Disturbance set
Hd = [1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1];
Kd = [Ts/20 Ts/20 Ts^2*v_bar/20*3 Ts^2*v_bar/20*3  Ts*v_bar/20 Ts*v_bar/20]';
Kd =Kd+d;
%% Constraint set of S, W
fA = Ad+Bd*F;
Hs = [Hx;Hu*F];
Ks = [Kx;Ku];
S = Polyhedron('A',Hs,'B',Ks);
W = Polyhedron('A',1*Hd,'B',real(1*Kd));

%% Compute mRPI


OmegaK =W; % Set the initial Omega
maxiter = 100; % Set the maximum iterations  
for i = 1:maxiter 

    Reach = fA*OmegaK;
    OmegaKplusOne = Reach+W;
    if OmegaKplusOne == OmegaK
        fprintf(['Completed at iteration ' num2str(i) ' \n'])
        break 
    end
    % the OmegaKplusOne we've got in this iteration will be the OmegaK for 
    % the next iteration
    OmegaK = OmegaKplusOne;
end
mRPI = OmegaKplusOne;

%% Plot Oinf
% figure
% plot(mRPI, 'color', 'r', 'alpha', 0.2, W, 'color', 'b', 'alpha', 0.8, mRPI.outerApprox(), 'color', 'g', 'alpha', 0.01)
% title('mRPI')
% axis equal
E=mRPI;
%% Restore warning state
warning(Wstate);

% end