function [S,c] = D_set1(d,x,L, gamm, UBx, LBx)
T=size(d,2);
nd=size(d,1);
nx=size(x,1);
S=sdpvar(nd, nd);
c=sdpvar(nd,1);

[Px,Psigx,~]=svd(x(:,1:T)-1*mean(x(:,1:T)')');
% [Px,Psigx,~]=svd(x(:,1:T));
% Px=Px*Psigx(:,1:nx);
% Px=eye(nx);
[Pd,Psigd,~]=svd(d-1*mean(d')');
% [Pd,Psigd,~]=svd(d);
% Pd=Pd*Psigd(:,1:nd);
% Pd=eye(nd);
Q=blkdiag(8*gamm^2, 2*L^2*eye(nx),-eye(nd));
% Q=blkdiag(8*gamm^2, 2*L^2*Px'*Px,-Pd'*Pd);


B_x=sdpvar(nx, nx, 'full');

Mx=[-2*LBx'*B_x*UBx UBx'*B_x'*Px+LBx'*B_x*Px zeros(1,nd);
     Px'*B_x*UBx+Px'*B_x'*LBx -Px'*B_x*Px-Px'*B_x'*Px zeros(nx, nd);
     zeros(nd, nx+nd+1)];
% UBx=Px\UBx;
% LBx=Px\LBx;
% Mx=[-2*LBx'*B_x*UBx UBx'*B_x'+LBx'*B_x zeros(1,nd);
%      B_x*UBx+B_x'*LBx -B_x-B_x' zeros(nx, nd);
%      zeros(nd, nx+nd+1)];

 
M11=[1;zeros(nx+nd,1)]*[1;zeros(nx+nd,1)]'-Mx;
t_d=sdpvar(T,1);

for i=1:T
    Md=(blkdiag(1,Px,Pd)-[0;x(:,i);d(:,i)]*[1 zeros(1,nx+nd)])'*Q*(blkdiag(1,Px,Pd)-[0;x(:,i);d(:,i)]*[1 zeros(1,nx+nd)]);
%     Md=(eye(1+nd+nx)-blkdiag(1,Px,Pd)\[0;x(:,i);d(:,i)]*[1 zeros(1,nx+nd)])'*Q*(eye(1+nd+nx)-blkdiag(1,Px,Pd)\[0;x(:,i);d(:,i)]*[1 zeros(1,nx+nd)]);
    M11=M11-t_d(i)*Md; 
end

% LMI=[M11 [-c';zeros(nx,nd);Pd'];
%      [-c';zeros(nx,nd);Pd']' S];
% 
LMI=[M11 [-c';zeros(nx,nd);eye(nd)];
     [-c';zeros(nx,nd);eye(nd)]' S];
%  
% constr= [ S>=0; LMI>=0; Px'*B_x*Px-0.0000*ones(nx)>=0;t_d-0.0000*ones(T,1)>=0];
constr= [ S>=0.000*eye(nd); LMI>=0; B_x>=0;t_d>=0];
% ops = sdpsettings('solver', 'scs','verbose',1);
ops = sdpsettings('verbose',0);
% ops.scs.rho_x=0.09;
% ops.scs.alpha=.05;
ops.scs.rho_x=0.001;
ops.scs.alpha=0.1;
ops.scs.max_iters=20000;

% ops.sedumi.beta=0.8;
% ops.sedumi.theta=0.1;
% ops.sedumi.beta=0.8;
% ops.sedumi.theta=0.5;

% sol=optimize(constr, 10000*trace(S), ops);
sol=optimize(constr, 1000*S(1,1), ops);
S=abs(Pd*value(S)*Pd');
% S=value(S)+0.0000000*eye(nd);
c=Pd*value(c);

% value(B_x)
% value(t_d)


end