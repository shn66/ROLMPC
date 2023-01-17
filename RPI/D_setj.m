function [S,c] = D_setj(d,x,L, gamm, UBx, LBx, Sj,cj)

T=min(size(d,2),800);
nx=size(x,1);
nd=size(d,1);
S=sdpvar(nd, nd);
c=sdpvar(nd,1);
[Px,Psigx,~]=svd(x(:,1:T)-1*mean(x(:,1:T)')');
% Px=Px*Psigx(:,1:nx);
% Px=eye(nx);
[Pd,Psigd,~]=svd(d-1*mean(d')');
% [Pd,Psigd,~]=svd(d);
% Pd=Pd*Psigd(:,1:nd);
% Pd=eye(nx);
Q=blkdiag(8*gamm^2, 2*L^2*eye(nx),-eye(nd));

B_x=sdpvar(nx, nx, 'full');

Mx=[-2*LBx'*B_x*UBx UBx'*B_x'*Px+LBx'*B_x*Px zeros(1,nd);
     Px'*B_x*UBx+Px'*B_x'*LBx -Px'*B_x*Px-Px'*B_x'*Px zeros(nx,nd);
     zeros(nd, nx+nd+1)];


 
M11=[1;zeros(nd+nx,1)]*[1;zeros(nd+nx,1)]'-Mx;
iSj=inv(Sj);
% Mj_cont=[iSj -iSj*cj eye(nx);-cj'*iSj cj'*iSj*cj -cj'; eye(nx) -cj Sj];
% t_j1=sdpvar(1);
% Mj_cont=[1-cj'*iSj*cj zeros(1,nx) cj'*iSj*Pd; 
%           zeros(nx, nd+nx+1);
%           Pd'*iSj*cj  zeros(nd,nx) -Pd'*iSj*Pd];

M11=[1;zeros(nd+nx,1)]*[1;zeros(nd+nx,1)]'-Mx; %-0*t_j1*Mj_cont;
t_d=sdpvar(T,1);
for i=1:T
    Md=(blkdiag(1,Px,Pd)-[0;x(:,i);d(:,i)]*[1 zeros(1,nx+nd)])'*Q*(blkdiag(1,Px,Pd)-[0;x(:,i);d(:,i)]*[1 zeros(1,nd+nx)]);
    M11=M11-t_d(i)*Md; 
end

LMI_1=[M11 [-c';zeros(nx,nd);eye(nd)];
     [-c';zeros(nx,nd);eye(nd)]' S];
LMI_2=(Sj-Pd*S*Pd');
% t_j2=sdpvar(1);

constr= [ S>=0; LMI_1>=0; LMI_2>=-0.001*eye(nd); B_x>=0;t_d>=0];
% ops = sdpsettings('solver', 'scs','verbose',0);
ops = sdpsettings('verbose',0);

% ops.scs.rho_x=0.01;
% ops.scs.alpha=.1;
ops.scs.max_iters=20000;

% ops.scs.rho_x=0.01;
% ops.scs.alpha=0.09;
% ops.scs.max_iters=15000;

ops.sedumi.beta=0.8;
ops.sedumi.theta=0.1;

sol=optimize(constr, 100*trace(S), ops);
if isnan(value(S)) || sol.problem==1
   S=0.9*Sj;
   c=0.99*cj;
else
S=abs(Pd*value(S)*Pd')+0.0000001;
% S=inv(Pd)*value(S)*inv(Pd)'+0.0000000*eye(nd);
c=min(Pd*value(c),2);
ct=0.5*(min(sqrt(S)+c,sqrt(Sj)+cj)+max(-sqrt(S)+c,-sqrt(Sj)+cj));
S=(0.5*(min(sqrt(S)+c,sqrt(Sj)+cj)-max(-sqrt(S)+c,-sqrt(Sj)+cj)))^2;
c=ct;
end


end