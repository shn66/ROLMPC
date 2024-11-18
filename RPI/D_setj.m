function [S,c] = D_setj(d,z,L, gamm, UBz, LBz, Sj,cj)

% z contains (x,u) data

T=min(size(d,2),800);
nz=size(z,1);
nd=size(d,1);
S=sdpvar(nd, nd);
c=sdpvar(nd,1);
[Pz,Psigz,~]=svd(z(:,1:T)-1*mean(z(:,1:T)')');
% Pz=Pz*Psigz(:,1:nz);
% Pz=eye(nz);
[Pd,Psigd,~]=svd(d-1*mean(d')');
% [Pd,Psigd,~]=svd(d);
% Pd=Pd*Psigd(:,1:nd);
% Pd=eye(nd);
Q=blkdiag(8*gamm^2, 2*L^2*eye(nz),-eye(nd));

B_z=sdpvar(nz, nz, 'full');

Mz=[-2*LBz'*B_z*UBz UBz'*B_z'*Pz+LBz'*B_z*Pz zeros(1,nd);
     Pz'*B_z*UBz+Pz'*B_z'*LBz -Pz'*B_z*Pz-Pz'*B_z'*Pz zeros(nz,nd);
     zeros(nd, nz+nd+1)];

 
M11=[1;zeros(nd+nz,1)]*[1;zeros(nd+nz,1)]'-Mz;
iSj=inv(Sj);


M11=[1;zeros(nd+nz,1)]*[1;zeros(nd+nz,1)]'-Mz;
t_d=sdpvar(T,1);
for i=1:T
    Md=(blkdiag(1,Pz,Pd)-[0;z(:,i);d(:,i)]*[1 zeros(1,nz+nd)])'*Q*(blkdiag(1,Pz,Pd)-[0;z(:,i);d(:,i)]*[1 zeros(1,nd+nz)]);
    M11=M11-t_d(i)*Md; 
end

LMI_1=[M11 [-c';zeros(nz,nd);eye(nd)];
     [-c';zeros(nz,nd);eye(nd)]' S];
LMI_2=(Sj-Pd*S*Pd');


constr= [ S>=0; LMI_1>=0; LMI_2>=0.00001*eye(nd); B_z>=0;t_d>=0];
ops = sdpsettings('solver', 'sdpt3','verbose',0);
% ops = sdpsettings('verbose',0);

ops.scs.rho_x=0.01;
ops.scs.alpha=.1;
ops.scs.max_iters=20000;

% ops.scs.rho_x=0.01;
% ops.scs.alpha=0.09;
% ops.scs.max_iters=15000;

% ops.sedumi.beta=0.8;
% ops.sedumi.theta=0.1;

sol=optimize(constr, 1000*trace(S), ops);
if isnan(value(S)) || sol.problem==1
   S=0.99*Sj;
   c=cj;
else
S=abs(inv(Pd)*value(S)*inv(Pd'));%+0.0000001;
% S=inv(Pd)*value(S)*inv(Pd)'+0.0000000*eye(nd);
c=inv(Pd)*value(c);
ct=0.5*(min(sqrt(S)+c,sqrt(Sj)+cj)+max(-sqrt(S)+c,-sqrt(Sj)+cj));
S=(0.5*(min(sqrt(S)+c,sqrt(Sj)+cj)-max(-sqrt(S)+c,-sqrt(Sj)+cj)))^2;
c=ct;
end


end