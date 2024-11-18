function [S,c] = D_set1(d,z,L, gamm, UBz, LBz)
% z contains (x,u) data

T=size(d,2);
nd=size(d,1);
nz=size(z,1);
S=sdpvar(nd, nd);
c=sdpvar(nd,1);

[Pz,Psigz,~]=svd(z(:,1:T)-1*mean(z(:,1:T)')');
% Pz=Pz*Psigz(:,1:nz);
% Pz=eye(nz);
[Pd,Psigd,~]=svd(d-1*mean(d')');
% Pd=Pd*Psigd(:,1:nd);
% Pd=eye(nd);
Q=blkdiag(8*gamm^2, 2*L^2*eye(nz),-eye(nd));
% Q=blkdiag(8*gamm^2, 2*L^2*Pz'*Pz,-Pd'*Pd);


B_z=sdpvar(nz, nz, 'full');

Mz=[-2*LBz'*B_z*UBz UBz'*B_z'*Pz+LBz'*B_z*Pz zeros(1,nd);
     Pz'*B_z*UBz+Pz'*B_z'*LBz -Pz'*B_z*Pz-Pz'*B_z'*Pz zeros(nz, nd);
     zeros(nd, nz+nd+1)];

% UBz=Pz\UBz;
% LBz=Pz\LBz;
 
M11=[1;zeros(nz+nd,1)]*[1;zeros(nz+nd,1)]'-Mz;
t_d=sdpvar(T,1);

for i=1:T
    Md=(blkdiag(1,Pz,Pd)-[0;z(:,i);d(:,i)]*[1 zeros(1,nz+nd)])'*Q*(blkdiag(1,Pz,Pd)-[0;z(:,i);d(:,i)]*[1 zeros(1,nz+nd)]);
%     Md=(eye(1+nd+nz)-blkdiag(1,Pz,Pd)\[0;x(:,i);d(:,i)]*[1 zeros(1,nz+nd)])'*Q*(eye(1+nd+nz)-blkdiag(1,Pz,Pd)\[0;x(:,i);d(:,i)]*[1 zeros(1,nz+nd)]);
    M11=M11-t_d(i)*Md; 
end

% LMI=[M11 [-c';zeros(nz,nd);Pd'];
%      [-c';zeros(nz,nd);Pd']' S];
% 
LMI=[M11 [-c';zeros(nz,nd);eye(nd)];
     [-c';zeros(nz,nd);eye(nd)]' S];

constr= [ S>=0.000*eye(nd); LMI>=0; B_z>=0;t_d>=0];
% ops = sdpsettings('solver', 'scs','verbose',1);
ops = sdpsettings('solver', 'sdpt3','verbose',0);
% ops.scs.rho_x=0.09;
% ops.scs.alpha=.05;
% ops.scs.rho_x=0.001;
% ops.scs.alpha=0.1;
% ops.scs.max_iters=20000;

% ops.sedumi.beta=0.8;
% ops.sedumi.theta=0.1;
% ops.sedumi.beta=0.8;
% ops.sedumi.theta=0.5;

% sol=optimize(constr, 10000*trace(S), ops);
sol=optimize(constr, 1000*trace(S), ops);
S=abs(inv(Pd)*value(S)*inv(Pd)');

c=Pd\value(c);



end