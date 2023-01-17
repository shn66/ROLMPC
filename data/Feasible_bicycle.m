% clear all
% clc
x_0=[0;0;0];
x_0g=[0;0;0];

x=[];
x=[x x_0+[0;1;0]];
xb=x;
xg=[];
xg=[xg x_0g+[0;1;0]];

L=2;
R=20;
K_est=1/R;


xc=[];
xc=[xc x_0];
xcg=[];
xcg=[xcg x_0g];
xcgu=xcg+[0;3.5;0];
xcgd=xcg-[0;3.5;0];
i=1;
while xc(1,i)<70
    K=2/R/pi*(atan(100-.5*xc(1,i)^2));
    Uc=[1;atan(K*L)];
    x_next=state_update_bicycle(xc(:,i),Uc);
    x_nextg=state_update_bicycle_global(xcg(:,i),Uc);
    x_nextgu=x_nextg+[cos(x_nextg(3)+pi/2)*4;sin(x_nextg(3)+pi/2)*4;x_nextg(3)];
    x_nextgd=x_nextg+[cos(x_nextg(3)-pi/2)*4;sin(x_nextg(3)-pi/2)*4;x_nextg(3)];
    xc=[xc x_next];
%     x_next
    xcg=[xcg x_nextg];
    xcgu=[xcgu x_nextgu];
    xcgd=[xcgd x_nextgd];
    i=i+1;
end
u_feasible=[];
d=[];
i=1;
while x(1,i)<70
    K=2/R/pi*(atan(100-.5*x(1,i)^2));
    [ind,dist]=dsearchn(xcg(1:2,:)', xg(1:2,i)');
    s=xc(1,ind);
    ephi=xg(3,i)-xcg(3,ind);
    if ind~=1
       diff=xcg(1:2,ind)-xcg(1:2,ind-1);
       diff_e=xg(1:2,i)-xcg(1:2,ind);
    else
       diff=xcg(1:2,ind+1)-xcg(1:2,ind);
       diff_e=xg(1:2,i)-xcg(1:2,ind);
    end
    ey=-dist*sign(det([diff_e diff]));
    d=[d [s;ey;ephi]-x(:,i)];
    x(:,i)=[s;ey;ephi];

    U=[8;atan(L*K-.09*x(2,i)-.1*x(3,i))];
    x_next=state_update_bicycle(x(:,i),U);

    xb_next=state_update_bicycle(xb(:,i),U);
    x_nextg=state_update_bicycle_global(xg(:,i),U);

    x=[x x_next];
    xb=[xb xb_next];
    xg=[xg x_nextg];
    u_feasible=[u_feasible U];
    i=i+1;
end

x_feasible=x;
xn_feasible=xb;

