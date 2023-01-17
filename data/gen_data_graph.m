% clear all
% clc
L=2;
R=20;
x_0=[0;0;0];
x_0g=[0;0;0];

xc=[];
xc=[xc x_0];
xcg=[];
xcg=[xcg x_0g];
xcgu=xcg+[0;6.5;0];
xcgd=xcg-[0;6.5;0];
i=1;
while xc(1,i)<60
    K=2/R/pi*atan(100-.5*xc(1,i)^2);
    Uc=[1.;atan(K*L)];
    x_next=state_update_bicycle(xc(:,i),Uc);
    x_nextg=state_update_bicycle_global(xcg(:,i),Uc);
    x_nextgu=x_nextg+[cos(x_nextg(3)+pi/2)*6.5;sin(x_nextg(3)+pi/2)*6.5;x_nextg(3)];
    x_nextgd=x_nextg+[cos(x_nextg(3)-pi/2)*6.5;sin(x_nextg(3)-pi/2)*6.5;x_nextg(3)];
    xc=[xc x_next];
    xcg=[xcg x_nextg];
    xcgu=[xcgu x_nextgu];
    xcgd=[xcgd x_nextgd];
    i=i+1;
end
u_feasible=[];
d_data=[];
x_data=[];
iters=5;
I=0;
while I<=iters
x=[];
d=[];
x=[x x_0+[0;2*cos((iters-I)*pi/iters);0]];
xb=x;
xg=[];
xg=[xg x_0g+[0;2*cos((iters-I)*pi/iters);0]];
i=1;
while x(1,i)<60
%     K=2/R*(-0.5+1/(1+exp(5*x(1,i)-100)))+0.02*(1-exp(-x(1,i)));
    K=2/R/pi*atan(100-.5*x(1,i)^2);
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
%     ind
    d=[d [s;ey;ephi]-x(:,i)];
    x(:,i)=[s;ey;ephi];
    U=[15*exp((10-x(1,i))/100);asin(L*K-0.1*x(2,i)-0.9*x(3,i))];
    
    x_next=state_update_bicycle(x(:,i),U);
    xb_next=state_update_bicycle(xb(:,i),U);
    x_nextg=state_update_bicycle_global(xg(:,i),U);

    x=[x x_next];
    xb=[xb xb_next];
    xg=[xg x_nextg];
    u_feasible=[u_feasible U];
    i=i+1;
end
x_data=[x_data x(:,1:end-1)];
d_data=[d_data d];
I=I+1;
end
%%
lip1=sdpvar(1);
gamm1=sdpvar(1);
lip2=sdpvar(1);
gamm2=sdpvar(1);
lip3=sdpvar(1);
gamm3=sdpvar(1);
cost1=1000*lip1^2+10*gamm1^2;
constr1=[lip1>=0;gamm1>=0];
cost2=1000*lip2^2+10*gamm2^2;
constr2=[lip2>=0;gamm2>=0];
cost3=1000*lip3^2+10*gamm3^2;
constr3=[lip3>=0;gamm3>=0];
for i=1:1:size(d_data,2)-1
    for j=i+1:1:size(d_data,2)
      constr1=[constr1; norm(d_data(1,i)-d_data(1,j))^2<=2*lip1*norm(x_data(1,i)-x_data(1,j))^2+8*gamm1];  
      constr2=[constr2; norm(d_data(2,i)-d_data(2,j))^2<=2*lip2*norm(x_data(1,i)-x_data(1,j))^2+8*gamm2];  
      constr3=[constr3; norm(d_data(3,i)-d_data(3,j))^2<=2*lip3*norm(x_data(1,i)-x_data(1,j))^2+8*gamm3];  
    
    end
end

ops=sdpsettings('verbose',0,'solver','ipopt');
sol1=optimize(constr1, cost1,ops)
sol2=optimize(constr2, cost2,ops)
sol3=optimize(constr3, cost3,ops)
