function x_next= state_update_bicycle(x,u)
R=20;
L=2;
dt=0.2;
K=2/R/pi*(atan(100-.5*x(1)^2));
x_next=x+dt*[u(1)*cos(x(3))/(1-x(2)*K);u(1)*sin(x(3));u(1)*(sin(u(2))/L-K*cos(x(3))/(1-x(2)*K))];
if x_next(3)>=pi
    x_next(3)=x_next(3)-2*pi;
   
elseif x_next(3)<=-pi
    x_next(3)=2*pi+x_next(3);
    
end
x_next(3)=max(-pi/2,min(pi/2,x_next(3)));