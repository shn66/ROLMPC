function [P, C,U, X, Q, R, N, Ny] = DefineSystem_bicycle()
% Set controller horizon
N = 8;

% Set cost function
Q = 100;
  
% R = 2*[1 0;0 1];
R = 10;
P=10;

C=[1 0 0;0 1 0];
% P=10*[1 0;0 1];
Ny=2;
% Fu= @(Y) [1 -2 1];
% Define the input constraint set U
u1_max = 18.0;
u1_min = 0.0;
u2_max=pi/2;

U = Polyhedron([u1_max u2_max;u1_max -u2_max;u1_min u2_max;u1_min -u2_max]);

% Define the Recovery input constraint set R
% Define the state constraint set X
X=Polyhedron('A',[1 0 0;-1 0 0;0 1 0;0 -1 0; 0 0 1; 0 0 -1], 'b', [60;1.0;5;5;pi/2-pi/60;pi/2-pi/60]);
% x_max = 15;
% v_max = 15;
% X = Polyhedron([x_max v_max; x_max -v_max; -x_max -v_max; -x_max v_max]);
end