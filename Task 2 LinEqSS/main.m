% Define the variables for the calculations
% g is gravitational constant
g = 9.8;
% m is the mass of the pendulum rod
m = 0.314;
% l is the length of the rod
l = 0.64;
% mu is the viscous friction
mu = 0.05;
% I is the moment of inertia of a rod pendulum
I = (m * l^2)/12;
% b0 is the the constant for the control input variable
b0 = m*l/(I + m * l^2);
% b1 = 0;
% a0 = 0;
% a1 is the constant of x2 state
a1 = mu/(I + m * l^2);
% a2 is the constant for x1 state
a2 = m*l*g/(I + m * l^2);

% define the linearize state space matrices from differential equation

A = [0 1; -a2 -a1;]
B = [b0; -(a1*b0);]
C = [1 0;]
D = 0
