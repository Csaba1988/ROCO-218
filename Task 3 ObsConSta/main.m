clear all
close all 
clc

% Define the variables for the calculations
% g is gravitational constant
g = 9.8;
% m is the mass of the pendulum rod
m = 0.314;
% l is the length of the rod
l = 0.64;
% mu is the viscous friction
mu = 0.05;
% I is the moment of inertia of a simple pendulum
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

A = [0 1; -a2 -a1;];
B = [b0; -(a1*b0);];
C = [1 0;];
D = 0;

%For OBSERVABILITY I need to create the observability matrix Mxo = [C;
%CA;..CA^(n-1);]n corresponds to the size of the matrix.
%If its rank is full the system is observable.Full rank means we cannot
%reduce any of the MXo matrix's line to 0s by algebraicly manipulating the
%rows.
disp('observability matrix for 2x2 A system matrix');
CA = C * A;
MXo = [C; CA;]
disp('rank of MXo');
rank(MXo)

%For CONTROLLABILITY I need to create the controllability matrix Mc = [B
%BA..BA^n]If its rank is full the system is observable.
disp('controllability matrix for 2x2 A system matrix');
AB = A * B;
MXc = [B AB]
disp('rank of MXc');
rank(MXc)

% For STABILITY I need to examine the eigenvalues of the the A matrix

% A_inverted is the A matrix when the angle theta(x1) is pi. The value in
% a2 is the value when the pendulum is in equilibrium at angle 0.So by
% making a2 negative we acquire the value that we would have if we would
% have linearized around pi.
A_inverted = [0 1; -a2 -a1;];
disp('A_inverted eigenvalues');
eig(A_inverted)

% The real part of the eigenvalues are -0.1794. They are marginally negetive hence the
% system is marginally stable. It will decay very slowly

% A_noninverted leaves a2 unchanged so this system matrix represents the
% matrix for when the pendulum is 'dangling' down theta(x1) is = to 0. And
% its eigenvalues show that it is not stable
A_noninverted = [0 1; a2 -a1;];
disp('A_noninverted eigenvalues');
eig(A_noninverted)

