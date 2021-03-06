% simulate force on uncontrolled pendulum on cart using SS model of the
% linearized system equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean up matlab before launching script
clear all
close all 
clc

% pendulum point mass
m = 0.314;

% cart mass
M = 2;

% pendulum length
L = 0.64;

% acceleration due to gravity
g = -9.8;

% damping
d = 0.05;


% time steps of 100ms for integration
timeStep = 0.01;
% total time for simulation
totalTime = 1;
% build timepoint vector
tspan = 0:timeStep:totalTime;


titleMessage = 'SFC on inverted pendulum with Euler formula';
disp(titleMessage)


% initial conditions
% located at x=0
% velocity =0
% angle  pi (inverted)
% angular velocity = 0.5 rads-1
x0 = [0; 0.8;];
I = (m*L^2)/12;
den = I + (m*L^2);

A = [0 1;(-(m*g*L)/den) (-1*(d/den))];
B = [((m*L)/den); ((d*m*L)/den)];
C = [1; 0;];
D = 0;

% SFC gain calculation

PX = 8*[-1 -1.1];
K = place(A, B, PX);

% Luenberger gain

p = 20*[-1,-1.2];
Lu = place(A, C, p);


% use Euler's formula to solve the state space equation and generate values
% for a certain time with small time steps\difference in between values. 
[xhatout, yObserver] = SSMSolverEuler(A, B, C, D, K, Lu, tspan, x0);
% for all time point animate the results
range=1;
len = length(tspan);
kickFlag = zeros(1,len);

% get variables states
x = zeros(len, 1);    % cart positon
theta = xhatout(1, :);   % pendulum angle
th = transpose(theta);
% plot(th);
% return
% animate pendulum
figure
AnimatePendulumCart(th+pi,  x, L/2, tspan, range, kickFlag, titleMessage);
plot(tspan,th,'r-');
hold on;
plot(tspan, x, 'b-');
xlabel('time (in seconds)');
ylabel('Pendulum angle (Theta) and Cart velocity');
title('\fontsize{20}{\color{red}Change of angle theta against time}');
legend('{\color{red}Theta}', '{\color{blue}Cart velocity}');
grid on;
