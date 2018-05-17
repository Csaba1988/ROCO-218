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
timeStep = 0.1;
% total time for simulation
totalTime = 10;
% build timepoint vector
tspan = 0:timeStep:totalTime;


titleMessage = 'uncontrolled linear sim of FUC pendulum on cart';
disp(titleMessage)


% initial conditions
% located at x=0
% velocity =0
% angle  pi (inverted)
% angular velocity = 0.5 rads-1
y0 = [pi; 0;];
I = (m*L^2)/12;
den = I + (m*L^2);

A = [0 1;(-(m*g*L)/den) (-1*(d/den))];
B = [((m*L)/den); ((d*m*L)/den)];

% SFC gain calculation

PX = 50*[-1 -1.1];
K = place(A, B, PX);


% use ode to solve with FCPendOnCart with no control force input u
% representing a force controlled pendulum on a cart
% model introduces slight amount of noise to wont stay balanced
[t,y] = ode45(@(t,y)SSSimulate(y, A, B, K),tspan,y0);

% for all time point animate the results
range=1;
len = length(tspan);
kickFlag = zeros(1,len);

% get variables states
x = 0 * y(:, 2);    % cart positon
th = y(:, 1);   % pendulum angle
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
