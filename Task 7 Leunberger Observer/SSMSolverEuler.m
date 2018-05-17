function [xhatout, yObserver] = SSMSolverEuler (A, B, C, D, K, Lu, tspan, x0)
% Solves State Space model by Euler integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%% System %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get signal length
len = length(tspan);

%initialize the output matrix for the real system
y = zeros(1, len);

%put the first value (initial value for state) into variable x
x = x0;

yObserver = zeros(1, len);
xhatout = zeros(2, len);

%put the first value from xout into x0 (initial state)
xhatout(:, 1)=x0;
xhat = x0;

% Input U is -K xhat the estimated state variables
u(1)= -K(1)*xhat(1) - K(2)*xhat(2);

%y is the output of the real system which is the angle of the pendulum. The
%first value is the starting position set by the initial codition.
y(1)= C(1)*x(1) + C(2)*x(2) + D(1) * u(1);

% This loop is calculating the rest of the values for the angle 
for idx = 2:len
    %%%%%%%%%%%%%%%%%%%%%%%%%% Real System %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate the input for each iteration based on current state
    u(idx) = -K(1)*xhat(1) - K(2)*xhat(2);
    
    % Calculate step size
    h = tspan(idx) - tspan(idx-1);
    
    % It is effectively = to Xdot = AX + BU but written in C language
    % style
    xdot(1) = A(1,1)*x(1) + A(1,2)*x(2) + B(1)*u(idx);
    xdot(2) = A(2,1)*x(1) + A(2,2)*x(2) + B(2)*u(idx);
    
    % Calculate the next state variable based on the current one using
    % Euler integration
    x(1)= x(1)+h*xdot(1);
    x(2)= x(2)+h*xdot(2);

    % Record current output theta.
    y(idx)= C(1)*x(1) + C(2)*x(2) + D(1) * u(idx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Leunberger observer %%%%%%%%%%%%%%%%%%%%%
    
    % Calculate the error between the estimated output and the real systems
    % output
    errorCorrection = y(idx)- C(1)*xhat(1) - C(2)*xhat(2);
    
    % Calculate the estimated values of xdot => xhatdot. And correct the
    % estimated outuput by the error * a gain(Luenberger gain)
    xhatdot(1) = A(1,1)*xhat(1) + A(1,2)*xhat(2) + B(1)*u(idx) + Lu(1)*errorCorrection; 
    xhatdot(2) = A(2,1)*xhat(1) + A(2,2)*xhat(2) + B(2)*u(idx) + Lu(2)*errorCorrection;
    
    % Calculate the next state based on current state using Euler
    % integration
    xhat(1)= xhat(1)+h*xhatdot(1);
    xhat(2)= xhat(2)+h*xhatdot(2);
    
    % Record the states in xhatout
    xhatout(:,idx) = xhat;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
end


end