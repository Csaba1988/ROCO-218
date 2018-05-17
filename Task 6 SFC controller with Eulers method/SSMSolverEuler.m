function [xout, y] = SSMSolverEuler (A, B, C, D, K,tspan, x0)


%get signal length
len = length(tspan);

%initialize the matrices for position and angle
y = zeros(1, len);
xout = zeros(2, len);

%put the first value from xout into x0 (initial state)
xout(:, 1)=x0;
x = x0

% The input of the system is -KX. This line is using the initial conditions
% to calculate the first value
u(1)= -K(1)*x(1) - K(2)*x(2);

%y is the output of the system which is the angle of the pendulum. The
%first value is the starting position set by the initial codition.
y(1)= C(1)*x(1) + C(2)*x(2) + D(1) * u(1);

% This loop is calculating the rest of the values for the angle
for idx = 2:len
    % Calculate the input for each iteration based on current state
    u(idx) = -K(1)*x(1) - K(2)*x(2);
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
    
    % Record the states
    xout(:,idx) = x;
    % Record current output theta.
    y(idx)= C(1)*x(1) + C(2)*x(2) + D(1) * u(idx);
end

end