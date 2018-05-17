function [xhatout, yObserver] = SSMSolverEuler (A, B, C, D, K, Lu, tspan, x0)
% Solves State Space model by Euler integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%% System %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get signal length
len = length(tspan);

%initialize the matrices for position and angle
y = zeros(1, len);

%put the first value (initial value for state) into variable x
x = x0;

% Initialize the output and state variable for the observer
yObserver = zeros(1, len);
xhatout = zeros(3, len);

%put the first value from xout into x0 (initial state)
xhatout(:, 1)=x0;
xhat = x0;

% Input U is -K xhat the estimated state variables
u(1)= -K(1)*xhat(1) - K(2)*xhat(2)- K(3)*xhat(3);

%y is the output of the real system which is the angle of the pendulum. The
%first value is the starting position set by the initial codition.
y(1)= C(1)*x(1) + C(2)*x(2) + C(3)*x(3)+ D(1) * u(1);

% This loop is calculating the rest of the values for the angle 
for idx = 2:len
    %%%%%%%%%%%%%%%%%%%%%%%%%% Real System %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate the input for each iteration based on current state
    u(idx) = -K(1)*xhat(1) - K(2)*xhat(2)- K(3)*xhat(3);
    
    % Calculate step size
    h = tspan(idx) - tspan(idx-1);
    
    % It is effectively = to Xdot = AX + BU but written in C language
    % style. I have included the hole calculations but the multiplication
    % with the A matrix in xdot3 can be omitted as it is = 0.
    xdot(1) = A(1,1)*x(1) + A(1,2)*x(2) + A(1,3)*x(3) + B(1)*u(idx);
    xdot(2) = A(2,1)*x(1) + A(2,2)*x(2) + A(2,3)*x(3) + B(2)*u(idx);
    xdot(3) = A(3,1)*x(1) + A(3,2)*x(2) + A(3,3)*x(3) + B(3)*u(idx);
    
    % Calculate the next state variable based on the current one using
    % Euler integration
    x(1)= x(1)+h*xdot(1);
    x(2)= x(2)+h*xdot(2);
    x(3)= x(3)+h*xdot(3);
    
    % Record current output theta.
    y(idx)= C(1)*x(1) + C(2)*x(2) + C(3)*x(3) + D(1) * u(idx);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Leunberger observer %%%%%%%%%%%%%%%%%%%%
    
    % Calculate the error between the estimated output and the real systems
    % output
    errorCorrection = y(idx)- C(1)*xhat(1) - C(2)*xhat(2)- C(3)*xhat(3);
    
    % Calculate the estimated values of the real system xhatdot. And correct the
    % estimated outuput by the error * a gain(Luenberger gain). As position
    % is only simulated and y only includes the angle I did not have a
    % Luenberger gain for xhatdot3. As the 3rd row of A matrix is 0 the
    % multiplication with A matrix can be omitted from xhatdot3
    xhatdot(1) = A(1,1)*xhat(1) + A(1,2)*xhat(2) + A(1,3)*xhat(3) + B(1)*u(idx) + Lu(1)*errorCorrection; 
    xhatdot(2) = A(2,1)*xhat(1) + A(2,2)*xhat(2) + A(2,3)*xhat(3) + B(2)*u(idx) + Lu(2)*errorCorrection;
    xhatdot(3) = A(3,1)*xhat(1) + A(3,2)*xhat(2) + A(3,3)*xhat(3) + B(3)*u(idx)
    
    % Calculate the next state based on current state using Euler
    % integration
    xhat(1)= xhat(1)+h*xhatdot(1);
    xhat(2)= xhat(2)+h*xhatdot(2);
    xhat(3)= xhat(3)+h*xhatdot(3);
    
    % Record the states in xhatout
    xhatout(:,idx) = xhat;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


end