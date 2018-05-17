function xDot = SSSimulate1 (X, A, B, K)
%Function describes the state space model

% x is state
% A,B state matrices
% u control input
% xdot is the returned 1st derivative of state

xDot = A*X + B*(-K*X) + [0; 0.01 *rand;];
end