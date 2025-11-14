function dx_vect = dynamics_function(t, x_vect, theta, input, param, reference)

% INPUTS:
%   ~ = [1 x N] time vector (not used)
%   x_vect = [5 x N] state vectory
%       x_vect(1,:) = x = horizontal position
%       x_vect(2,:) = y = vertical position
%       x_vect(3,:) = dx = horizontal speed
%       x_vect(4,:) = dy = vertical speed
%       x_vect(5,:) = m = mass of the launcher
%   theta: steering angle
%   input: struct containing all the input data
%   param: struct containing basic parameters
%   reference: struct containing reference values for computation
% OUTPUTS:
%   dx_vect = d(x_vect)/dt = time-derivative of the state
%

% Time derivative of the state:
dx_vect = zeros(size(x_vect));
dx_vect(1:2,:) = x_vect(3:4,:);

% Compute the speed 
vx = x_vect(3,:);
vy = x_vect(4,:);
v = sqrt(vx.*vx + vy.*vy); 

% Mass 
m = input.M0-input.m_dot*t;

% Computation of the density
rho = reference.rho * exp(- x_vect(2, :) / reference.h);

% Drag
D = 0.5 * rho * v^2 * input.A * input.Cd;

% Computation of flight path angle
gamma = atan2(vy,vx);

% Compute the force vectors
fx = input.F / m * cos(theta(t)) - D / m * cos(gamma);
fy = input.F / m * sin(theta(t)) - D / m * sin(gamma) - param.g;  

% Record the accelerations (derivative of velocity states)
dx_vect(3,:) = fx;  
dx_vect(4,:) = fy;

% Mass flow rate
end