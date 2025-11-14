function dxVect = dynFunctRocket(t, stateVect, thrustAngleFun, rocketInput, envParam, atmRef)
% dynFunctRocket
% Computes the time derivative of the rocket state.
%
% State vector structure:
%   stateVect(1,:) = x-position [m]
%   stateVect(2,:) = y-position [m]
%   stateVect(3,:) = x-velocity [m/s]
%   stateVect(4,:) = y-velocity [m/s]
%
% Inputs:
%   t               - time [s]
%   stateVect       - state vector (4 x N)
%   thrustAngleFun  - function handle: theta(t) [rad]
%   rocketInput     - struct with fields:
%                       .M0   : initial mass [kg]
%                       .mDot : mass flow rate [kg/s]
%                       .F    : thrust magnitude [N]
%                       .A    : reference area [m^2]
%                       .Cd   : drag coefficient [-]
%   envParam        - struct with fields:
%                       .g    : gravitational acceleration [m/s^2]
%   atmRef          - struct with fields:
%                       .rho  : air density at reference level [kg/m^3]
%                       .h    : scale height [m]

    % Preallocate derivative of the state
    dxVect = zeros(size(stateVect));

    % Unpack state for readability
    posX = stateVect(1, :);   % x-position
    posY = stateVect(2, :);   % y-position
    velX = stateVect(3, :);   % x-velocity
    velY = stateVect(4, :);   % y-velocity

    % Position derivative = velocity
    dxVect(1, :) = velX;
    dxVect(2, :) = velY;

    % Compute speed magnitude
    speed = sqrt(velX.^2 + velY.^2);

    % Current mass (simple linear decrease with time)
    mass = rocketInput.M0 - rocketInput.mDot * t;

    % Atmospheric density (exponential model)
    airDensity = atmRef.rho .* exp(- posY ./ atmRef.h);

    % Aerodynamic drag magnitude (opposes velocity)
    dragForce = 0.5 .* airDensity .* speed.^2 .* rocketInput.A .* rocketInput.Cd;

    % Flight path angle (direction of velocity vector)
    flightPathAngle = atan2(velY, velX);

    % Thrust direction from user-supplied function (inertial frame angle)
    thrustAngle = thrustAngleFun(t);

    % Accelerations in x and y
    accelX = (rocketInput.F ./ mass) .* cos(thrustAngle) ...
           - (dragForce ./ mass) .* cos(flightPathAngle);

    accelY = (rocketInput.F ./ mass) .* sin(thrustAngle) ...
           - (dragForce ./ mass) .* sin(flightPathAngle) ...
           - envParam.g;

    % Record acceleration components (derivative of velocity)
    dxVect(3, :) = accelX;
    dxVect(4, :) = accelY;

end
