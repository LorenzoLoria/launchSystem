
function dsdt = launcherDynamics(t, x, thrustData, mission)
% LAUNCHERDYNAMICS  3D launcher equations of motion.
%   This function computes the time derivative of the state vector for a 
%   multistage rocket in a 3D Cartesian coordinate system.
%   Assumption: Flat Earth approximation.
%
%   State Vector x (7x1):
%     x(1) = pos_x   [m]   Position North/East relative
%     x(2) = pos_y   [m]   Position Cross-range
%     x(3) = pos_z   [m]   Vertical Position 
%     x(4) = vx      [m/s] Velocity x-component
%     x(5) = vy      [m/s] Velocity y-component
%     x(6) = vz      [m/s] Velocity z-component
%     x(7) = m       [kg]  Instantaneous total mass
%
%   thrustData:
%     .F                     = @(t) Thrust magnitude [N]
%     .thrustVectoringAngle  = @(t) Pitch angle [rad] (Angle from vertical z-axis or horizon depending on convention. Here: 0 = Vertical/Up)
%     .psi                   = @(t) Yaw angle   [rad] (Angle in x-y plane, 0 = along x-axis)
%
%   mission.capsule:
%     .Area   = Reference cross-sectional area [m^2]
%     .Cd     = Drag coefficient [-] 
%
%   mission.environment:
%     .altRange = Vector of altitudes [m]
%     .rhoVal   = Vector of air densities [kg/m^3]
%     .g0       = Standard gravity [m/s^2]
%
%   mission.prop:
%     .Isp      = Specific Impulse [s]

    % Unpack State Vector and Mission Data
    pos_x = x(1);
    pos_y = x(2);
    pos_z = x(3);
    vx    = x(4);
    vy    = x(5);
    vz    = x(6);
    m     = x(7);

    A   = mission.capsule.Area;
    Cd  = mission.capsule.Cd;
    g0  = mission.environment.g0;
    Isp = mission.prop.Isp;   % questo ancora non esiste

    % Interpolate air density based on current altitude
    h = max(pos_z, 0);  
    rho = interp1(mission.environment.altRange, mission.environment.rhoVal, h, 'linear', 'extrap');

    % Velocity magnitude
    vMag = sqrt(vx^2 + vy^2 + vz^2);

    % Get current thrust and steering angles from function handles
    T_mag = thrustData.F(t);                                      % [N] Thrust Magnitude
    thrustVectoringAngle = thrustData.thrustVectoringAngle(t);    % [rad] (0 = Vertical, pi/2 = Horizontal)
    psi   = thrustData.psi(t);                                    % [rad] Yaw angle (Azimuth)

    % Decompose Thrust Vector
    Tx = T_mag * sin(thrustVectoringAngle) * cos(psi);
    Ty = T_mag * sin(thrustVectoringAngle) * sin(psi);
    Tz = T_mag * cos(thrustVectoringAngle);

    % Drag contribution
    D_mag = 0.5 * rho * (vMag^2) * A * Cd; 
    
    Dx = -D_mag * (vx / vMag);
    Dy = -D_mag * (vy / vMag);
    Dz = -D_mag * (vz / vMag);

    % Gravity contribution
    Gx = 0;
    Gy = 0;
    Gz = -m * g0; 

    % mass flow rate
    mDot = -T_mag / (g0 * Isp); 

    % Equation of motion
    dsdt = zeros(7,1);
    
    % Velocity derivatives (Position rates)
    dsdt(1) = vx;   % dx/dt
    dsdt(2) = vy;   % dy/dt
    dsdt(3) = vz;   % dz/dt
    
    % Acceleration derivatives (Velocity rates)
    dsdt(4) = (Tx + Dx + Gx) / m;  
    dsdt(5) = (Ty + Dy + Gy) / m;  
    dsdt(6) = (Tz + Dz + Gz) / m; 
    
    % Mass derivative
    dsdt(7) = mDot; 

end