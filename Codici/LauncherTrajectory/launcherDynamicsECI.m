
function dsdt = launcherDynamicsECI(t, x,thrustData, mission)
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
    r     = x(1:3);
    v     = x(4:6);
    m     = x(7);

    A   = mission.capsule.Area;
    Cd  = mission.capsule.Cd;
    g0  = mission.environment.g0;
    Isp = mission.launcher.engines{1}.isp;             % questo ancora non esiste
    GM  = mission.environment.GM;

    % Interpolate air density based on current altitude
    h   = norm(r)-mission.environment.rEarth;  
    rho = interp1(mission.environment.altRange, mission.environment.rhoVal, h, 'linear', 'extrap');

    T = thrustData(t); % thrustdata dovrebbe essere una funzione vettoriale con Tx,Ty,Tz

    % Drag contribution
    D = - 0.5 .* rho .* norm(v) .* A .* Cd .* v; 

    % Gravity contribution
    
    G = - GM * r /norm(r)^3;

    % mass flow rate
    mDot = - norm(T) / (g0 * Isp); 

    % Equation of motion
    dsdt = zeros(7,1);
    
    % Velocity derivatives (Position rates)

    dsdt(1:3) = v;

    % Acceleration derivatives (Velocity rates)
    dsdt(4:6) = (T + D + G) / m;  

    
    % Mass derivative
    dsdt(7) = mDot; 

end