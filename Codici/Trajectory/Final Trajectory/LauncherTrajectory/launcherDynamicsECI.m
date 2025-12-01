
function dsdt = launcherDynamicsECI(t, x,thrustData, mission,stageNumber,opt,option2D)

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
%
%     ThrustData(:,1) = Tx
%     ThrustData(:,2) = Ty
%     ThrustData(:,3) = Tz
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
    Cd  = mission.capsule.supersonicCD;
    g0  = mission.environment.g0; 
    GM  = mission.environment.GM;

    % Interpolate air density based on current altitude
    h   = norm(r)-mission.environment.rEarth;  
    %rho = interp1(mission.environment.altRange, mission.environment.rho, h, 'linear', 'extrap');
    rho = mission.environment.rhoFun(h);
    
    optVar = thrustData(t); % thrustdata dovrebbe essere una funzione vettoriale con Tx,Ty,Tz
    

if option2D == 1
    gammaGimball =deg2rad(optVar(2)) ;
    thetaGimball =0;
else
    gammaGimball =deg2rad(optVar(2));
    thetaGimball =deg2rad(optVar(3));
end

percVec = optVar(1);



ThrustBRF = percVec * opt.stage{stageNumber}.nEngines *opt.stage{stageNumber}.engine.thrust* [cos(thetaGimball)*cos(gammaGimball); cos(thetaGimball)*sin(gammaGimball); sin(thetaGimball)];


ThrustIRF = mission.Rfinal'*ThrustBRF;


    % Drag contribution
    D = - 0.5 .* rho .* norm(v) .* A .* Cd .* v; 

    % Gravity contribution
    G = - GM * r /norm(r)^3;

    % mass flow rate
    mDot = - norm(ThrustIRF) / (g0 * opt.stage{stageNumber}.engine.isp); 

    % Equation of motion
    dsdt = zeros(7,1);
    
    % Velocity derivatives (Position rates)

    dsdt(1:3) = v;

    % Acceleration derivatives (Velocity rates)
    dsdt(4:6) = (ThrustIRF + D ) / m + G;  

    
    % Mass derivative
    dsdt(7) = mDot; 

  
end