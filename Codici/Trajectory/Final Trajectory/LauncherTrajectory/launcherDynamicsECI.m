
function dxdt = launcherDynamicsECI(t, x,thrustData, mission,stageNumber,opt)

% LAUNCHERDYNAMICS  3D launcher equations of motion.
%   This function computes the time derivative of the state vector for a 
%   multistage rocket in a 3D Cartesian coordinate system.
%   Assumption: Flat Earth approximation.
%
%
%   x axis coincident with radius connecting earth center and initial
%   position
%   
%   State Vector x (7x1):
%     x(1) = pos_x   [m]   x position
%     x(2) = pos_y   [m]   y position
%     x(3) = vx   [m]   velocity wrt to x axis
%     x(4) = vy      [m/s] velocity wrt to y axis
%     x(5) = theta      [m/s] attitude angle wrt x-axis
%     x(6) = omega      [m/s] angular velocity
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
    r     = x(1:2) ;
    v     = x(3:4) ;
    theta = x(5) ;
    omega = x(6) ;
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
    
   
throttling = optVar(1);
gammaGimbal =deg2rad(optVar(2));
nEngines = opt.stage{stageNumber}.nEngines ;
nominalThrust = opt.stage{stageNumber}.engine.thrust ;
theta = x(5) ;
BRFtoIRF = [cos(theta) -sin(theta) ; sin(theta) cos(theta)] ;
IRFtoBRF = BRFtoIRF' ;


    % alpha computation

    vBRF = IRFtoBRF * v ;
    alpha = atan2(vBRF(2), vBRF(1)) ;

    % M = norm(v) / soundSpeed
    % cd = f(M,alpha)
    % cl = f(M,alpha) 
    % xCP = f(M,alpha)
    % xCG = f(m)
    % inertiaCG = f(m)




thrustBRF = throttling * nEngines * nominalThrust * [cos(gammaGimbal) ; sin(gammaGimbal)] ; 

thrustIRF =  BRFtoIRF * thrustBRF;


    % Drag contribution
    dragIRF = - 0.5 .* rho .* norm(v) .* A .* Cd .* v; 

    % Lift contribution
    ClAlpha = 2*pi ; % Momentaneamente finchè non abbiamo funzione di Lucrezia
    alpha0 = 0 ; % Momentaneamente finchè non abbiamo funzione di Lucrezia

    liftIRF = 0.5 .* rho .* norm(v) .* A .* ClAlpha .* (alpha-alpha0) .* cross([v;0] , [0; 0; 1]);

    % liftIRF = 0.5 .* rho .* norm(v) .* A .* Cd .* v;
    
    % Gravity contribution
    gravityIRF = - GM * r /norm(r)^3;

    % mass flow rate
    mDot = - norm(ThrustIRF) / (g0 * opt.stage{stageNumber}.engine.isp); 

    % Equation of motion
    dxdt = zeros(7,1);
    
    % Velocity derivatives (Position rates)

    dxdt(1:2) = v;

    % Acceleration derivatives (Velocity rates)
    dxdt(3:4) = (ThrustIRF + dragIRF + liftIRF ) / m + gravityIRF;  

    
    dxdt(5) = omega ;



    % dxdt(6) = -norm(liftIRF) * (xCP - xCG) * cos(alpha) - norm(dragIRF) *
    % (xCP - xCG) * sin(alpha) - thrustBRF(2) * xCG ;
    % dxdt(6) = dxdt(6) / inertiaCG ;
    dxdt(6) = 0;



    % Mass derivative
    dxdt(7) = mDot; 

  
end