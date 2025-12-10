function dxdt = launcherDynamicsAndControlECI(t, x,thrustData, mission, configuration, mer, launcher, stageNumber, guidancePoints, gains)
% LAUNCHERDYNAMICSOL  3D launcher equations of motion.
%   This function computes the time derivative of the state vector for a 
%   multistage rocket in a 3D Cartesian coordinate system.
%
%   x axis coincident with radius connecting earth center and initial
%   position
%   
%   State Vector x (7x1):
%     x(1) = pos_x   [m]   x position
%     x(2) = pos_y   [m]   y position
%     x(3) = vx      [m]   velocity wrt to x axis
%     x(4) = vy      [m/s] velocity wrt to y axis
%     x(5) = theta   [m/s] attitude angle wrt x-axis
%     x(6) = omega   [m/s] angular velocity
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
%   mission.environment:
%     .altRange = Vector of altitudes [m]
%     .rhoVal   = Vector of air densities [kg/m^3]
%     .g0       = Standard gravity [m/s^2]
%
%   guidancePoint: desire position from the optimal trajectory

    % Unpack State Vector and Mission Data
      r     = x(1:2) ;
      v     = x(3:4) ;
      theta = x(5) ;
      omega = x(6) ;
      m     = x(7);
       
     rMag = sqrt(r'*r);
     vMag = sqrt(v'*v); 
        
      A   = mission.capsule.Area;
      g0  = mission.environment.g0; 
      GM  = mission.environment.GM;
        
    totalStageNumber = launcher(1) ;

    inertia =  InertaEvaluation(mission, configuration, mer, launcher, totalStageNumber );

    vTarget = guidancePoints(3:4, end) ;

    % Interpolate air density based on current altitude
    h   = rMag - mission.environment.rEarth;  
    rho = mission.environment.rhoFun(h);
    
    dynamicPressure = 0.5 * rho * vMag^2;
    [soundspeed] = mission.aerodynamics.soundspeedFun(h);
    Mach = vMag/soundspeed;

    BRFtoIRF = [cos(theta) -sin(theta) ; sin(theta) cos(theta)] ;
    IRFtoBRF = BRFtoIRF' ;

    % alpha computation
    vBRF = IRFtoBRF * v ;
    alpha = atan2(vBRF(2), vBRF(1)) ;  
  
    
    % Evaluation of the Aerodynamics coefficients
    if Mach == 0
        Cd = 0.01;
    else
       [Cl,Cd,~,~] = CLCDcomputation(Mach,alpha,dynamicPressure,1,mission,stageNumber,configuration);
       % ClAlpha = 2*pi ; 
    end

    optVar = thrustData(t); 
   
    throttling  = optVar(1);
    gammaGimbal = deg2rad(optVar(2));
    nEngines    = configuration.stage{stageNumber}.nEngines ;
   
    if stageNumber == 1 % Sea-Level Parameters
      
        staticContribution = (101325-mission.environment.pressure(h))*configuration.stage{stageNumber}.engine.effAreaZero;
        nominalThrust = configuration.stage{stageNumber}.engine.thrustZero ;
        iSp = configuration.stage{stageNumber}.engine.ispZero ;
    else % Vaccum Parameters

       staticContribution = 0;
       nominalThrust = configuration.stage{stageNumber}.engine.thrustVacuum ;
       iSp = configuration.stage{stageNumber}.engine.ispVacuum ;
    end


    % xCG = f(m) (xCG = 30)
    xCG = computeXCG(mission, configuration) ;
    
    % xCP = f(M,alpha) (xCP = 10)
    xCP = computeXcp(mission, configuration) ;
 
    % inertiaCG = f(m) - I_xx = 1/12 * m * (3*r^2 + h^2) (inertiaCG = 27 * 1e6)
    % inertiaCG = (1/12) * m * (3*launcherRadius^2 + launcherHeight^2);

    thrustBRF = throttling * nEngines * (nominalThrust + staticContribution)  * [cos(gammaGimbal) ; sin(gammaGimbal)] ;
    thrustIRF =  BRFtoIRF * thrustBRF;

    % Drag contribution
    dragIRF = - 0.5 .* rho .* vMag .* A .* Cd .* v; 

    % Lift contribution Reference
    alpha0 = 0 ;

    liftIRF = 0.5 .* rho .* vMag .* A .* (Cl .* (alpha-alpha0).* (rad2deg(alpha)<15) ).* cross([v;0] , [0; 0; 1]);
    liftIRF = liftIRF(1:2) ;
    
    % Gravity contribution
    gravityIRF = - GM * r / rMag^3 ;

    % mass flow rate
    mDot = - norm(thrustIRF) / (g0 *  iSp); 
    
    % Guidance e Control: Simple PD guidance to get a desired acceleration vector
    aCommandECI = gains.Kp_pos * (guidancePoints(1:2) - r) + ...
                  gains.Kd_vel * (guidancePoints(3:4) - v) - gravityIRF;
    
    % desired theta: the LV shall be point in order to release the capsule
    % in the proper direction
    thetaDesired = atan2(vTarget(2), vTarget(1));
     
    torqueCommand = gains.Kp_theta * (thetaDesired - theta) + gains.Kp_omega* (0 - omega) ;

    % For Guidance 
    thrustGuidanceBR = IRFtoBRF * (m * aCommandECI) ;
    
    % For Attitude
    thrustCommandBR = - torqueCommand / xCG ; 

    thrustBRF = [thrustGuidanceBR; thrustCommandBR];

    % Saturate to available thrust
    if norm(thrustBRF) > nominalThrust  && norm(thrustBRF) > 0
        thrustBRF = thrustBRF / norm(thrustBRF) * nominalThrust ;
    end
    
    % Need thrust in the Inertial Frames
    thrustIRF = BRFtoIRF * thrustBRF ;

    % Equation of motion
    dxdt = zeros(7,1);
    
    % Velocity derivatives
    dxdt(1:2) = v;

    % Acceleration derivatives
    dxdt(3:4) = (thrustIRF + dragIRF + liftIRF ) / m + gravityIRF;  
    
    dxdt(5) = omega ;

    dxdt(6) = -norm(liftIRF) * (xCP - xCG) * cos(alpha) - norm(dragIRF) *...
    (xCP - xCG) * sin(alpha) - thrustBRF(2) * xCG ;
    
    dxdt(6) = dxdt(6) / inertia ;
    
    % Mass Derivative
    dxdt(7) = mDot ;

end