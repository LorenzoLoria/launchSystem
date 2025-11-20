function [mDot] = massIntegration(theta,m, vSpline, position, mission)
        % MASSINTEGRATION Computes the mass flow rate for trajectory integration.
    %
    %   [mDot] = MASSINTEGRATION(s, m, vSpline, position, mission) calculates the 
    %   instantaneous change in mass (dm/ds) required to follow the specified 
    %   trajectory kinematics. It solves the inverse dynamics problem to find 
    %   the required Thrust, and converts this to mass consumption based on Isp.
    %
    %   INPUTS:
    %       s        : scalar, Current independent variable (Theta/Angle in radians)
    %       m        : scalar, Current mass of the vehicle [kg]
    %       vSpline  : struct, Contains function handles for velocity profile
    %       position : struct, Contains function handles for pos, vel, acc vectors
    %       mission  : struct, Contains environmental and vehicle parameters
    %
    %   OUTPUTS:
    %       mDot     : scalar, Mass flow rate dm/ds [kg/rad]
    
    % Extrapolation of Data 
    

    vModule    = vSpline.profile ; 
    dvModule   = vSpline.dProfile ;
    pos        = position.position;
    dPosition  = position.dPosition; 
    ddPosition = position.ddPosition ; 

    surface = mission.capsule.Area;
    CD      = mission.capsule.supersonicCD;
    g0      = mission.environment.g0;
    iSp     = mission.launcher.engines{1}.isp;  
    GM      = mission.environment.GM;
    rEarth  = mission.environment.rEarth;

    thetaVec = vModule(:,1) ;
    finalIdx = thetaVec(end)==theta ;
    idx = sum((thetaVec-theta)<=0) - finalIdx ;
    thetaInitial = thetaVec(idx) ; 
    thetaFinal = thetaVec(idx+1) ; 
    
    pos = (pos(:,idx+1) - pos(:,idx)) / (thetaFinal - thetaInitial) * (theta-thetaInitial) + pos(:,idx) ;
    dPosition = (dPosition(:,idx+1) - dPosition(:,idx)) / (thetaFinal - thetaInitial) * (theta-thetaInitial) + dPosition(:,idx) ;
    ddPosition = (ddPosition(:,idx+1) - ddPosition(:,idx)) / (thetaFinal - thetaInitial) * (theta-thetaInitial) + ddPosition(:,idx) ;
    vModule = (vModule(idx+1) - vModule(idx)) / (thetaFinal - thetaInitial) * (theta-thetaInitial) + vModule(idx) ;
    dvModule = (dvModule(idx+1) - dvModule(idx)) / (thetaFinal - thetaInitial) * (theta-thetaInitial) + dvModule(idx) ;




    % Atmospheric Density (Exponential approximation for better physics)
    h = norm(pos) - rEarth; 
    rho = 1.225 * exp(-h/7200); 
    
    % % Gravity Vector
    % gVector = - GM * pos(s) /norm(pos(s))^3;
    % 
    % % Velocity Direction
    % velocityUnitVector = dPosition(s) ./ norm(dPosition(s)) ; 
    % 
    % % Acceleration Direction
    % accelerationVector = (dvModule(s) * vModule(s) / norm(dPosition(s))) * velocityUnitVector + (vModule(s) * vModule(s) / norm(dPosition(s))) * ( (ddPosition(s) * norm(dPosition(s)) - dPosition(s) * (dot(dPosition(s), ddPosition(s)) / norm(dPosition(s)))) / (norm(dPosition(s))^2) );
    % 
    % % Evaluation of Thrst
    % thrust = m * accelerationVector - m * gVector + 0.5 * rho * (vModule(s)).^2 * surface * CD .* velocityUnitVector ;


    % Gravity Vector
    gVector = - GM .* pos ./ norm(pos).^3 ; 


    % Velocity Direction
    velocityUnitVector = dPosition ./ norm(dPosition) ; 

    % Acceleration Direction
    accelerationVector = (dvModule * vModule / norm(dPosition)) * velocityUnitVector + (vModule * vModule / norm(dPosition)) * ( (ddPosition * norm(dPosition) - dPosition * (dot(dPosition, ddPosition) / norm(dPosition))) / (norm(dPosition)^2) );
    
    % Evaluation of Thrst
    thrust = m * accelerationVector - m * gVector + 0.5 * rho * (vModule).^2 * surface * CD .* velocityUnitVector ;

    % Mass Flow rate
    mDot = - 1 / iSp / g0 * (norm(thrust));



end


