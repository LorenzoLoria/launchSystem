function [mDot] = massIntegration(s,m,vSpline, position, mission)
    
    % Extrapolation of Data 
    vModule    = vSpline.polynomial ; 
    dvModule   = vSpline.firstDer ;
    pos        = position.position;
    dPosition  = position.dPosition; 
    ddPosition = position.ddPosition ; 

    surface = mission.capsule.Area;
    CD      = mission.capsule.supersonicCD;
    g0      = mission.environment.g0;
    iSp     = mission.launcher.engines{1}.isp;  

    % Interpolate air density based on current altitude
    h   = norm(pos(s))-mission.environment.rEarth;  
    rho = mission.environment.gridInterp(h);
    
    % Gravity contribution
    gVector = - GM * pos(s) /norm(pos(s))^3;
    
    velocityUnitVector = dPosition(s) ./ norm(dPosition(s)) ; 
    accelerationVector = (dvModule(s) * vModule(s) / norm(dPosition(s))) * velocityUnitVector + (vModule(s) * vModule(s) / norm(dPosition(s))) * ( (ddPosition(s) * norm(dPosition(s)) - dPosition(s) * (dot(dPosition(s), ddPosition(s)) / norm(dPosition(s)))) / (norm(dPosition(s))^2) );
    
    % Evaluation of Thrst
    thrust = m * accelerationVector - m * gVector + 0.5 * rho * (vModule(s)).^2 * surface * CD .* velocityUnitVector ;
    
    % Mass Flow rate
    mDot = - 1 / iSp / g0 * (norm(thrust));

end


