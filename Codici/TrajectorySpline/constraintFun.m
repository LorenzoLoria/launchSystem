function [C,Ceq] = constraintFun(x,xCoord,x0,x_target, mission)
    
    % Estrapolation Data from Mission Struct
    g0  = mission.environment.g0;
    earthRadius = mission.environment.rEarth;

    % No Equality Constrains
    Ceq = [] ;

    % Flag for trajectoryGeneration
    fitnessFlag = 0 ; 
    
    [output] = trajectoryGeneration(x,xCoord, x0,x_target,mission,fitnessFlag);
    
    tol = 0.25;
    earthRadius =  earthRadius + tol; 
    earthCenter = [0 ; 0 ; 0];
    
    % Output Data from trajectoryGeneration
    velocityUnitVector    = output.velocityUnitVector ;
    accelerationMagnitude = output.acceleration ; 
    thrust   = output.thrust ; 
    position = output.position ; 
    s = output.s ;  
    
    pos = zeros(3,length(s)) ; 
    acc = zeros(length(s),1) ;
    distance = zeros(length(s),1) ;  
    steeringAngle   = zeros(length(s),1) ; 
    thrustMagnitude = zeros(length(s),1) ; 
    
    for ii = 1:length(s)
        
        pos(:,ii) = position(s(ii));
        distance(ii,1) = sqrt(sum((pos(:,ii) - earthCenter).^2)) ;
        acc(ii,1) = accelerationMagnitude(ii)';
        steeringAngle(ii,1) = acosd(dot(thrust(:,ii) /norm(thrust(:,ii)), velocityUnitVector(:,ii))) ;
        thrustMagnitude(ii,1) = norm(thrust(:,ii));
        
    end
  
    maxAcc = 5 * g0 ; 
    maxSteeringAngle = 45 ; 
    maxThrust = 2500 ; 
    
    % Inequality Constrains
    C = [(earthRadius - min(distance))/earthRadius ;...
        (max(acc) - maxAcc)/g0 ;...
        (max(thrustMagnitude) - maxThrust)/maxThrust;...
        (max(steeringAngle) - maxSteeringAngle)/maxSteeringAngle ];
    

end
