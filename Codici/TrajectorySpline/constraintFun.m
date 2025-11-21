function [C,Ceq] = constraintFun(x,thetaCoord,x0,xf, mission)
    
    % Estrapolation Data from Mission Struct
    g0  = mission.environment.g0;
    earthRadius = mission.environment.rEarth;

    % No Equality Constrains
    Ceq = [] ;

    % Flag for trajectoryGeneration
    fitnessFlag = 0 ; 
    
    [output] = trajectoryGenerationSpherical(x,thetaCoord, x0,xf,mission,fitnessFlag);
    
    tol = 10 ;
    earthRadius =  earthRadius + tol; 
    earthCenter = [0 ; 0 ; 0];
    
    % Output Data from trajectoryGeneration
    velocityUnitVector    = output.velocityUnitVector ;
    accelerationMagnitude = output.acceleration ; 
    thrust   = output.thrust ; 
    position = output.position ; 
    theta = output.theta ;  
    
    pos = zeros(3,length(theta)) ; 
    acc = zeros(length(theta),1) ;
    distance = zeros(length(theta),1) ;  
    steeringAngle   = zeros(length(theta),1) ; 
    thrustMagnitude = zeros(length(theta),1) ; 
    
    for ii = 1:length(theta)
        
        pos(:,ii) = position(:,ii);
        distance(ii,1) = sqrt(sum((pos(:,ii) - earthCenter).^2)) ;
        acc(ii,1) = accelerationMagnitude(ii)';
        steeringAngle(ii,1) = acosd(dot(thrust(:,ii) /norm(thrust(:,ii)), velocityUnitVector(:,ii))) ;
        thrustMagnitude(ii,1) = norm(thrust(:,ii));
        
    end
  
    maxAcc = 5 * g0 ; 
    maxSteeringAngle = 45 ; 
    maxThrust = mission.launcher.engines{1}.thrust * 4 ; 
    
    % Inequality Constrains
    % C = [(earthRadius - min(distance))/earthRadius ;...
    %     (max(acc) - maxAcc)/g0 ;...
    %     (max(thrustMagnitude) - maxThrust)/maxThrust;...
    %     (max(steeringAngle) - maxSteeringAngle)/maxSteeringAngle ];


    C = [(earthRadius - min(distance))/earthRadius ;...
        (max(acc) - maxAcc)/g0 ;...
        (max(thrustMagnitude) - maxThrust)/maxThrust];


    

end
