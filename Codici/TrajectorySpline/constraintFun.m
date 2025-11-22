function [C,Ceq] = constraintFun(x,x0,target, mission)
    
    nOptPoints = (length(x)-1)/3 ;
    thetaCoord = linspace(x0(2), x(end), nOptPoints+1) ;

    % Estrapolation Data from Mission Struct
    g0  = mission.environment.g0;
    earthRadius = mission.environment.rEarth;

    % No Equality Constrains
    Ceq = [] ;

    % Flag for trajectoryGeneration
    fitnessFlag = 0 ; 
    
    [output] = trajectoryGenerationSpherical(x,thetaCoord, x0,target,mission,fitnessFlag);
    
    if output.unpheasibleFlag
        C = output.value ;
        return
    end

    tol = 10 ;
    earthRadius =  earthRadius + tol; 
    earthCenter = [0 ; 0 ; 0];
    
    % Output Data from trajectoryGeneration
    velocityUnitVector    = output.velocityUnitVector ;
    accelerationMagnitude = output.acceleration ; 
    thrust   = output.thrust ; 
    position = output.position ; 
    sVec = output.sVec ;  
    
    pos = zeros(3,length(sVec)) ; 
    acc = zeros(length(sVec),1) ;
    distance = zeros(length(sVec),1) ;  
    steeringAngle   = zeros(length(sVec),1) ; 
    thrustMagnitude = zeros(length(sVec),1) ; 
    
    for ii = 1:length(sVec)
        
        pos(:,ii) = position(:,ii);
        distance(ii,1) = sqrt(sum((pos(:,ii) - earthCenter).^2)) ;
        acc(ii,1) = accelerationMagnitude(ii)';
        steeringAngle(ii,1) = acosd(dot(thrust(:,ii) /norm(thrust(:,ii)), velocityUnitVector(:,ii))) ;
        thrustMagnitude(ii,1) = norm(thrust(:,ii));
        
    end
  
    maxAcc = 5 * g0 ; 
    maxSteeringAngle = 45 ; 
    maxThrust = mission.launcher.engines{1}.thrust * 4 ; 
    minVel = 0 ;
    
    % Inequality Constrains
    % C = [(earthRadius - min(distance))/earthRadius ;...
    %     (max(acc) - maxAcc)/g0 ;...
    %     (max(thrustMagnitude) - maxThrust)/maxThrust;...
    %     (max(steeringAngle) - maxSteeringAngle)/maxSteeringAngle ];


    C = [(earthRadius - min(distance))/earthRadius ;...
        (max(acc) - maxAcc)/g0 ;...
        (max(thrustMagnitude) - maxThrust)/maxThrust ;...
        - min(output.vel) + minVel];


    

end
