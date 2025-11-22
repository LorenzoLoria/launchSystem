function [objective] = fitnessFun(x,x0,target, mission)
   

    nOptPoints = (length(x)-1)/3 ;
    thetaCoord = linspace(x0(2), x(end), nOptPoints+1) ;


    % flag for trajectoryGeneration or trajectoryGenerationSpherical
    fitnessFlag = 1 ;
    
    %[output] = trajectoryGeneration(x,xCoord,x0,x_target,mission,fitnessFlag);
    [output] = trajectoryGenerationSpherical(x, thetaCoord, x0, target, mission, fitnessFlag);

    if output.unpheasibleFlag
        objective = output.value ;
        return
    end



    % Data Estrapolation
     rCoord             = output.rCoord;
     phiCoord           = output.phiCoord;
     vel                = output.vel ; 
     sVec               = output.sVec ;

   % meanVel = (vel(1:end-1) + vel(2:end))./2 ;
    % Evalutation of the total final time
    tFinal = trapz(sVec,1./(vel));
   
    % Evalutation of oscillation on y and z axis
    rOscillation = sqrt(sum(diff(rCoord).^2));
    phiOscillation = sqrt(sum(diff(phiCoord).^2));
    vOscillation = sqrt(sum(diff(vel).^2));



    rEarth = mission.environment.rEarth ;

    weightPhiOscillations = max(phiCoord) - 0.999 * min(phiCoord) ; 
    weightVOscillations = max(vel) - 0.999 * min(vel) ; 

    %objective = 1/m(end) + tFinal + yOscillation / 15 + zOscillation / 8 ; 
    objective = tFinal + rOscillation / rEarth + phiOscillation / weightPhiOscillations + vOscillation / weightVOscillations ; 


end
