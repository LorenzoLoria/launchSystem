function [objective] = fitnessFun(x,thetaCoord,x0,xf, mission)
   
    % flag for trajectoryGeneration or trajectoryGenerationSpherical
    fitnessFlag = 1 ;
    
    %[output] = trajectoryGeneration(x,xCoord,x0,x_target,mission,fitnessFlag);
    [output] = trajectoryGenerationSpherical(x, thetaCoord, x0, xf, mission, fitnessFlag);

    % Data Estrapolation
     rCoord             = output.rCoord;
     phiCoord           = output.phiCoord;
     theta              = output.theta; 
     vel                = output.vel ; 
     pos                = output.position' ; 
     trajectoryLength   = output.trajectoryLength ;
    

   % meanVel = (vel(1:end-1) + vel(2:end))./2 ;
    % Evalutation of the total final time
    tFinal = trapz(trajectoryLength,1./(vel));
   
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
