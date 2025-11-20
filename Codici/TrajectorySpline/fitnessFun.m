function [objective] = fitnessFun(x,thetaCoord,x0,xf, mission)
   
    % flag for trajectoryGeneration or trajectoryGenerationSpherical
    fitnessFlag = 1 ;
    
    %[output] = trajectoryGeneration(x,xCoord,x0,x_target,mission,fitnessFlag);
    [output] = trajectoryGenerationSpherical(x, thetaCoord, x0, xf, mission, fitnessFlag);

    % Data Estrapolation
     rCoord   = output.rCoord;
     phiCoord = output.phiCoord;
     theta    = output.theta; 
     vel      = output.vel ; 
    
    % Evalutation of the total final time
    tFinal = trapz(theta,abs(1./vel));
   
    % Evalutation of oscillation on y and z axis
    rOscillation = sqrt(sum(diff(rCoord).^2));
    phiOscillation = sqrt(sum(diff(phiCoord).^2));
    vOscillation = sqrt(sum(diff(vel).^2));

    
    %objective = 1/m(end) + tFinal + yOscillation / 15 + zOscillation / 8 ; 
    objective = tFinal + rOscillation / 5 + phiOscillation / 3 + vOscillation / 10 ; 


end
