function [objective] = fitnessFun(x,sCoord,x0,x_target, mission)
   
    % flag for trajectoryGeneration or trajectoryGenerationSpherical
    fitnessFlag = 1 ;
    
    %[output] = trajectoryGeneration(x,xCoord,x0,x_target,mission,fitnessFlag);
    [output] = trajectoryGenerationSpherical(x, sCoord, x_target, x0, x_target, mission, fitnessFlag);

    % Data Estrapolation
      rCoord   = output.rCoord;
      phiCoord = output.phiCoord;
      s        = output.s; 
      vel      = output.vel ; 
    
    % Evalutation of the total final time
    tFinal = trapz(s,abs(1./vel));
   
    % Evalutation of oscillation on y and z axis
    yOscillation = sqrt(sum(diff(yCoord).^2));
    zOscillation = sqrt(sum(diff(zCoord).^2));
    
    %objective = 1/m(end) + tFinal + yOscillation / 15 + zOscillation / 8 ; 
    objective = tFinal + yOscillation / 5 + zOscillation / 3 ; 


end
