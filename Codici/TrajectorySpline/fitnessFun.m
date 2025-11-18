function [objective] = fitnessFun(x,xCoord,x0,x_target, mission)
   
    % flag for trajectoryGeneration
    fitnessFlag = 1 ;
    
    [output] = trajectoryGeneration(x,xCoord,x0,x_target,mission,fitnessFlag);
    
    vel = output.vel ; 
    yCoord = output.yCoord ; 
    zCoord = output.zCoord ; 
    s = output.s ;
    
    % Evalutation of the total final time
    tFinal = trapz(s,abs(1./vel));
   
    % Evalutation of oscillation on y and z axis
    yOscillation = sqrt(sum(diff(yCoord).^2));
    zOscillation = sqrt(sum(diff(zCoord).^2));
    
    %objective = 1/m(end) + tFinal + yOscillation / 15 + zOscillation / 8 ; 
    objective = tFinal + yOscillation / 5 + zOscillation / 3 ; 


end
