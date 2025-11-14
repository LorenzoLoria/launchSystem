function []=launcherGeneration(x,mission)
% x(1) : engineType
% x(2) : tBurn 
% x(3:end) : throttlePercentage vector

engine = mission.launcher.engines{1,x(1)} ;
tBurn = x(2) ; 



m0 = 200 * 1e3 ;
g0 = 9.81;
if m0*g0 > engine.thrust
    nEngines = 1 ;
else
    nEngines = ceil(m0*g0/engine.thrust);
end



% Da spostare in un fmincon dedicato alla traiettoria
throttlePercentage = x(3:end);
if iscolumn(throttlePercentage)
    throttlePercentage = throttlePercentage';
end

throttleThroughTime = zeros(2,length(throttlePercentage)+1);
throttleThroughTime(1,:) = linspace(0, tBurn, length(throttleThroughTime(1,:)));
throttleThroughTime(2,:) = [1, throttlePercentage];

thrustProfile = @(t) nEngines .* engine.thrust .* interp1(throttleThroughTime(1,:), throttleThroughTime(2,:), t);

tStep = 1e-1 ; 
tSample = 0:tStep:tBurn ;

propellantMass = trapz(tSample,thrustProfile(tSample)/g0/engine.isp);

% Da spostare, diventa un vettore quando ci saranno più stadi
eps0 = 0.13 ;
structuralMass = eps0 * propellantMass / (1 - eps0) ;
totalMassInitialGuess = mission.capsule.weight + structuralMass + propellantMass ;

% serve guess iniziale mati e Tano, non guess iniziale massa che viene già
% trovata in aeroStructTrajAnalysis



[output] = fmincon(@(x) funzioneStructuralTrajAero)


thrustVec = output(1:length(output)/2) ;
thrustVectoringAngle = output(length(output)/2 + 1 : end) ;

[trajectory] = funzioneDiMatiETano(output) ;

[struttura] = funzioneDiFra(output) ;

[aerodinamica] = funzioneDiLucrezia(output) ;



% Da trajectory e strutture avremo tempi e masse per definire la nostra
% costFun per ga
tTotEMassTot = ;

end