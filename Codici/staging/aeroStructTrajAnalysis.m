function [] = aeroStructTrajAnalysis(x, engine, tBurn )

throttlePercentage = x(1:length(x)/2) ;
thrustVectoringAngle = x(length(x)/2+1 : end) ;

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


[] = funzioneDiMatiETano()


[] = funzioneDiFra()


[] = funzioneDiLucrezia()
















end