function  [objective] = objFunGlobalGA(launcher,mission,settings)

nlconFlag = 0;

[output] = launcherSimulation(launcher,mission,settings,nlconFlag);

if ~isstruct(output)
    keyboard
end



totalMass = output.launcherMass ; 
tof = output.tof ; 


objective = totalMass/200000 + tof/3600 ; 

if isnan(totalMass) || isnan(tof)
    keyboard
end

if isinf(totalMass) || isinf(tof)
    keyboard
end

if isempty(totalMass) || isempty(tof)
    keyboard
end







end
