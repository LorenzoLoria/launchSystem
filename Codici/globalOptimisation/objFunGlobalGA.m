function  [objective] = objFunGlobalGA(launcher,mission,settings)

nlconFlag = 0;

[output] = launcherSimulation(launcher,mission,settings,nlconFlag);

if ~isstruct(output)
    keyboard
end


totalMass = output.launcherMass ; 
tof = output.tof ; 

if isempty(output)
    keyboard
end

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

if any(size(objective)<1)
    keyboard
end






end
