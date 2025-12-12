function  [objective] = objFunGlobalGA(launcher,mission,settings)

nlconFlag = 0;

[outputOBJ] = launcherSimulation(launcher,mission,settings,nlconFlag);

if ~isstruct(outputOBJ)
    keyboard
end


totalMass = outputOBJ.launcherMass ; 
tof = outputOBJ.tof ; 

if isempty(outputOBJ)
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
