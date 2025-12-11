function [cin,ceq] = nlconGlobalGA(launcher,mission,settings)

ceq = [];

nlconFlag = 1; 

[~,output] = launcherSimulation(launcher,mission,settings,nlconFlag);

if isnan(output)
    keyboard
end
if isinf(output)
    keyboard
end
if ~isreal(output)
    keyboard
end


cin = [output];





