function [cin,ceq] = nlconGlobalGA(launcher,mission,settings)

ceq = [];

nlconFlag = 1; 

[~,outputNLC] = launcherSimulation(launcher,mission,settings,nlconFlag);

if isnan(outputNLC)
    keyboard
end
if isinf(outputNLC)
    keyboard
end
if ~isreal(outputNLC)
    keyboard
end


cin = [outputNLC];





