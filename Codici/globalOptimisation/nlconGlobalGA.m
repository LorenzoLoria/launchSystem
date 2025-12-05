function [cin,ceq] = nlconGlobalGA(launcher,mission,settings)

ceq = [];

nlconFlag = 1; 

[output] = launcherSimulation(launcher,mission,settings,nlconFlag);

cin = [output];





