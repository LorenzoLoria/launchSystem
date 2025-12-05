function [cin,ceq] = nlconGlobalGA(launcher,mission,opt,settings,option2D)

ceq = [];

nlconFlag = 1; 

[output] = launcherSimulation(launcher,mission,opt,settings,option2D,nlconFlag);

cin = [output];





