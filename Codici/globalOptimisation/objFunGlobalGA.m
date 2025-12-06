function  [objective] = objFunGlobalGA(launcher,mission,opt,settings,option2D)

nlconFlag = 0;

[output] = launcherSimulation(launcher,mission,opt,settings,option2D,nlconFlag);


    output.launcherMass = totalMass;
    output.tof = fvalFMCTraj; 

objective = output.totalMass/200000 + output.tof/3600 ; 











end
