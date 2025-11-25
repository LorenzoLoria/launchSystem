clear;
clc;

mission.launcher.diameter = 4; 
mission.launcher.length = 30; 
mission.launcher.nx = 4;           
mission.environment.g0 = 9.81; 
mission.launcher.structures.SF = 1.5; 
mission.launcher.structures.ultimate = 480e6; 
mission.launcher.structures.E = 72e9;
mission.launcher.engines.oxDens = 1143;
mission.launcher.structures.rho = 2840; 
mission.launcher.tankPressure = 50000; 
mission.launcher.structures.N = 845e3;
mission.launcher.structures.T = 2e3; 
mission.launcher.structures.Mb = 2e4; 

[massStruct, t, tMax, stressMatrix] = thicknessFunction(mission, 1)