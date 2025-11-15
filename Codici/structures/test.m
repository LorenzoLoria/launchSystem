clear
clc

% --- Data
launcher.M = 50000;
launcher.h = 8;
launcher.hCM = 4;
launcher.R = 2.5 / 2;
launcher.nx = 2.5; % depends on the most critical load condition
launcher.nz = 0.2; % depends on the most critical load condition
g0 = 9.81;
launcher.SF = 1.25; % typically between 1.1 and 1.5
launcher.sigmaAllowable = 1034e6;
launcher.shearAllowable = launcher.sigmaAllowable / 2;
launcher.E = 207e9;
launcher.rhoProp = 1400;
launcher.rhoMaterial = 0;
launcher.tankPressure = 50000;
t = pressurizedTanks(launcher)