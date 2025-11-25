function [mer] = MER(mission)

% Computation of the Mass Estimation Relations
% --- INPUTS
% mission = structure containing all the launcher data
% --- OUTPUTS
% mission.launcher.mer = structure containing all Mass Estimation Relations results

mp = mission.launcher.mPropellant1; % mass of propellant [kg]
OF = mission.launcher.engines{1}.OF; % OF ratio
rhoOx = mission.launcher.engines{1}.oxDens; % density of oxidizer [kg/m^3]
rhoFu = mission.launcher.engines{1}.fuelDens; % density of fuel [kg/m^3]
areaTankLOX = mission.launcher.areaTankOx; % area of the LOX tank [m^2]
areaTankFu = mission.launcher.areaTankFu; % area of the fu tank [m^2]
A = 1.3; % ranges between 1.3-2.6 ---> for turbopump mass
b = 0.6; % ranges between 0.6-0.666 ---> for turbopump mass
pumpRotationalSpeed = mission.launcher.engines{1}.pumpRotationalSpeed; 
requiredPower = mission.launcher.engines{1}.requiredPower;  

volumeTankLOX = mp * OF / (OF + 1) / rhoOx * 1.055; % volume occupied by the ox [m^3] w/ margin
volumeTankFu = mp / (OF + 1) / rhoFu * 1.055; % volume occupied by the fuel [m^3] w/ margin 

% margins : 3% for ullage, 0.5% for residuals, 2% for loading margins

% ============================ SOLUTION ===================================

% LOX Tank Mass
mer.tankMassLOX = 12.2 * volumeTankLOX + 255.2;

% LH2 Tank Mass
mer.tankMassLH2 = 9.08 * volumeTankFu + 100.09;

% LOX Cryogenic Insulation Mass 
mer.cryoInsuLOX = 1.12 * areaTankLOX;

% LH2 Cryogenic Insulation Mass 
mer.cryoInsuFu = 2.88 * areaTankFu;

% Turbopumps Mass
mer.turbopumps = A * (requiredPower / pumpRotationalSpeed)^b;

% Fairing and Shroud Mass
mer.fairing = 4.95 * fairingArea^1.15;

% Avionics Mass
mer.avionics = 10 * (M0)^0.361;

% Wiring Mass
mer.wiring = 1.058 * sqrt(M0)^0.25;

% Thrust structure
mer.thrustStructure = 2.55e-4 * thrust;

% Mass of the engine
mer.Engine = 150 + 0.086 * (thrustVacuum)^0.86; % value in lb
mer.Engine = mer.Engine * 0.45359237; % value in kg