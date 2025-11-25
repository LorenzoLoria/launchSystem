function [t, tMax, stressMatrix, mStruct] = pressurizedTanks(mission, nComponents)

% Function required to size the launcher thickess when tanks are
% pressurized during the flight. 

% --- INPUTS
% mission = struct containing LV data
% nComponents = number of components in which the LV is discretized

% --- OUTPUT
% t = vector containing thickness of each component [m]
% tMax = maximum thickness of the launcher [m]
% stressMatrix = (6 x nComponents) matrix. Each column represents a
% component of the LV (e.g. Thrust structure, 1st stage etc.), each row
% represents a different kind of stress
% mStruct = vector containing mass of each structure [kg]

% ============================== DATA =====================================
M              = mission.launcher.mass; % total mass of the element considered [kg]
r              = mission.launcher.diameter / 2; % radius of cylindrical section [m]
h              = mission.launcher.length; % length of the cylinder [m]
hCM            = mission.launcher.hCM; % position of the center of mass [m]
nx             = launcher.nx; % longitudinal load factor
nz             = launcher.nz; % lateral load factor
g0             = mission.environment.g0; 
SF             = mission.launcher.structures.SF; % safety factor
ultimateStress = mission.launcher.structures{ii}.ultimate; % ultimate stress for chosen material [Pa]
E              = mission.launcher.structures{ii}.E; % Young Modulus for chosen material [Pa]
shearAllowable = ultimateStress / 2; % ultimate shear stress for chosen material [Pa]
rhoOX          = mission.launcher.engines{ii}.oxDens; % density of the oxidizer [kg/m^3]
rhoFu          = mission.launcher.engines{ii}.fuDens; % density of the fuel [kg/m^3]
rhoMaterial    = mission.launcher.structures{ii}.rho; % density of the chosen material [kg/m^3] 
p              = mission.launcher.tankPressure; % pressure of tanks [kg/m^3]
longitudinalLoad = mission.launcher.structures.N; % Axial Load [N] (from loadFinder)
lateralLoad    = mission.launcher.structures.T; % Shear Load [N] (from loadFinder)
bendingMoment  = mission.launcher.structures.Mb; % Bending Moment [Nm] (from loadFinder)
% ============================ Solution ===================================

% % Loads
% longitudinalLoad = nx .* M * g0; % longitudinal load vector
% lateralLoad = nz .* M .* g0; % lateral force vector
% bendingMoment = lateralLoad .* hCM; % bending moment vector

% Minimum Allowable Thicknesses
tAxial = abs(- longitudinalLoad / (2 * pi * r) - bendingMoment / (pi * r^2) + p * r / 2 + rhoOX * nx * g0 * r / 2) / ultimateStress * SF;
tShear = lateralLoad / (2 * pi * r * shearAllowable) * SF;

for ii = 1:length(tAxial)
    if tAxial(ii) > tShear(ii)
        t(ii) = tAxial(ii);
    else
        t(ii) = tShear(ii);
    end
end

% Hydrostatic pressure
pHydro = rhoOX .* nx .* g0 .* h;

% Area 
A = pi * (r^2 - (r-t).^2);

% Inertia (valid for a cylinder)
I = pi / 4 * (r^4 - (r - t).^4);

% Mass computation
volume = pi .* h .* (r^2 - (r - t).^2);
mStruct = volume .* rhoMaterial;
[tMax] = max(t);

% Buckling for Pressurized Cylinders
k0 = 9 * (t / r).^(0.6) + 0.16 .* (r / h).^(1.3) .* (t / r).^(0.3);
kp = 0.191 * (p / E) .* (r / t).^2;

for ii = 1:length(kp)
    if kp(ii) < 0.229
    else
        kp(ii) = 0.229;
    end
end

sigmaBuckling = ((k0 + kp) * E .* t) / r;

% Stresses (MAGNITUDE)
longitudinalStress = longitudinalLoad / A * SF; 
bendingStress = bendingMoment ./ I .* r * SF;
shearStress = lateralLoad ./ A * SF;
hoopStress = (p .* r + rhoOX .* g0 .* nx .* h .* r) ./ t;
axialStress = hoopStress / 2;

stressMatrix = [longitudinalStress, bendingStress, shearStress, hoopStress, axialStress, sigmaBuckling];

end





% --- Data
% mass = 50000;
% height = 8;
% hCM = height / 2;
% radius = 2.5 / 2;
% thickness = 1e-3;
% loadFactorX = 2.5; % depends on the most critical load condition
% loadFactorZ = 0.2; % depends on the most critical load condition
% g0 = 9.81;
% safeyFactor = 1.25; % typically between 1.1 and 1.5
% sigmaAllowableSteel = 1034e6;
% sigmaAllowableAl = 448e6;
% shearAllowableSteel = sigmaAllowableSteel / 2;
% shearAllowableAl = sigmaAllowableAl / 2;
% stiffnessSteel = 207e9;
% stiffnessAl = 69e9;