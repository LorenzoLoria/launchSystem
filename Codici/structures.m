function [structuralMass] = structures(launcher)

% Function required to size the launcher thickness

% --- INPUTS
% M = total mass of the launcher;
% r = radius of the launcher;
% h = height of the launcher;
% hCM = height of the center of mass of the launcher;
% nx = load factor in direction x;
% nz = load factor in direction z;
% g0 = gravity acceleration;
% SF = safety factors;
% sigmaAllowable = maximum allowable stress;
% E = stiffness;
% p = tank pressure;
% rho = density of the propellant;
% rhoMaterial = density of the material chosen for the structure

% --- OUTPUT
% t = thickness
% sigmaCritical = critical axial stress
% maxStress = maximum axial stress

M = launcher.M0;
r = launcher.R;
h = launcher.h;
hCM = launcher.hCM;
nx = launcher.nx;
nz = launcher.nz;
g0 = launcher.g0;
SF = launcher.SF;
sigmaAllowable = launcher.sigmaAllowable;
E = launcher.E;
p = launcher.p;
rho = launcher.rho;
rhoMaterial = launcher.rhoMaterial;

% --- Solution ------------------------------------------------------------
% Loads
P = nx * M * g0; % axial load
Fb = nz * M * g0; % lateral force
bendingMoment = Fb * hCM; % bending moment
pressureLoad = pi * r^2 * p; % pressure load

% Shear Stress
shearAllowable = sigmaAllowable / 2; % typical value

% Minimum Allowable Thickness
tAxialUnpressurized = ((P / (2 * pi * r)) + (bendingMoment / (pi * r^2))) / sigmaAllowable * SF;
tAxialPressurized = abs(- P / (2 * pi * r) - bendingMoment / (pi * r^2) + p * R / 2 + rho * nx * g0 * r / 2) / sigmaAllowable * SF;
tShear = Fb / (2 * pi * r * shearAllowable) * SF;

% Choice of the maximum thickness (represents heaviest scenario)
if tAxialUnpressurized > tShear && tAxialUnpressurized > tAxialPressurized
    t = tAxialUnpressurized;
elseif tAxialPressurized > tShear && tAxialPressurized > tAxialUnpressurized
    t = tAxialPressurized;
else
    t = tShear;
end

% Area 
A = pi * (r^2 - (r - t)^2);

% Inertia
I = pi * r^3 * t;

% Structural Mass Computation
volume = pi * (r^2 - (r - t)^2) * h;
structuralMass = volume * rhoMaterial;



% Critical Axial Stress Computation
sigmaCriticalUnpressurized = E * (9 * (t / r)^1.6 + 0.16 * (t / h)^1.3);

k0 = 9 * (t / r)^0.6 + 0.16 * (r / h)^1.3 * (t / r)^0.3;
kp = 0.191 * (p / E) * (r / t)^2;
if kp < 0.229
else
    kp = 0.229;
end

sigmaCriticalPressurized = ((k0 + kp) * E * t) / r;

if sigmaCriticalPressurized > sigmaCriticalUnpressurized
    sigmaCritical = sigmaCriticalPressurized;
else
    sigmaCritical = sigmaCriticalUnpressurized;
end

% Hydrostatic pressure
hydrostaticPressure = rho * loadFactorX * g0 * height;

% Stresses
sigmaAxial = P / A; % - sign since it is compressing (+ if it is expanding)
sigmaBending = bendingMoment / I * r;
sigmaPressure = p * r / (2 * t);
sigmaHydro = rho * nx * g0 * h * r / (2 * t);

maxStressUnpressurized = - abs(sigmaAxial) - abs(sigmaBending); % sign is because of compression
maxStressPressurized = - abs(sigmaAxial) - abs(sigmaBending) + abs(sigmaPressure) + abs(sigmaHydro);

