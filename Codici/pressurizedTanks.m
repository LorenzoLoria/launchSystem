function [sigmaCritical] = pressurizedTanks(M, r, h, hCM, nx, nz, g0, SF, sigmaAllowable, E, p, rho)

% Function required to size the launcher thickess

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

% --- OUTPUT
% t = thickness

% --- Solution ------------------------------------------------------------
% Loads
P = nx * M * g0; % axial load
Fb = nz * M * g0; % lateral force
bendingMoment = Fb * hCM; % bending moment
Faxial = pi * r^2 * p; % pressure load

% Minimum Allowable Thickness
tAxial = abs(- P / (2 * pi * r) - bendingMoment / (pi * r^2) + p * R / 2 + rho * nx * g0 * r / 2) / sigmaAllowable * SF;
tShear = Fb / (2 * pi * r * shearAllowable) * SF;

if tAxial > tShear
    t = tAxial;
else
    t = tShear;
end

% Critical Stress for Pressurized Cylinders

k0 = 9 * (t / r)^0.6 + 0.16 * (r / h)^1.3 * (t / r)^0.3;
kp = 0.191 * (p / E) * (r / t)^2;
if kp < 0.229
else
    kp = 0.229;
end

sigmaCritical = ((k0 + kp) * E * t) / r;

% Hydrostatic pressure
pHydro = rho * loadFactorX * g0 * height;

% Area 
A = 2 * pi * r * t;

% Inertia
I = pi * r^3 * t;

% Stresses
sigmaAxial = - axial / A * safeyFactor; % - sign since it is compressing (+ if it is expanding)
sigmaBending = bendingMoment / I * radius * safeyFactor;


maxStress = abs(sigmaAxial) + abs(sigmaBending);

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