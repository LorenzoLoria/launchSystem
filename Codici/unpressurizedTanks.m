function [sigmaCritical, maxStress] = unpressurizedTanks(M, r, h, hCM, nx, nz, g0, SF, sigmaAllowable, E)

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

% --- OUTPUT
% t = thickness

% --- Solution ------------------------------------------------------------
% Loads
P = nx * M * g0; % axial load
Fb = nz * M * g0; % lateral force
bendingMoment = Fb * hCM; % bending moment

% Minimum Allowable Thickness
tAxial = ((P / (2 * pi * r)) + (bendingMoment / (pi * r^2))) / sigmaAllowable * SF;
tShear = Fb / (2 * pi * r * shearAllowable) * SF;

if tAxial > tShear
    t = tAxial;
else
    t = tShear;
end


% Critical Axial Stress Computation (Empirical Correlation) 
sigmaCritical = E * (9 * (t / r)^1.6 + 0.16 * (t / h)^1.3);

% Area 
A = 2 * pi * r * t;

% Inertia
I = pi * r^3 * t;

% Stresses
sigmaAxial = - axial / A * safeyFactor; % - sign since it is compressing (+ if it is expanding)
sigmaBending = bendingMoment / I * radius * safeyFactor;

maxStress = abs(sigmaAxial) + abs(sigmaBending);

% shearStress = lateralForce / A * safeyFactor;
% percentage = maxStress / sigmaAllowable * 100;