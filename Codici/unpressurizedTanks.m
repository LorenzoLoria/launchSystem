function [sigmaCritical, maxStress] = unpressurizedTanks(launcher)

% Function required to size the launcher thickess for unpressurized
% tanks. This represents a critical condition when the vehicle is still on
% land

% --- INPUTS
% M = total mass of the launcher [kg];
% r = radius of the launcher [m];
% h = height of the launcher [m];
% hCM = height of the center of mass of the launcher [m];
% nx = load factor in direction x;
% nz = load factor in direction z;
% g0 = gravity acceleration [ms^-2];
% SF = safety factors;
% sigmaAllowable = maximum allowable stress [Pa];
% E = stiffness [Pa];
% nStages = number of stages;

% --- OUTPUT
% t = thickness [m]
% mStruct = mass of the structure [kg]

M              = launcher.M;
r              = launcher.R;
h              = launcher.h;
hCM            = launcher.hCM;
nx             = launcher.nx;
nz             = launcher.nz;
g0             = launcher.g0;
SF             = launcher.SF;
sigmaAllowable = launcher.sigmaAllowable;
E              = launcher.E;
nStages        = launcher.nStages;

% --- Solution ------------------------------------------------------------

for i = 1:nStages
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