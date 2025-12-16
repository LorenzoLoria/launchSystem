function [tMax] = unpressurizedTanks(launcher)

% Function required to size the launcher thickess for unpressurized
% tanks. This represents a critical condition when the vehicle is still on
% land. Must be used on components which are not subject to internal
% pressurization

% --- INPUTS
% M = vector containing the mass of the stacks [kg];
% r = radius of the launcher [m];
% h = height of the launcher [m];
% hCM = height of the center of mass of the launcher [m];
% nx = load factor in direction x;
% nz = load factor in direction z;
% g0 = gravity acceleration [ms^-2];
% SF = safety factors;
% sigmaAllowable = yielding stress [Pa];
% E = stiffness [Pa];
% nStages = number of stages;

% --- OUTPUT
% t = thickness [m]
% mStruct = mass of the structure [kg]

M              = launcher.M0;
r              = launcher.R;
h              = launcher.h;
hCM            = launcher.hCM;
nx             = launcher.nx;
nz             = launcher.nz;
g0             = launcher.g0;
SF             = launcher.SF;
sigmaAllowable = launcher.sigmaAllowable;
E              = launcher.E;
shearAllowable = sigmaAllowable / 2;

% --- Solution ------------------------------------------------------------

% Loads
longitudinalLoad = nx .* M * g0; % longitudinal load vector
lateralLoad = nz .* M .* g0; % lateral force vector
bendingMoment = lateralLoad .* hCM; % bending moment vector

% Area A and Inertia I are initially computed using the following
% approximation:
% A = 2 * pi * r * t;
% I = pi * r^3 * t;

% Minimum Allowable Thicknesses
tAxial = abs(- longitudinalLoad / (2 * pi * r) - bendingMoment / (pi * r^2)) / sigmaAllowable * SF;
tShear = lateralLoad / (2 * pi * r * shearAllowable) * SF;

for ii = 1:length(tAxial)
    if tAxial(ii) > tShear(ii)
        t(ii) = tAxial(ii);
    else
        t(ii) = tShear(ii);
    end
end

% Area 
A = pi * (r^2 - (r-t).^2);

% Inertia
I = pi / 4 * (r^4 - (r - t).^4);

% Stresses (MAGNITUDE)
longitudinalStress = longitudinalLoad / A * SF; 
bendingStress = bendingMoment ./ I .* r * SF;
shearStress = lateralLoad ./ A * SF;

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

% Exctraction of the most critical result
[tMax, tMaxPos] = max(t);
% volume = pi * h(1) * (r^2 - (r - tMax)^2);
% mStruct = volume * rhoMaterial;
[bucklingMax, bucklingMaxPos] = max(sigmaBuckling);