%% Structural Analysis
%% Unpressurized tank
clear;
clc;
close all;

% --- Data
mass = 50000;
height = 8;
hCM = height / 2;
radius = 2.5 / 2;
thickness = 1e-3;
loadFactorX = 2.5; % depends on the most critical load condition
loadFactorZ = 0.2; % depends on the most critical load condition
g0 = 9.81;
safeyFactor = 1.25; % typically between 1.1 and 1.5
sigmaAllowableSteel = 1034e6;
sigmaAllowableAl = 448e6;
shearAllowableSteel = sigmaAllowableSteel / 2;
shearAllowableAl = sigmaAllowableAl / 2;
stiffnessSteel = 207e9;
stiffnessAl = 69e9;

% --- Solution

% Cylinder Surface (approximated formulation)
A = 2 * pi * radius * thickness;

% Inertia of the main body (approximated formulation)
I = pi * radius^3 * thickness;

% Loads
axial = loadFactorX * mass * g0;
lateralForce = loadFactorZ * mass * g0;
bendingMoment = lateralForce * hCM;

% Stresses
sigmaAxial = - axial / A * safeyFactor; % - sign since it is compressing (+ if it is expanding)
sigmaBending = bendingMoment / I * radius * safeyFactor;

maxStress = abs(sigmaAxial) + abs(sigmaBending);
% fprintf('The maximum stress is: %f MPa \n', maxStress*1e-6)

shearStress = lateralForce / A * safeyFactor;
% fprintf('The shear stress is: %f MPa \n', shearStress*1e-6)

% Comparison between materials
percentageSteel = maxStress / sigmaAllowableSteel * 100;
percentageAl = maxStress / sigmaAllowableAl * 100;

% Minimum Allowable Thickness
tSteelSigma = ((axial / (2 * pi * radius)) + (bendingMoment / (pi * radius^2))) / sigmaAllowableSteel;
tAlSigma = ((axial / (2 * pi * radius)) + (bendingMoment / (pi * radius^2))) / sigmaAllowableAl;

tSteelShear = lateralForce / (2 * pi * radius * shearAllowableSteel);
tAlShear = lateralForce / (2 * pi * radius * shearAllowableAl);

% Critical Axial Stress Computation (Empirical Correlation) 
criticalSigmaSteel = stiffnessSteel * (9 * (tSteelSigma / radius)^1.6 + 0.16 * (tSteelSigma / height)^1.3);
criticalSigmaAl = stiffnessAl * (9 * (tAlSigma / radius)^1.6 + 0.16 * (tAlSigma / height)^1.3);
% Note: the minimum thickness has been used, even though the one of the
% main body for steel is actually difficult to manufacture

%% Pressurized Tank
clear;
clc;
close all;

% --- Data
mass = 50000;
height = 8;
hCM = height / 2;
radius = 2.5 / 2;
thickness = 1e-3;
loadFactorX = 2.5; % depends on the most critical load condition
loadFactorZ = 0.2; % depends on the most critical load condition
g0 = 9.81;
safeyFactor = 1.25; % typically between 1.1 and 1.5
pressure = 50 / 14.504 *1e5;
sigmaAllowableSteel = 1034e6;
sigmaAllowableAl = 448e6;
shearAllowableSteel = sigmaAllowableSteel / 2;
shearAllowableAl = sigmaAllowableAl / 2;
stiffnessSteel = 207e9;
stiffnessAl = 69e9;


% --- Solution
% Cylinder Surface (approximated formulation)
A = 2 * pi * radius * thickness;

% Inertia of the main body (approximated formulation)
I = pi * radius^3 * thickness;

% Loads
axial = loadFactorX * mass * g0;
lateralForce = loadFactorZ * mass * g0;
bendingMoment = lateralForce * hCM;
axialForce = pi * radius^2 * pressure;

% Stresses
sigmaAxial = - axial / A * safeyFactor; % - sign since it is compressing (+ if it is expanding)
sigmaBending = bendingMoment / I * radius * safeyFactor;

% Critical Stress for Pressurized Cylinders
sigmaCritical;