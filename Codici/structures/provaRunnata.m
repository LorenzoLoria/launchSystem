%% ======================= Prova Runnata ================================== 

clear all;
clc
close all

% ==========================  DATI ========================================
[mission, opt] = dataStruct;

thrustData(:, :, 1) =     [0.992951743036793	26.174503353658913
0.918245484987595	5.524888452757455
0.949721966745758	4.816713623391522
0.917770387441014	73.262523254569814
0.962817207168307	42.233780805019514];
thrustData(:, :, 2) =    [0.609126028473365	8.838970511617315
0.949237174404871	70.435231583491500
0.702233911900803	1.003603356352884e+02
0.845617612957625	52.146897489107950
0.413941032471940	99.388409631182256];
    
[timeCollocation, stateCollocation] = totalTrajectory(mission,opt,thrustData,1);
alpha = deg2rad(0);

[q,dMaxQ,lMaxQ,aMaxQ,tMaxQ, massMaxQ, g] = externalLoads(timeCollocation,stateCollocation,mission,opt,thrustData, alpha);

nComponents = 8;
nPointsPerComponent = 100;
m1Stage = massMaxQ - mission.capsule.weigth - opt.stage{2}.mStage;
launcher.mass = [mission.capsule.weigth, opt.stage{2}.mStage, m1Stage];
launcher.drag = dMaxQ;
launcher.lift = lMaxQ;
launcher.acceleration = aMaxQ;
launcher.g0 = g;
launcher.elementLength = [4,2,3,10,11,10,3,12];
launcher.dragFins = [0 0];
launcher.liftFins = [0 0];

[N, T, M, A] = loadsFinder(nComponents, launcher);

%% ============================== PLOTS ====================================

x_all = [];
N_all = [];
T_all = [];
M_all = [];

x_coordinates = cumsum([0, launcher.elementLength]); % defines the coordinates of
% start and finish of each component

% Required for interpolation
M_end_values = [M(2:end); 0];

for i = 1:nComponents
    
    x_start = x_coordinates(i);
    x_end   = x_coordinates(i+1);
    x_current = linspace(x_start, x_end, nPointsPerComponent);
    
    % Axial Load
    N_current = ones(1, nPointsPerComponent) * N(i);
    
    % Shear Load
    T_current = ones(1, nPointsPerComponent) * T(i);
    
    % Bending Moment (changes linearly)
    M_start = M(i);
    M_end   = M_end_values(i);
    M_current = linspace(M_start, M_end, nPointsPerComponent);
    
    x_all = [x_all, x_current];
    N_all = [N_all, N_current];
    T_all = [T_all, T_current];
    M_all = [M_all, M_current];
end

figure(1); 
plot(x_all, N_all, 'LineWidth', 1.5, 'Color', 'b');
hold on
yline(0, 'LineWidth', 1.5)
grid on;
xlabel('x [m]');
ylabel('Axial Load [N]');
xlim([0,x_all(end)])

figure(2); 
plot(x_all, T_all, 'LineWidth', 1.5, 'Color', 'b');
hold on
yline(0, 'LineWidth', 1.5)
grid on;
xlabel('x [m]');
ylabel('Shear Load [N]');
xlim([0,x_all(end)])

figure; 
plot(x_all, M_all, 'LineWidth', 1.5, 'Color', 'b');
hold on
yline(0, 'LineWidth', 1)
grid on;
xlabel('x [m]');
ylabel('Bending Moment [Nm]');
xlim([0,x_all(end)])