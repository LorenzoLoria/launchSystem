%% ======================= Prova Runnata ================================== 

clear all;
clc
close all

% ==========================  DATI ========================================
[mission, opt] = dataStruct;

thrustData(:,:,1) =     [0.9880  -17.5905    3.0108
    0.9617   24.7776    6.2794
    0.9691    2.4299    9.6882
    0.9000   -8.2757   34.2836
    0.9253    5.9692   59.8887];
thrustData(:,:,2) =    [0.6357   80.9894   44.2025
    0.6141   -5.6012   77.1215
    0.8019   30.9378   72.1261
    0.6585   10.2303   85.1640
    0.3508   29.7187   66.9287];
    
[timeCollocation, stateCollocation] = totalTrajectory(mission,opt,thrustData);

[q,aCC_qMax,F,D,angle,gamma, mQmax,g, vQmax] = externalLoads(timeCollocation,stateCollocation,mission,opt,thrustData);

nComponents = 8;
nPointsPerComponent = 100;
m1Stage = mQmax - mission.capsule.weigth - opt.stage{2}.mStage;
launcher.mass = [mission.capsule.weigth, opt.stage{2}.mStage, m1Stage];
launcher.drag = D;
launcher.lift = 0;
% launcher.accelerationAxial = an;
% launcher.accelerationNormal = at;
launcher.gamma = gamma;
launcher.alpha = 0;
launcher.g0 = g;
launcher.elementLength = [4,2,3,10,11,10,3,12];
launcher.dragFins = 0;
launcher.liftFins = 0;

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