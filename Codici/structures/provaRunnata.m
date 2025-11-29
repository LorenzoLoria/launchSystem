%% ======================= Prova Runnata ================================== 

clear all;
clc

% ==========================  DATI ========================================
[mission, opt] = dataStruct;

thrustData(:,:,1) =     [0.9956   -2.7415   13.0559; ...
    0.9243   -2.3127    4.7079; ...
    0.9720   -2.4832   27.0284; ...
    0.9607   13.3110   22.9276; ...
    0.9868  -33.8890   83.4509];
thrustData(:,:,2) =    [0.7495   26.9887   29.4877; ...
    0.7732  -22.4865   93.7944; ...
    0.8558   42.2771   77.4467; ...
    0.8221   30.4242   77.8149; ...
    0.4019   35.0922   87.7839];
[timeCollocation, stateCollocation] = totalTrajectory(mission,opt,thrustData);

[q,an,at,F,D,angle,gamma] = externalLoads(timeCollocation,stateCollocation,mission,opt,thrustData);

nComponents = 8;
nPointsPerComponent = 100;
launcher.mass = [mission.capsule.weigth, opt.stage{2}.mStage, opt.stage{1}.mStage];
launcher.drag = D;
launcher.lift = 0;
launcher.accelerationAxial = an;
launcher.accelerationNormal = at;
launcher.gamma = gamma;
launcher.alpha = 0;
launcher.g0 = mission.environment.g0;
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