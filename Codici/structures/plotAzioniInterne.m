clear
clc
close all

% ============================= DATA ======================================
nComponents = 4;
launcher.mass = ones(nComponents,1) * 50e3;
launcher.drag = ones(nComponents,1) * 50000;
launcher.acceleration = 4 * 9.81;
launcher.alpha = deg2rad(45);
launcher.theta = deg2rad(4);
launcher.thrust = 4*845e3;
launcher.lift= ones(nComponents,1) * 50000;
launcher.g0 = 9.81;

% Launcher single component dimension
launcher.firstStage  = 30;
launcher.interStage  = 2;
launcher.secondStage = 15;
launcher.fairing     = 10;

stagesDimensions = [launcher.firstStage, launcher.interStage, launcher.secondStage, launcher.fairing];

x1 = linspace(0, stagesDimensions(1), 100);
x2 = linspace(stagesDimensions(1), stagesDimensions(1)+stagesDimensions(2), 100);
x3 = linspace(stagesDimensions(1)+stagesDimensions(2), stagesDimensions(1)+stagesDimensions(2)+stagesDimensions(3), 100);
x4 = linspace(stagesDimensions(1)+stagesDimensions(2)+stagesDimensions(3), stagesDimensions(1)+stagesDimensions(2)+stagesDimensions(3)+stagesDimensions(4), 100);

% ============================ SOLUTION ===================================
loads = @(x) loadsFinder(nComponents, x, launcher);
loadsResults = loads(stagesDimensions);

% Axial Loads
N = loadsResults(1:3:end);
% Shear Loads
T = loadsResults(2:3:end);
% Bending Moment Loads
M = loadsResults(3:3:end);

% ============================== PLOTS ====================================

N1 = ones(100,1) * N(1);
N2 = ones(100,1) * N(2);
N3 = ones(100,1) * N(3);
N4 = ones(100,1) * N(4);

T1 = ones(100,1) * T(1);
T2 = ones(100,1) * T(2);
T3 = ones(100,1) * T(3);
T4 = ones(100,1) * T(4);

M1 = linspace(M(1),M(2), 100);
M2 = linspace(M(2),M(3), 100);
M3 = linspace(M(3),M(4), 100);
M4 = linspace(M(4),0, 100);

x_all = [x1, x2, x3, x4];
N_all = [N1.', N2.', N3.', N4.'];   % qui metto il .' per avere riga
T_all = [T1.', T2.', T3.', T4.'];
M_all = [M1, M2, M3, M4];


figure; plot(x_all, N_all, 'LineWidth', 1.5);
grid on;
xlabel('x');
ylabel('N');
title('N(x)');

figure; plot(x_all, T_all, 'LineWidth', 1.5);
grid on;
xlabel('x');
ylabel('T');
title('T');

figure; plot(x_all, M_all, 'LineWidth', 1.5);
grid on;
xlabel('x');
ylabel('M');
title('M(x)');