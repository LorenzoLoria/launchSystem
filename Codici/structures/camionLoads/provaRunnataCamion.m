%% ================== Prova Runnata Sul Camion ============================ 

clear all;
clc;
close all

% ==========================  DATI ========================================

addpath(genpath('..\..\'))

[mission, settings] = dataStructGlobal;

mission.structure.alphaQmax = 0;

% launcher = [nStages, nMotore1, nMotore2, nMotore3, %massa1, %massa2,
% %massa3];
launcher = [2,2,3,4,0.459952176990556, 0.753370531158904, 0.634795741885559];

for i = 1:launcher(1)
    configuration.stage{i}.engine = mission.engines{launcher(1+i)};
end

[mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);
thrustDataGA = load('thrustdataVecTraj.mat','xGATraj');

thrustDataVecFMC(:,:,1) = [0.902082365568723	1.480898931628005
                            0.999984156345040	23.253294859564580
                            0.900002678098914	52.979241033086943
                            0.900000000000007	59.571815331701984
                            0.903941814015555	55.058714159781090];


thrustDataVecFMC(:,:,2 ) = [0.400917809388214	65.122710138507202
                            0.964494359624014	79.658359202140389
                            0.975968800776448	91.801043507018605
                            0.992714640706230	89.085172454227390
                            0.993244065056187	99.345740209598944];



[timeCollocation, stateCollocation] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustDataVecFMC);

%% ============================ SOLUTION ==================================

[landLoads] = camionLoads(mission, configuration, mer, launcher, staging);

N = landLoads.N;
T = landLoads.T;
M = landLoads.M;

%% ============================== PLOTS ===================================

if launcher(1) == 1
    nElements = 14;
elseif launcher(1) == 2
    nElements = 18;
elseif launcher(1) == 3
    nElements = 22;
end

% Length of the element used for structural analysis
h = [2 / 3 * mission.capsule.height, 1 /3 * mission.capsule.height]; % centroid of cone

for ii = launcher(1):-1:1
    h = [ h, configuration.geometry.stage{ii}.interstage.length/2, ...
    configuration.geometry.stage{ii}.interstage.length/2, configuration.stage{ii}.fuelTankH/2, ...
    configuration.stage{ii}.fuelTankH/2, configuration.stage{ii}.oxTankH/2, ...
    configuration.stage{ii}.oxTankH/2, configuration.geometry.stage{ii}.thrustFrame/2,...
    configuration.geometry.stage{ii}.thrustFrame/2];
end

nPointsPerComponent = 100;

x_all = [];
N_all = [];
T_all = [];
M_all = [];

x_coordinates = cumsum([0, h]); % defines the coordinates of
% start and finish of each component

% Required for interpolation
M_end_values = [M(2:end); 0];

for i = 1:nElements
    
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

% --- Figura 2: Shear Load ---
figure(2); 
ar2 = area(x_all, T_all);
ar2.FaceColor = 'b';
ar2.FaceAlpha = 0.3;
ar2.EdgeColor = 'b';
ar2.LineWidth = 1.0;

hold on
yline(0, 'LineWidth', 1.5, 'Color', 'k')
grid on;
xlabel('x [m]');
ylabel('Shear Load [N]');
xlim([0, x_all(end)])
% xline([0, cumsum(h)], 'LineStyle','--')

% --- Figura 3: Bending Moment ---
figure(3); 
ar3 = area(x_all, M_all);
ar3.FaceColor = 'b';
ar3.FaceAlpha = 0.3;
ar3.EdgeColor = 'b';
ar3.LineWidth = 1.0;

hold on
yline(0, 'LineWidth', 1.5, 'Color', 'k')
grid on;
xlabel('x [m]');
ylabel('Bending Moment [Nm]');
xlim([0, x_all(end)])