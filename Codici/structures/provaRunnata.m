%% ======================= Prova Runnata ================================== 

clear all;
clc
close all

% ==========================  DATI ========================================

addpath(genpath('..\..\'))

[mission, opt] = dataStruct;

X = load('..\Trajectory\Final Trajectory\InitialPop3stages.mat','X');
X = X.X(:,end);
thrustData(:,:,1) = [X(1:5),X(6:10)];
thrustData(:,:,2) = [X(11:15),X(16:20)];
thrustData(:,:,3) = [X(21:25),X(26:30)];

[timeCollocation, stateCollocation] = totalTrajectory(mission,opt,thrustData,1);

%%
mission.structure.alphaQmax = deg2rad(3.4);

[mission] = externalLoads(timeCollocation,stateCollocation,mission,opt,thrustData, mission.structure.alphaQmax);

if opt.nStages == 1
    mission.structure.nComponents = 8;
elseif opt.nStages == 2
    mission.structure.nComponents = 12;
elseif opt.nStages == 3
    mission.structure.nComponents = 16;
end

mission.structure.nNodes = mission.structure.nComponents + 1;

% Nodes
mission.structure.loadNodes = [2,3:2:mission.structure.nNodes-1,mission.structure.nNodes-1];

% Length of the components of the LV
mission.structure.componentLength = [mission.capsule.height];

for ii = opt.nStages:-1:1
    mission.structure.componentLength = [ mission.structure.componentLength, mission.structures{ii}.lengthInterstage, opt.stage{ii}.length];
end

% Computation of position of xCp and xCp_a (fins)
mission.structure.launcherLength = cumsum(mission.structure.componentLength);
mission.structure.launcherLength = mission.structure.launcherLength(end);

mission.diameter = mission.structure.diameter;
xcp = computeXcp(mission, opt);
xcp_a = mission.aerodynamics.rootChord - computeFinXcp(mission); % compute position starting from the launcher bottom
xcg = computeXCG(mission, opt);

% Length of the element used for structural analysis
mission.structure.elementLength = [mission.capsule.height/2,xcp-mission.capsule.height/2,mission.capsule.height-xcp];

for ii = opt.nStages:-1:1
    mission.structure.elementLength = [ mission.structure.elementLength, mission.structures{ii}.lengthInterstage/2, mission.structures{ii}.lengthInterstage/2, opt.stage{ii}.length/2,opt.stage{ii}.length/2];
end

mission.structure.elementLength(end) = opt.stage{1}.length / 2 - xcp_a;
mission.structure.elementLength(end+1) = xcp_a;

% ====================== CALCOLO AZIONI INTERNE ===========================
[mission] = loadsFinder(mission);

N = mission.structure.N;
T = mission.structure.T;
M = mission.structure.M;


% ==================== SPESSORE E MASSA STRUTTURA =========================

engineUsed = 1;
mission    = thicknessFunction(mission, engineUsed);

thick_mm = mission.structure.thickness * 1e3;   % [mm]
fprintf('Thicknesses of the structures starting from the nose are:\n');
fprintf('  %.3f mm\n', thick_mm);  % stampa un valore per riga
mStruct_ton = mission.structure.mStruct * 1e-3; % [ton]
fprintf('Masses of the structures starting from the nose are:\n');
fprintf('  %.3f tons\n', mStruct_ton);

%% ============================== PLOTS ===================================

nPointsPerComponent = 100;

x_all = [];
N_all = [];
T_all = [];
M_all = [];

x_coordinates = cumsum([0, mission.structure.elementLength]); % defines the coordinates of
% start and finish of each component

% Required for interpolation
M_end_values = [M(2:end); 0];

for i = 1:mission.structure.nComponents
    
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

% --- Figura 1: Axial Load ---
figure(1); 
ar1 = area(x_all, N_all); 
% Proprietà grafiche
ar1.FaceColor = 'b';       % Colore riempimento (blu)
ar1.FaceAlpha = 0.3;       % Trasparenza (0 = invisibile, 1 = solido)
ar1.EdgeColor = 'b';       % Colore della linea superiore
ar1.LineWidth = 1.0;       % Spessore linea (meno spessa di 1.5)

hold on
yline(0, 'LineWidth', 1.5, 'Color', 'k') % Linea zero nera per contrasto
grid on;
xlabel('x [m]');
ylabel('Axial Load [N]');
xlim([0, x_all(end)])
xline([0, cumsum(mission.structure.elementLength)], 'LineStyle','--')

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
xline([0, cumsum(mission.structure.elementLength)], 'LineStyle','--')

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
xline([0, cumsum(mission.structure.elementLength)], 'LineStyle','--')
xline(xcg, 'k',  'LineWidth',1.5)
%%

Ftfins = (- mission.structure.tMaxQ(2) * (mission.structure.launcherLength - xcg) + ( mission.structure.dMaxQ(2) + mission.structure.lMaxQ(2)) * (xcg - xcp)) / (mission.structure.launcherLength - xcp_a - xcg)
CLnew = abs(Ftfins) / (mission.structure.dynamicPressure * mission.aerodynamics.bodyGeom.Aref)
