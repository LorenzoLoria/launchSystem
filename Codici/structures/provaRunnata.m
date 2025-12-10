%% ======================= Prova Runnata ================================== 

clear all;
clc
close all

% ==========================  DATI ========================================

addpath(genpath('..\..\'))

[mission, settings] = dataStructGlobal;

mission.structure.alphaQmax = 0;

% launcher = [nStages, nMotore1, nMotore2, nMotore3, %massa1, %massa2,
% %massa3];
launcher = [2,1,3,4,0.5605,0.5595,0.7];

for i = 1:launcher(1)
    configuration.stage{i}.engine = mission.engines{launcher(1+i)};
end

[mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);
thrustDataGA = load('thrustdataVecTraj.mat','xGATraj');

thrustData(:,:,1) = [thrustDataGA.xGATraj(1:5)',thrustDataGA.xGATraj(6:10)'];
thrustData(:,:,2) = [thrustDataGA.xGATraj(11:15)',thrustDataGA.xGATraj(16:20)'];



[timeCollocation, stateCollocation] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustData);

%%

[mission] = externalLoads(timeCollocation,stateCollocation,mission,configuration,thrustData, launcher, mer, 0);

if launcher(1) == 1
    nElements = 8;
elseif launcher(1) == 2
    nElements = 12;
elseif launcher(1) == 3
    nElements = 16;
end

nNodes = nElements + 1;

% Nodes
loadNodes = [2,3:2:nNodes-1,nNodes-1];

% Length of the components of the LV
mission.structure.componentLength = [mission.capsule.height];

for ii = launcher(1):-1:1
    if ii == 1
        mission.structure.componentLength = [ mission.structure.componentLength, configuration.geometry.stage{ii}.interstage.length, configuration.geometry.stage{ii}.length - configuration.stage{1}.engine.length];
    else
    mission.structure.componentLength = [ mission.structure.componentLength, configuration.geometry.stage{ii}.interstage.length, configuration.geometry.stage{ii}.length];
    end
end

% Computation of position of xCp and xCp_a (fins)
mission.structure.launcherLength = cumsum(mission.structure.componentLength);
mission.structure.launcherLength = mission.structure.launcherLength(end);

%mission.diameter = mission.structure.diameter;
xcp = computeXcp(mission, configuration,launcher);
xcp_a = mission.aerodynamics.rootChord - computeFinXcp(mission);
% xcg = computeXCG(mission, configuration, launcher, mer);

% Length of the element used for structural analysis
mission.structure.elementLength = [mission.capsule.height/2, xcp - mission.capsule.height/2,mission.capsule.height - xcp];

for ii = launcher(1):-1:1
    mission.structure.elementLength = [ mission.structure.elementLength, configuration.geometry.stage{ii}.interstage.length/2, configuration.geometry.stage{ii}.interstage.length/2, configuration.geometry.stage{ii}.length/2,configuration.geometry.stage{ii}.length/2];
end

mission.structure.elementLength(end) = configuration.geometry.stage{1}.length/2 - xcp_a;
mission.structure.elementLength(end+1) = xcp_a;

% ====================== CALCOLO AZIONI INTERNE ===========================
% mission.structure.Ftfins = (- mission.structure.tMaxQ(2) * (mission.structure.launcherLength - xcg) + ( mission.structure.dMaxQ(2) + mission.structure.lMaxQ(2)) * (xcg - xcp)) / (mission.structure.launcherLength - xcp_a - xcg);
% clAlpha = 4;
% alphaMax = 1/360*2*pi;
% 
% newAreaFins = abs(mission.structure.Ftfins / (mission.structure.dynamicPressure * clAlpha * alphaMax) );

[mission] = loadsFinder(mission, nElements, loadNodes);

N = mission.structure.N;
T = mission.structure.T;
M = mission.structure.M;

% ==================== SPESSORE E MASSA STRUTTURA =========================

engineUsed = 1;

% Creation of radius vector --> for now we are considering same radius for
% interstage and stage
mission.structure.radius = [mission.capsule.radius];

for ii = launcher(1):-1:1
    mission.structure.radius = [mission.structure.radius configuration.geometry.stage{ii}.radius configuration.geometry.stage{ii}.radius];
end

% Pressure vector creation
nComponents = length(mission.structure.componentLength);

mission.structure.pressurization = zeros(nComponents, 1);

for ii = 3:2:nComponents
    mission.structure.pressurization(ii) = mission.structure.tankPressure;
end

mission    = thicknessFunction(mission, engineUsed);

thick_mm = mission.structure.thickness * 1e3;   % [mm]
fprintf('Thicknesses of the structures starting from the nose are:\n');
fprintf('  %.3f mm\n', thick_mm);  % stampa un valore per riga
mStruct_ton = mission.structure.mStruct * 1e-3; % [ton]
fprintf('Masses of the structures starting from the nose are:\n');
fprintf('  %.3f tons\n', mStruct_ton);

%%
xcg = computeXCG(mission, opt);
equatMoment = @(alphafins) mission.structure.tMaxQ(2) * (mission.structure.launcherLength - ...
    xcg) + (mission.structure.dMaxQ(2)+mission.structure.lMaxQ(2)) * (xcg - ...
    xcp_a) + (norm(mission.structure.liftFinsMaxQ) * cos(alphafins) - ...
    norm(mission.structure.dragFinsMaxQ) * sin(alphafins)) * (xcp_a - xcg); 
sol = fzero(equatMoment, mission.structure.alphaQmax)

alphas = deg2rad(linspace(-180, 180, 400)); % in gradi per esempio
plot(rad2deg(alphas), equatMoment(alphas)), grid on
xlabel('\alpha_{fins} [deg]')
ylabel('M(\alpha)')


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