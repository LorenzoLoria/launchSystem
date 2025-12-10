%% ======================= Prova Runnata ================================== 

clear all;
clc
close all

% ==========================  DATI ========================================

addpath(genpath('..\..\'))

[mission, opt] = dataStruct;

X = load('..\Trajectory\Final Trajectory\InitialPop1stages.mat','X');
X = X.X(:,end);
thrustData(:,:,1) = [X(1:5),X(6:10)];
%thrustData(:,:,2) = [X(11:15),X(16:20)];
%thrustData(:,:,3) = [X(21:25),X(26:30)];

[timeCollocation, stateCollocation] = totalTrajectory(mission,opt,thrustData,1);
EarthPlot(mission.environment.rEarth)
hold on
plot3(stateCollocation(1,:,1),stateCollocation(2,:,1),stateCollocation(3,:,1),'r')
plot3(stateCollocation(1,:,2),stateCollocation(2,:,2),stateCollocation(3,:,2),'g')
% plot3(stateCollocation(1,:,3),stateCollocation(2,:,3),stateCollocation(3,:,3),'b')
% plot3(stateCollocation(1,:,4),stateCollocation(2,:,4),stateCollocation(3,:,4),'g')
plot3(mission.target.initialPointECI(1),mission.target.initialPointECI(2),mission.target.initialPointECI(3),'bo')
targetFinalLat = mission.target.latInitial ; 
targetFinalLon = mission.target.lonInitial + mission.target.omega * timeCollocation(end,end); 
targetFinalPosECI = 6371000*[cos(targetFinalLat)*cos(targetFinalLon); cos(targetFinalLat)*sin(targetFinalLon); sin(targetFinalLat) ];
plot3(targetFinalPosECI(1),targetFinalPosECI(2),targetFinalPosECI(3), 'ob')
title("Trajectory")
hold off


%norm(stateCollocation(1:3,end,end) - mission.target.initialPointECI)
%%
mission.structure.alphaQmax = deg2rad(0);

if opt.nStages == 1
    mission.structure.nElements = 8;
elseif opt.nStages == 2
    mission.structure.nElements = 12;
elseif opt.nStages == 3
    mission.structure.nElements = 16;
end

mission.structure.nNodes = mission.structure.nElements + 1;

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

[mission] = externalLoads(timeCollocation,stateCollocation,mission,opt,thrustData, mission.structure.alphaQmax);

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
mission.structure.nComponents = length(mission.structure.componentLength);

% ====================== CALCOLO AZIONI INTERNE ===========================
mission.structure.Ftfins = (- mission.structure.tMaxQ(2) * (mission.structure.launcherLength - xcg) + ( mission.structure.dMaxQ(2) + mission.structure.lMaxQ(2)) * (xcg - xcp)) / (mission.structure.launcherLength - xcp_a - xcg);
[mission] = loadsFinder(mission);

N = mission.structure.N;
T = mission.structure.T;
M = mission.structure.M;


% ==================== SPESSORE E MASSA STRUTTURA =========================
engineUsed = 1;

% Pressure vector creation
mission.structure.pressurization = zeros(mission.structure.nElements, 1);

for ii = 3:2:mission.structure.nComponents
    mission.structure.pressurization(ii) = mission.structure.tankPressure;
end

mission = thicknessFunction(mission, engineUsed);

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

for i = 1:mission.structure.nElements
    
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

Mach = mission.structure.machNumber;
be  = mission.aerodynamics.finsGeom.be;
Se  = mission.aerodynamics.finsGeom.Se;
Sref = mission.aerodynamics.bodyGeom.Aref;
cmac = mission.aerodynamics.finsGeom.cmac;
delta_le = mission.aerodynamics.finsGeom.delta_le;
lambda_le = mission.aerodynamics.finsGeom.lambda_le;
b = mission.aerodynamics.finsGeom.b;
tmac = mission.aerodynamics.finsGeom.tmac;
Nfins = mission.aerodynamics.finsGeom.Nfins;
q = mission.structure.dynamicPressure;
alphaFin = deg2rad(10);
A = be^2/Se;
M_ale = Mach * cosd(lambda_le);


% --- Normal force coefficient
if Mach > sqrt(1 + (8/(pi*A))^2)
    CN_surf = ((4*abs(sin(alphaFin)*cos(alphaFin)) / sqrt(Mach^2 - 1)) + 2*sin(alphaFin)^2) * Se / Sref;
else
    CN_surf = ((pi*A/2*abs(sin(alphaFin)*cos(alphaFin)) + 2*sin(alphaFin)^2) * Se / Sref);
end

% --- CD0 surface friction
CD0_surf_friction = 0.0133 * (Mach / (q*cmac))^0.2 * 2 * Se / Sref;

% --- CD0 surface wave

if M_ale < 1
    CD0_surf_wave = 0;
else
    CD0_surf_wave = (1.429 / M_ale^2) * ((1.2*M_ale^2)^3.5 * (2.4/(2.8*M_ale^2 - 0.4))^2.5 - 1) * (sin(deg2rad(delta_le))^2 * cos(deg2rad(lambda_le)) * tmac * b) / Sref;
end

CN_fins_tot = Nfins * CN_surf;
CD0_fins_tot = Nfins * (CD0_surf_friction + CD0_surf_wave);
CA_fins_tot = CD0_fins_tot * cos(alphaFin)^2;
% CL, CD Fins
finsCL = CN_fins_tot * cos(alphaFin) - CA_fins_tot * sin(alphaFin);
finsCD = CA_fins_tot * cos(alphaFin) + CN_fins_tot * sin(alphaFin);

Lfin = mission.structure.dynamicPressure * mission.aerodynamics.bodyGeom.Aref * finsCL * 2;
Dfin = mission.structure.dynamicPressure * mission.aerodynamics.bodyGeom.Aref * finsCD *2;

LfinBody = [Lfin * sin(alphaFin) ; Lfin * cos(alphaFin)];
DfinBody = [Dfin * cos(alphaFin) ; -Dfin * sin(alphaFin)];

% thrust angle such that the torque is balanced
mission.structure.delta = asind((-(LfinBody(2) + DfinBody(2)) * (mission.structure.launcherLength - xcp_a - xcg) + (mission.structure.dMaxQ(2) + mission.structure.lMaxQ(2)) * (xcg - xcp))/(norm(mission.structure.tMaxQ)*(mission.structure.launcherLength - xcg)));


%%
momentFins = Ftfins * xcp_a
CLnew = abs(Ftfins) / (mission.structure.dynamicPressure * mission.aerodynamics.bodyGeom.Aref)
