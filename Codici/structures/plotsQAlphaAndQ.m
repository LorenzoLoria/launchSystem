%% ======================= Q and Q alpha ==================================

clear all;
clc;
close all

% ============================= DATI ======================================

addpath(genpath('..\..\'))

[mission, settings] = dataStructGlobal;

mission.structure.alphaQmax = 5 * pi / 180;

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

vel = stateCollocation(4:6,:,1:end-1)-stateCollocation(4:6,1,1);
vel = mission.target.Rfinal* vel(1:3,:);
absVel = sqrt ( vel(1,:).^2+ vel(2,:).^2 + vel(3,:).^2 );

stageNumber = 1;


pos = stateCollocation(1:3,:,1:end-1);
pos = pos(1:3,:);
absH = sqrt ( pos(1,:).^2+ pos(2,:).^2 + pos(3,:).^2 )-mission.environment.rEarth;

vxWind = mission.environment.windXFun(absH/1000);
vyWind = mission.environment.windYFun(absH/1000);
 
wind = [vxWind;vyWind;zeros(1,length(vxWind))];

totalVel = vel - wind;

for i = 1:length(totalVel)
    alpha(i) = acos (dot(totalVel(:,i) , vel(:,i)) / norm(vel(:,i))/norm(totalVel(:,i) ));
end

mass = stateCollocation(7, :, 1:end-1);
mass = mass(:);

rhoVec = mission.environment.rhoFun(absH);

q = 0.5*rhoVec.*(absVel).^2;

timeStage = timeCollocation(:,1:end-1);
timeStage = timeStage(:);

%% =========================== PLOTS =======================================
q = q * 1e-3; % kPa
qAlpha = q .* alpha * 180 / pi; % kPa * deg

qPlots = figure(1);

yyaxis left
plot(timeStage, q, 'Color', settings.color.blu, 'LineWidth', 1.5)
ylabel('$q$ [kPa]', 'Interpreter', 'latex', 'FontSize', 12, 'Color', settings.color.blu)
yyaxis right
hold on
plot(timeStage, qAlpha, 'Color', settings.color.orange, 'LineWidth', 1.5)
ylabel('$q\alpha$ [kPa $\cdot$ deg]', 'Interpreter', 'latex', 'FontSize', 12, 'Color', settings.color.orange)
grid on
xlabel('t [s]')
[qmax,idx1] = max(q);
[qAlphaMax, idx2] = max(qAlpha);

% xline(timeStage(idx1), 'Color', settings.color.green)
% xline(timeStage(idx2), 'Color', settings.color.green)

% yline(qmax, 'Color', settings.color.green)
% yline(qAlphaMax, 'Color', settings.color.green)

xlim([0, timeStage(115)])
legend('q [kPa]', 'q$\alpha$ [kPa $\cdot$ deg]')

ax = gca;
ax.YAxis(1).Color = settings.color.blu;       
ax.YAxis(2).Color = settings.color.orange; 

setPlotSettings(title(''))
exportStandardizedFigure(qPlots,'qPlots',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','..\..\figures\structures')
ax = gca;
ax.YAxis(1).Color = settings.color.blu;       
ax.YAxis(2).Color = settings.color.orange; 