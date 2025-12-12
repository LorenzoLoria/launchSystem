clear all
clc
close all

%% Path Directory

addpath(genpath("..\..\"))

%% Upload Mission Struct
[mission,settings] = dataStructGlobal;

%% GA


[x_ga, fval_ga] = ga(@(x) objFunGlobalGA(x,mission,settings),...
                    settings.globalGAoptVariables,[],[],[],[],...
                    settings.lowerBoundsGlobalGA,settings.upperBoundsGlobalGA,...
                    @(x) nlconGlobalGA(x,mission,settings),...
                    settings.intconGlobalGA,settings.globalGAOptions);

%% Skip GA

% launcher = [2	1	3	3	0.591385004532601	0.854996724580398	0.601814212277586];
% for i = 1:launcher(1)
%     configuration.stage{i}.engine = mission.engines{launcher(1+i)};
% end
% 
% [mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);
% 
% [xGATraj, fvalGATraj] = ga( @(x) objFunGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings), ...
%                         launcher(1)*2*settings.nOptPointsTraj,...
%                         [],[],[],[],settings.lowerBoundsGA,settings.upperBoundsGA, ...
%                         @(x) nlconGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings),settings.gaTrajOptions);
% 
% xGATrajRS = reshape(xGATraj,settings.nOptPointsTraj,2,launcher(1));
% %%
% localLowerBoundsFMC = settings.lowerBoundsFMC(:,:,1:launcher(1)) ;
%     localUpperBoundsFMC = settings.upperBoundsFMC(:,:,1:launcher(1)) ;
% 
% 
% [thrustDataVecFMC,fvalFMCTraj,~,maxContrViol] = fmincon ( @(x)objFunFMCTraj(x,launcher,configuration,mission,settings),...
%             xGATrajRS,[],[],[],[],...
%             localLowerBoundsFMC-eps,localUpperBoundsFMC+eps,...
%             @(x) nlconFMCTraj(x,launcher,configuration,mission,settings),...
%             settings.fminconTrajOptions);
% 
% 
% thrustDataVecFMC = maxContrViol.bestfeasible.x;
% thrustData(:,:,1) = [thrustDataVecFMC(1:5)',thrustDataVecFMC(6:10)'];
% thrustData(:,:,2) = [thrustDataVecFMC(11:15)',thrustDataVecFMC(16:20)'];
% 
% 
% %% Plots
% 
% thrustData = reshape(thrustDataVecFMC,settings.nOptPointsTraj,2,2);
% 
% [timeCollocation, stateCollocation] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustData);
% figure
% EarthPlot(mission.environment.rEarth)
% hold on
% plot3(stateCollocation(1,:,1),stateCollocation(2,:,1),stateCollocation(3,:,1),'r')
% plot3(stateCollocation(1,:,2),stateCollocation(2,:,2),stateCollocation(3,:,2),'y')
% plot3(stateCollocation(1,:,3),stateCollocation(2,:,3),stateCollocation(3,:,3),'g')
% plot3(mission.target.initialPointECI(1),mission.target.initialPointECI(2),mission.target.initialPointECI(3),'bo')
% targetFinalLat = mission.target.latInitial ; 
% targetFinalLon = mission.target.lonInitial + mission.environment.omega * timeCollocation(end,end); 
% targetFinalPosECI = 6371000*[cos(targetFinalLat)*cos(targetFinalLon); cos(targetFinalLat)*sin(targetFinalLon); sin(targetFinalLat) ];
% plot3(targetFinalPosECI(1),targetFinalPosECI(2),targetFinalPosECI(3), 'ob')
% title("Trajectory")
% hold off
% 
% 
% %%
% 
% 
% 
% figure
% plot([thrustData(:,2,1);thrustData(:,2,2)])
% title("Angle1")
% 
% figure
% plot([thrustData(:,1,1);thrustData(:,1,2)])
% title("Angle1")
