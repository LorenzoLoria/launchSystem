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

% launcher = [2,2,3,3,0.399,0.788,0.478];
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
% thrustData(:,:,1) = [xGATraj(1:5)',xGATraj(6:10)'];
% thrustData(:,:,2) = [xGATraj(11:15)',xGATraj(16:20)'];
% 
% 
% %% Plots
% 
% thrustData = reshape(xGATraj,settings.nOptPointsTraj,2,2);
% 
% [timeCollocation, stateCollocation] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustData);
% figure(1)
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
% 
% 
% 
% 
