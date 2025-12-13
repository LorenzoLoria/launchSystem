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


% 2	2	3	3	0.598529436065056	0.760844390761651	0.751703499394933
% 2	2	3	3	0.582988078347329	0.771358165933520	0.558998974464695
% 2	4	3	3	0.556470803677211	0.923400327579640	0.790985847483056
% 2	2	3	3	0.581204806211464	0.912935133486741	0.576488534825573
% 2	2	3	3	0.592957115494177	0.812990434857997	0.473501153501329
% 2	2	3	3	0.457979725217222	0.876723067007092	0.447092534640540
%2	2	3	3	0.457979725217222	0.876723067007092	0.447092534640540
% 2	4	3	3	0.546875000000000	0.931250000000000	0.481250000000000



%% Skip GA

% launcher = [2	2	3	3	0.457979725217222	0.876723067007092	0.447092534640540];
% 
% for i = 1:launcher(1)
%     configuration.stage{i}.engine = mission.engines{launcher(1+i)};
% end
% 
% 
% [outputOBJ,outputNLC] = launcherSimulation(launcher,mission,settings,1);
% 
% 
% 
% %% Plots
% launcher = [2	2	3	3	0.457979725217222	0.876723067007092	0.447092534640540];
% for i = 1:launcher(1)
%     configuration.stage{i}.engine = mission.engines{launcher(1+i)};
% end
% 
% [mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);
% 
% thrustDataVecFMC = [0.999840298811896
% 0.997963997732943
% 0.920675105294849
% 0.996458841895991
% 0.994287377368226
% 2.45190950938795
% 19.9182379061273
% 41.8057234345297
% 59.5264606504011
% 52.3998780994523
% 0.400024967020885
% 0.999425045100713
% 0.997237143840542
% 0.992542430970901
% 0.999465476966899
% 51.4739605677812
% 76.4824021625324
% 83.6432540377504
% 92.6291287921728
% 93.2465223861152];
% 
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
