
% Initialization

clear all
clc
close all

addpath(genpath("..\..\"))

[mission,settings] = dataStructGlobal;

launcher = [2,1,4,4,0.4056,0.4016,0.7];

for i = 1:launcher(1)
    configuration.stage{i}.engine = mission.engines{launcher(1+i)};
end

[mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);

% Optimal nominal trajectory

[xGATraj, fvalGATraj] = ga( @(x) objFunGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings), ...
                        launcher(1)*2*settings.nOptPointsTraj,...
                        [],[],[],[],settings.lowerBoundsGA,settings.upperBoundsGA, ...
                        @(x) nlconGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings),settings.gaTrajOptions);

thrustData(:,:,1) = [xGATraj(1:5)',xGATraj(6:10)'];
thrustData(:,:,2) = [xGATraj(11:15)',xGATraj(16:20)'];

% Nominal Trajectory

[timeCollocation, stateCollocation] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustData);
%%

[timeCollocation, stateCollocation] = totalTrajectoryControlled(launcher,opt,mission,settings,windVelXFun,windVelYFun,stateCollocationRef,timeCollocationRef)