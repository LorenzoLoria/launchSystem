
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

%% Nominal Trajectory

[timeCollocation, stateCollocation] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustData);
figure
EarthPlot(mission.environment.rEarth)
hold on
plot3(stateCollocation(1,:,1),stateCollocation(2,:,1),stateCollocation(3,:,1))
plot3(stateCollocation(1,:,2),stateCollocation(2,:,2),stateCollocation(3,:,2))
plot3(stateCollocation(1,:,3),stateCollocation(2,:,3),stateCollocation(3,:,3))
plot3(mission.target.initialPointECI(1),mission.target.initialPointECI(2),mission.target.initialPointECI(3),'ro')
%%
clc
% Initialization of the variables
sizeMC = 10;
distanceFromTarget = zeros(1,sizeMC);
cumulativeMean = zeros(length(distanceFromTarget),1);
hVec = 0:100;

% Computation of the Wind Profiles

[meanWind, varWind] = GRAM07_HWM07_annual(hVec);
windUncertainty = sqrt(varWind) .* randn(sizeMC,1);
WindVelocityMag = meanWind + windUncertainty;
windAngVel = WindVelocityMag ./ (mission.environment.rEarth + hVec);
lonInitial = mission.launchBase.lonInitial;
montecarlo.vxWind = - windAngVel .* (mission.environment.rEarth + hVec) .* sin(lonInitial) ;
montecarlo.vyWind = windAngVel .* (mission.environment.rEarth + hVec) .* cos(lonInitial) ;
% Functions for wind profile on ECI (rotated inside the dynamics)
windVelXFun = griddedInterpolant(hVec,montecarlo.vxWind(1,:),'linear','linear');
windVelYFun = griddedInterpolant(hVec,montecarlo.vyWind(1,:),'linear','linear');

[t, x] = totalTrajectoryControlled(launcher,configuration,mission,settings,windVelXFun,windVelYFun,stateCollocation,timeCollocation,2000000,10,thrustData)