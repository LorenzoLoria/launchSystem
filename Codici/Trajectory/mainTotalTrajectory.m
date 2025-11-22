clear all
clc
close all

%% Path Directory

addpath(genpath("..\..\"))

%% Upload Mission Struct
tSpan = [0 600];

[mission,opt] = dataStruct;
thrustDataVec = [[1; 0.6 ; 0.4; 0.1; 0.05] , [0; 30; 30; 30; 30] ];
tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));
%thrustData = @(t) [interp1(tVec,thrustDataVec(:,1),t);interp1(tVec,thrustDataVec(:,2),t);interp1(tVec,thrustDataVec(:,3),t)];
F_thrust = griddedInterpolant(tVec, thrustDataVec, 'linear', 'none');  
thrustData = @(t) F_thrust(t).';  

[timeCollocation, stateCollocation] = totalTrajectory(mission,opt,thrustData);
%%

figure(1)

EarthPlot(mission.environment.rEarth)
hold on
plot3(stateCollocation(1,:,1),stateCollocation(2,:,1),stateCollocation(3,:,1))
plot3(stateCollocation(1,:,2),stateCollocation(2,:,2),stateCollocation(3,:,2))
%plot3(stateCollocation(1,:,3),stateCollocation(2,:,3),stateCollocation(3,:,3))

hold off


%%
stateCollocation(4:6,end,2)
