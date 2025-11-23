clear all
clc
close all

%% Path Directory

addpath(genpath("..\..\"))

%% Upload Mission Struct
tSpan = [0 500];

[mission,opt] = dataStruct;
%%

% Generate optimisation variables for general stages

thrustDataVec1 = [[1; 1 ;1; 1; 1] , [0; 0; 0; 0; 0] ];
thrustDataVec2 = [[1; 1 ; 1; 1; 1] , [0; 20; 30; 50; 70] ];

thrustData(:,:,1) =thrustDataVec1;
thrustData(:,:,2) =thrustDataVec2;

%tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec1,1));
%F_thrust1 = griddedInterpolant(tVec, thrustDataVec1, 'linear', 'none');  
%thrustDataFun1 = @(t) F_thrust1(t);  

%tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec1,1));
%F_thrust2 = griddedInterpolant(tVec, thrustDataVec2, 'linear', 'none');  
%thrustDataFun2 = @(t) F_thrust2(t);  


%thrustData{1} = thrustDataFun1;
%thrustData{2} = thrustDataFun2;
[timeCollocation, stateCollocation] = totalTrajectory(mission,opt,thrustData);



%%
figure(1)

EarthPlot(mission.environment.rEarth)
hold on
plot3(stateCollocation(1,:,1),stateCollocation(2,:,1),stateCollocation(3,:,1))
plot3(stateCollocation(1,:,2),stateCollocation(2,:,2),stateCollocation(3,:,2))
plot3(stateCollocation(1,:,3),stateCollocation(2,:,3),stateCollocation(3,:,3))
hold off


%%
figure
plot(timeCollocation(:),[stateCollocation(7,:,1),stateCollocation(7,:,2) ,stateCollocation(7,:,3)]/1000)

figure
plot(diff([stateCollocation(7,:,1),stateCollocation(7,:,2) ,stateCollocation(7,:,3)]))
%%
%stateCollocation(4:6,end,2)
