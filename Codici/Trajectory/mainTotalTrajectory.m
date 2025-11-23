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
thrustDataVec2 = [[1; 1 ; 1; 1; 1] , [0; 60; 70; 80; 90] ];

thrustData(:,:,1) =thrustDataVec1;
thrustData(:,:,2) =thrustDataVec2;

[timeCollocation, stateCollocation] = totalTrajectory(mission,opt,thrustData);

%%




% Initial Guess using GA
obj_ga = @(x) objFunGA( reshape(x,mission.optimisation.GA.variables,2), mission);
nonlcon_ga = @(x) nlconGA( reshape(x,mission.optimisation.GA.variables,2), mission );

lbGA = [0.1*ones(mission.optimisation.GA.variables,1);0*ones(mission.optimisation.GA.variables,1)];
ubGA = [ones(mission.optimisation.GA.variables,1);90*ones(mission.optimisation.GA.variables,1)];






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
