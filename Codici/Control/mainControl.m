clc
close all
clear all

% Optimal Solution form an old GA
launcher = [2,1,4,4,0.4056,0.4016,0.7];

[mission,settings] = dataStructGlobal;

for i = 1:launcher(1)
    configuration.stage{i}.engine = mission.engines{launcher(1+i)};
end

[mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);

[xGATraj, fvalGATraj] = ga( @(x) objFunGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings), ...
                        launcher(1)*2*settings.nOptPointsTraj,...
                        [],[],[],[],settings.lowerBoundsGA,settings.upperBoundsGA, ...
                        @(x) nlconGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings),settings.gaTrajOptions);


thrustData = reshape(xGATraj,settings.nOptPointsTraj,2,2);

% Upload state and time from Traj2D
load('stateCollocation.mat') ;
load('timeCollocation.mat') ;

% Manca il reshape di guidancePoints

% Manca Rotazione sul piano 2 degli stateCollocation

% Initial Condition
y0 = stateCollocation(:,1,1) ;

% Upload mission data
tSpan = [0 500];

Nstage = launcher(1) ;

%% GA for tuning gains

% lowerBounds = [10 1 10 10 10 10 100 100]'; 
% upperBounds = [1000 1 1000 1000 1000 1000 1000 1000]'; 
lowerBounds = 1 * ones(8,1) ; 
upperBounds = 1000 * ones(8,1) ;
nVars = length(lowerBounds) ;

gains0 = [50 0 10 1 0 60 100 1000]' ;
intCon = 1:8 ;

options_ga = optimoptions("ga", ...
    "Display","iter", ...
    "MaxGenerations",30, ...
    "PopulationSize",200,...
    "UseParallel",true,...
    "FunctionTolerance", 1e-4,...
    "MaxStallGenerations", 10,...
    'EliteCount',  6);


[gains] = ga(@(x) findGains(x,guidancePoints,y0,tSpan), nVars, [],[],[],[],lowerBounds,upperBounds, [] ,intCon, options_ga) ;


options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t, sol] = ode113(@(t,x) launcherDynamicsAndControlECI(t, x,thrustData, mission, configuration, launcher, stageNumber, guidancePoints, gains), tSpan, y0,options);




function [output] = findGains(x,guidancePoints,y0,tSpan)

gains = x ;

if ~iscolumn(gains)
    gains = gains' ;
end


options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t, sol] = ode113(@(t,x) launcherDynamicsAndControlECI(t, x,thrustData, mission, configuration, launcher, stageNumber, guidancePoints, gains), tSpan, y0,options);


output = sum(sol(:, 1:2)' - guidancePoints(1:2,:)) + sum(sol(:, 3:4)' - guidancePoints(3:4, :)) ;


end