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

Nstage = launcher(1) ;

% Upload state and time from Traj2D
load('stateCollocation.mat') ;
load('timeCollocation.mat') ;

guidancePoints = stateCollocation ;
guidanceTime = timeCollocation; 


%% GA for tuning gains

% lowerBounds = [10 1 10 10 10 10 100 100]'; 
% upperBounds = [1000 1 1000 1000 1000 1000 1000 1000]'; 
lowerBounds = 1 * ones(12,1) ; 
upperBounds = 1000 * ones(12,1) ;
nVars = length(lowerBounds) ;

gains0 = [50 0 10 1 0 60 100 10 4 2 8]' ;
intCon = 1:12 ;

options_ga = optimoptions("ga", ...
    "Display","iter", ...
    "MaxGenerations",30, ...
    "PopulationSize",200,...
    "UseParallel",false,...
    "FunctionTolerance", 1e-4,...
    "MaxStallGenerations", 10,...
    'EliteCount',  6);


[gains] = ga(@(x) findGains(x,mission, mer, configuration, settings, launcher, guidancePoints,guidanceTime, thrustData), nVars, [],[],[],[],lowerBounds,upperBounds, [] ,intCon, options_ga) ;


options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t, sol] = ode113(@(t,x) launcherDynamicsAndControlECI(t, x,thrustData, mission, configuration, launcher, stageNumber, guidancePoints, guidanceTime, gains), tSpan, y0,options);









function [objective] = findGains(x,mission, mer, configuration, settings,launcher, guidancePoints,guidanceTime, thrustDataVec)

gains = x ;

if ~iscolumn(gains)
    gains = gains' ;
end

[~, stateCollocationControlled] = totalTrajectoryControl(mission,mer,configuration, settings, launcher,thrustDataVec, guidancePoints, guidanceTime, gains);

objective = 0 ;

for ii = 1 : launcher(1)
    output = sum(abs(stateCollocationControlled(1:3, :, ii) - guidancePoints(1:3,:))) +...
             sum(abs(stateCollocationControlled(4:6, :, ii) - guidancePoints(4:6,:))) ;
    objective = objective + output ;
end

end