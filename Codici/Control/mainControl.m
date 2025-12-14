clc
close all
clear all

% Optimal Solution form an old GA
% launcher = [2,1,4,4,0.4056,0.4016,0.7];
launcher = [2,2,3,3,0.457979725,0.876723067, 0.5] ;

% Path Directory

addpath(genpath("..\..\"))

% Upload Mission Struct
[mission,settings] = dataStructGlobal;

for i = 1:launcher(1)
    configuration.stage{i}.engine = mission.engines{launcher(1+i)};
end

[mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);

stageNumber = 1;
m = 150e3;
[inertia,Xcg] = InertiaEvaluation(mission, configuration, mer, stageNumber, m, launcher);

%%
[xGATraj, fvalGATraj] = ga( @(x) objFunGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings), ...
                        launcher(1)*2*settings.nOptPointsTraj,...
                        [],[],[],[],settings.lowerBoundsGA,settings.upperBoundsGA, ...
                        @(x) nlconGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings),settings.gaTrajOptions);


thrustData = reshape(xGATraj,settings.nOptPointsTraj,2,2);
%%
[guidanceTime, guidancePoints] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustData);

figure
EarthPlot(mission.environment.rEarth)
hold on
plot3(guidancePoints(1,:,1),guidancePoints(2,:,1),guidancePoints(3,:,1))
plot3(guidancePoints(1,:,2),guidancePoints(2,:,2),guidancePoints(3,:,2))
plot3(guidancePoints(1,:,3),guidancePoints(2,:,3),guidancePoints(3,:,3))
plot3(mission.target.initialPointECI(1),mission.target.initialPointECI(2),mission.target.initialPointECI(3),'ro')

guidancePoints2D(:,:,1) = [mission.target.Rfinal*guidancePoints(1:3,:,1);mission.target.Rfinal*guidancePoints(4:6,:,1)];
guidancePoints2D(:,:,2) = [mission.target.Rfinal*guidancePoints(1:3,:,2);mission.target.Rfinal*guidancePoints(4:6,:,2)];
guidancePoints2D(:,:,3) = [mission.target.Rfinal*guidancePoints(1:3,:,3);mission.target.Rfinal*guidancePoints(4:6,:,3)];

%%
gains = ones(1,16);
gains(7)=100000;
[timeCollocationControlled, stateCollocationControlled] = totalTrajectoryControl(mission,mer,configuration,settings, launcher,thrustData, guidancePoints, guidanceTime, gains) ;
figure
EarthPlot(mission.environment.rEarth)
hold on
plot3(stateCollocationControlled(1,:,1),stateCollocationControlled(2,:,1),stateCollocationControlled(3,:,1))
plot3(stateCollocationControlled(1,:,2),stateCollocationControlled(2,:,2),stateCollocationControlled(3,:,2))
plot3(stateCollocationControlled(1,:,3),stateCollocationControlled(2,:,3),stateCollocationControlled(3,:,3))
plot3(mission.target.initialPointECI(1),mission.target.initialPointECI(2),mission.target.initialPointECI(3),'ro')

%% GA for tuning gains

% lowerBounds = [10 1 10 10 10 10 100 100]'; 
% upperBounds = [1000 1 1000 1000 1000 1000 1000 1000]'; 
lowerBounds = 0 * ones(16,1) ; 
upperBounds = inf * ones(16,1) ;
nVars = length(lowerBounds) ;

% gains0 = [50 0 10 1 0 60 100 10 4 2 8]' ;
% intCon = 1:12 ;

options_ga = optimoptions("ga", ...
    "Display","iter", ...
    "MaxGenerations",30, ...
    "PopulationSize",300,...
    "UseParallel",true,...
    "FunctionTolerance", 1e-4,...
    "MaxStallGenerations", 10,...
    'EliteCount',  10);


[gains] = ga(@(x) findGains(x, mission, mer, configuration, settings, launcher, guidancePoints, guidanceTime, thrustData), nVars, [],[],[],[],lowerBounds,upperBounds, [] ,[], options_ga) ;


%%
figure
% gains = 14.3472    0.6593   17.7579    1.3554    2.5457   17.2149    8.6441    4.2559  5.0148    0.3572   14.3324    8.4757   10.5220   20.7360   12.9394   10.9536
%gains = [16.5052    0.0999   16.4921    8.8607   11.2140   14.2543   26.0408   11.2452 12.0506    3.1178    7.7157   19.7268   15.6168   15.9040    3.5497    4.0199]
% gains = [6.5460    0.7318   11.9055   18.8923    9.3732   11.3010   19.4130    4.2449    1.3790    1.8875    6.7385   11.3978 10.9511   17.3521   18.0000   12.0162]
[timeCollocationControlled, stateCollocationControlled] = totalTrajectoryControl(mission,mer,configuration,settings, launcher,thrustData, guidancePoints, guidanceTime, gains);
figure
%EarthPlot(mission.environment.rEarth)
hold on

plot3(stateCollocationControlled(1,:,1),stateCollocationControlled(2,:,1),stateCollocationControlled(3,:,1),'g--','LineWidth',2)
plot3(stateCollocationControlled(1,:,2),stateCollocationControlled(2,:,2),stateCollocationControlled(3,:,2),'g--','LineWidth',2)
plot3(stateCollocationControlled(1,:,3),stateCollocationControlled(2,:,3),stateCollocationControlled(3,:,3))
%plot3(mission.target.initialPointECI(1),mission.target.initialPointECI(2),mission.target.initialPointECI(3),'ro')
hold on

plot3(guidancePoints(1,:,1),guidancePoints(2,:,1),guidancePoints(3,:,1),'r')
plot3(guidancePoints(1,:,2),guidancePoints(2,:,2),guidancePoints(3,:,2),'r')
plot3(guidancePoints(1,:,3),guidancePoints(2,:,3),guidancePoints(3,:,3))
axis equal









function [objective] = findGains(x,mission, mer, configuration, settings,launcher, guidancePoints,guidanceTime, thrustDataVec)

 gains = x ;

global angleErr
 [timeControlled, stateCollocationControlled] = totalTrajectoryControl(mission,mer,configuration, settings, launcher,thrustDataVec, guidancePoints, guidanceTime, gains);

 objective = sum(abs(stateCollocationControlled(1:6, end, end-1)  -  guidancePoints(1:6,end,end-1))+abs(stateCollocationControlled(1:6, end, end-2)  -  guidancePoints(1:6,end,end-2)))/1000 + abs(angleErr)*10 ;
 % error=zeros(1,launcher(1));
 % for i = 1:launcher(1)
 % 
 %     guidancePointsFun = @(t) interp1(guidanceTime(:,i),guidancePoints(1:6,:,i),'linear','extrap');
 %     error(i) = sum(abs(stateCollocationControlled(1:6,:,i)-guidancePointsFun(timeControlled)),"all");
 % end
%objective = sum(error);
 %error = vecnorm(stateCollocationControlled(1:3, :, :) - guidancePoints(1:3, :, :));
 %objective = sum(error, 'all') + 10 * finalError; % Weighted sum

end