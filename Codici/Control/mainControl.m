clear all
clc
close all

addpath(genpath("..\..\"))

[mission,settings] = dataStructGlobal;

launcher = [2	2	3	3	0.459952176990556	0.753370531158904	0.634795741885559];

for i = 1:launcher(1)
    configuration.stage{i}.engine = mission.engines{launcher(1+i)};
end

[mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);

 thrustDataVecFMC = [0.901784296794966
0.999967458028643
0.900002561600819
0.900000000000007
0.900000000000006
1.32376427751106
22.4721810110826
51.9025769669441
58.3043070770380
54.0837349292332
0.400884818399135
0.966674994308654
0.974049905676140
0.992054604106613
0.992888695383693
64.8416016397151
78.7571907487398
90.6394069314900
87.6489451551192
98.5947513930343];

thrustData = reshape(thrustDataVecFMC,settings.nOptPointsTraj,2,2);

    localLowerBoundsFMC = settings.lowerBoundsFMC(:,:,1:launcher(1)) ;
    localUpperBoundsFMC = settings.upperBoundsFMC(:,:,1:launcher(1)) ;
[thrustDataVecFMC,fvalFMCTraj,~,checkViolation] = fmincon ( @(x)objFunFMCTraj(x,launcher,configuration,mission,settings),...
            thrustData,[],[],[],[],...
            localLowerBoundsFMC-eps,localUpperBoundsFMC+eps,...
            @(x) nlconFMCTraj(x,launcher,configuration,mission,settings),...
            settings.fminconTrajOptions);

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
% intCon = 1:12 ;9

options_ga = optimoptions("ga", ...
    "Display","iter", ...
    "MaxGenerations",30, ...
    "PopulationSize",200,...
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
% gains = [8.4549    0.4456    6.1723   19.2523    5.1122   18.9453   20.5953    5.9918    4.0416    0.6793    1.3930   11.8240  16.6428   13.3837   17.3410   13.1817]
% gains = [10.9705    0.5071    4.4102   11.9522    4.8867    8.1187   20.5295   16.4141    2.9761    5.9335    5.6457   19.2647   20.1170   18.0131   14.9054  4.8741]
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