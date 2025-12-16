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






% 2	2	3	3	0.459952176990556	0.753370531158904	0.634795741885559


%% Skip GA

launcher = [3	2	3	3	0.532694822285913	0.339319975934292	0.727608416114238];

for i = 1:launcher(1)
    configuration.stage{i}.engine = mission.engines{launcher(1+i)};
end


[outputOBJ,outputNLC] = launcherSimulation(launcher,mission,settings,1);



%% Plots
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

%% PLOT
%Trajectory Plot
[timeCollocation, stateCollocation] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustDataVecFMC);

trajectory1 = figure(1);
setPlotSettings(title(""))
EarthPlot(mission.environment.rEarth)
hold on
plot3(stateCollocation(1,:,1),stateCollocation(2,:,1),stateCollocation(3,:,1),'Color',settings.color.terracotta,'LineWidth',1.5)
plot3(stateCollocation(1,:,2),stateCollocation(2,:,2),stateCollocation(3,:,2),'Color',settings.color.orange,'LineWidth',1.5)
plot3(stateCollocation(1,:,3),stateCollocation(2,:,3),stateCollocation(3,:,3),'Color',settings.color.gray,'LineWidth',1.5)
plot3(mission.target.initialPointECI(1),mission.target.initialPointECI(2),mission.target.initialPointECI(3))
view(240,10)
legend(["" "First Stage" "Second Stage" "Capsule" ""])
hold off


trajectory2 = figure(2);
setPlotSettings(title(""))
EarthPlot(mission.environment.rEarth)
hold on
plot3(stateCollocation(1,:,1),stateCollocation(2,:,1),stateCollocation(3,:,1),'Color',settings.color.terracotta,'LineWidth',1.5)
plot3(stateCollocation(1,:,2),stateCollocation(2,:,2),stateCollocation(3,:,2),'Color',settings.color.orange,'LineWidth',1.5)
plot3(stateCollocation(1,:,3),stateCollocation(2,:,3),stateCollocation(3,:,3),'Color',settings.color.gray,'LineWidth',1.5)
plot3(mission.target.initialPointECI(1),mission.target.initialPointECI(2),mission.target.initialPointECI(3))
view(100,10)
legend(["" "First Stage" "Second Stage" "Capsule" ""])
hold off

%%
%Throttling Plot 
throttling1 = thrustData(:,1,1);
throttling2 = thrustData(:,1,2);

time1 = timeCollocation(:,1);
time2 = timeCollocation(:,2);

tNodes1 = linspace(0,time1(end),5);
tNodes2 = linspace(time2(1),time2(end),5);

throttlingPlot1 = interp1(tNodes1,throttling1,time1);
throttlingPlot2 = interp1(tNodes2,throttling2,time2); 


throttling = figure(3);

plot(time1,throttlingPlot1*100, 'Color',settings.color.terracotta,'LineWidth',2)
hold on
plot(time2,throttlingPlot2*100, 'Color',settings.color.orange,'LineWidth',2)
hold off
ylim([0,100])
xlim([0,time2(end)])
legend(["Throttling Percentage First Stage" "Throttling Percentage Second Stage"])
xlabel("Time")
ylabel("Throttling Percentage")
setPlotSettings(title(""))


%% ANGLE

angle1 = thrustData(:,2,1);
angle2 = thrustData(:,2,2);

time1 = timeCollocation(:,1);
time2 = timeCollocation(:,2);

tNodes1 = linspace(0,time1(end),5);
tNodes2 = linspace(time2(1),time2(end),5);

anglePlot1 = interp1(tNodes1,angle1,time1);
anglePlot2 = interp1(tNodes2,angle2,time2); 


angle = figure(4);

plot(time1,anglePlot1, 'Color',settings.color.terracotta,'LineWidth',2)
hold on
plot(time2,anglePlot2, 'Color',settings.color.orange,'LineWidth',2)
hold off
ylim([0,110])
xlim([0,time2(end)])
legend(["Angle First Stage" "Angle Second Stage"])
xlabel("Time [s]")
ylabel("Angle [°]")
setPlotSettings(title(""))


%% Accelerations

[maxQData]=externalLoads(timeCollocation,stateCollocation,mission,configuration,launcher,mer,0,staging);
acc = [];
for i = 1:launcher(1)
    time      = linspace(timeCollocation(1,i),timeCollocation(end,i),length(timeCollocation)-1);
    v         = stateCollocation(4:6,:,i);
    vNorm     = sqrt( stateCollocation(4,:,i).^2 + stateCollocation(5,:,i).^2 +stateCollocation(6,:,i).^2 );
    accx      = diff(v(1,:))./(time(2)-time(1));
    accy      = diff(v(2,:))./(time(2)-time(1));
    accz      = diff(v(3,:))./(time(2)-time(1));
    acc1       = sqrt( accx.^2 + accy.^2 + accz.^2);
    acc = [acc,acc1];
end



accelerations = figure(5);
plot(timeCollocation(1:99,1),acc(1:99)./mission.environment.g0+1, 'Color',settings.color.terracotta,'LineWidth',2)
hold on
plot(timeCollocation(1:99,2),acc(100:198)./mission.environment.g0+1, 'Color',settings.color.orange,'LineWidth',2)
ylim([0,6])
xlim([0,time2(end)])
legend(["Load Factor First Stage" "Load Factor Second Stage"])
xlabel("Time [s]")
ylabel("Load Factor [-]")
setPlotSettings(title(""))







%%


timeCollocationGT = timeCollocation(:)';
rVec = stateCollocation(1:3,:,:);
rVec = rVec(1:3,:)';

[lonP, latP] = groundTrackFun(timeCollocationGT,rVec,0,0,1,0,titleStringPert,settings);


%%
figure
plot([thrustData(:,1,1);thrustData(:,1,2)])
title("Angle1")



%%




