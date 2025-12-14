
% Initialization

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

[timeCollocationRef, stateCollocationRef] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustDataVecFMC);
%% Create a wind different from nominal to get the gains from the GA

% Initialization of the variables
sizeMC = 1;
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


settings.gaControl = optimoptions("ga", ...
                        "Display","iter", ...
                        "MaxGenerations",20, ...
                        "PopulationSize",50,...
                        "UseParallel",true,...
                        "FunctionTolerance", 1e-4,...
                        "PlotFcn",{'gaplotbestf','gaplotbestindiv'}...
                        );

% set the GA
nVarsGA = launcher(1) * 6;
fun = @(gainGA) objGAGainsMonte(gainGA,launcher,configuration,mission,settings,windVelXFun,windVelYFun,stateCollocationRef,timeCollocationRef,10,thrustData);
lb = zeros(nVarsGA,1);
ub = inf * ones(nVarsGA,1);

% Run GA
[xga,fval,exitFlag,output,population,scores] = ga(fun,nVarsGA,[],[],[],[],lb,ub,[],[],settings.gaControl);

% extrapolate gains
gainGA = xga;

% Number of elements

sizeMC = 1000;

% Initialization of the variables

distanceFromTargetControlled = zeros(1,sizeMC);
cumulativeMeanControlled = zeros(length(distanceFromTargetControlled),1);
distanceFromTargetUncontrolled = zeros(1,sizeMC);
cumulativeMeanUncontrolled = zeros(length(distanceFromTargetUncontrolled),1);
hVec = 0:100;

% Computation of the Wind Profile Population

[meanWind, varWind] = GRAM07_HWM07_annual(hVec);
windUncertainty = sqrt(varWind) .* randn(sizeMC,1);
WindVelocityMag = meanWind + windUncertainty;
windAngVel = WindVelocityMag ./ (mission.environment.rEarth + hVec);
lonInitial = mission.launchBase.lonInitial;
montecarlo.vxWind = - windAngVel .* (mission.environment.rEarth + hVec) .* sin(lonInitial) ;
montecarlo.vyWind = windAngVel .* (mission.environment.rEarth + hVec) .* cos(lonInitial) ;

parfor parforiter = 1:sizeMC
    
    % Functions for wind profile on ECI (rotated inside the dynamics)
    windVelXFun = griddedInterpolant(hVec,montecarlo.vxWind(parforiter,:),'linear','linear');
    windVelYFun = griddedInterpolant(hVec,montecarlo.vyWind(parforiter,:),'linear','linear');

    % Integration of the trajectory
    [timeCollocation, stateCollocation] = totalTrajectoryControlled(launcher,configuration,mission,settings,windVelXFun,windVelYFun,stateCollocationRef,timeCollocationRef,10,thrustData,gainGA);

    % Error
    distanceFromTargetControlled(parforiter) = norm(stateCollocation(1:3,end,end)-mission.target.initialPointECI);

end

% Computation of the cumulative mean

k = 0;

for j = 1:1:length(distanceFromTargetControlled)
    k = k+1;
    cumulativeMeanControlled(k) = mean(distanceFromTargetControlled(1:j));
end

parfor parforiter = 1:sizeMC
    
    % Functions for wind profile on ECI (rotated inside the dynamics)
    windVelXFun = griddedInterpolant(hVec,montecarlo.vxWind(parforiter,:),'linear','linear');
    windVelYFun = griddedInterpolant(hVec,montecarlo.vyWind(parforiter,:),'linear','linear');

    % Integration of the trajectory
    [timeCollocation, stateCollocation] = totalTrajectoryMontecarlo(launcher,configuration,mission,settings,windVelXFun,windVelYFun,thrustData,parforiter);

    % Error
    distanceFromTargetUncontrolled(parforiter) = norm(stateCollocation(1:3,end,end)-mission.target.initialPointECI);

end

% Computation of the cumulative mean

k = 0;

for j = 1:1:length(distanceFromTargetUncontrolled)
    k = k+1;
    cumulativeMeanUncontrolled(k) = mean(distanceFromTargetUncontrolled(1:j));
end



%% ============================== PLOTS ===================================


figure
EarthPlot(mission.environment.rEarth)
hold on
plot3(stateCollocationRef(1,:,1),stateCollocationRef(2,:,1),stateCollocationRef(3,:,1))
plot3(stateCollocationRef(1,:,2),stateCollocationRef(2,:,2),stateCollocationRef(3,:,2))
plot3(stateCollocationRef(1,:,3),stateCollocationRef(2,:,3),stateCollocationRef(3,:,3))
plot3(mission.target.initialPointECI(1),mission.target.initialPointECI(2),mission.target.initialPointECI(3),'ro')
title('Reference Trajectory')
xlabel('X_{ECI}')
ylabel('Y_{ECI}')
zlabel('Z_{ECI}')
%%
% plot of the cumulative mean

figure
subplot(2,1,1)
plot(cumulativeMeanControlled)
title('Cumulative Mean for Controlled Trajectory')
subplot(2,1,2)
plot(cumulativeMeanUncontrolled)
title('Cumulative Mean for Uncontrolled Trajectory')

% Distribution of the 'distance' population

figure
distMeanControlled  = mean(distanceFromTargetControlled);
distStdControlled   = std(distanceFromTargetControlled);
distGaussControlled = linspace(min(distanceFromTargetControlled),max(distanceFromTargetControlled));
distPdfControlled   = normpdf(distGaussControlled, distMeanControlled, distStdControlled);
histogram(distanceFromTargetControlled,'Normalization','pdf','NumBins',2*ceil(log(length(distanceFromTargetControlled))/log(2)+1))
hold on
plot(distGaussControlled, distPdfControlled, 'r', 'LineWidth', 2)
hold off
title('Distribution of the d_{taget} for controlled trajectory')
legend('Real distribution', 'Estimated Gaussian Curve')

figure
distMeanUncontrolled  = mean(distanceFromTargetUncontrolled);
distStdUncontrolled   = std(distanceFromTargetUncontrolled);
distGaussUncontrolled = linspace(min(distanceFromTargetUncontrolled),max(distanceFromTargetUncontrolled));
distPdfUncontrolled   = normpdf(distGaussUncontrolled, distMeanUncontrolled, distStdUncontrolled);
histogram(distanceFromTargetUncontrolled,'Normalization','pdf','NumBins',2*ceil(log(length(distanceFromTargetUncontrolled))/log(2)+1))
hold on
plot(distGaussUncontrolled, distPdfUncontrolled, 'r', 'LineWidth', 2)
hold off
title('Distribution of the d_{taget} for uncontrolled trajectory')
legend('Real distribution', 'Estimated Gaussian Curve')

figure
distancePoints = distanceFromTargetControlled ;
plot(distancePoints,'k*')
yline(20e3,'r','LineWidth',2)
ylabel('Distance')
xlabel('iteration')

distPointGood =[]; %NaN * ones(length(distancePoints));
distPointNotGood =[]; % NaN * ones(length(distancePoints));
for k = 1:length(distancePoints)
    if distancePoints(k) <= 20e3
        distPointGood = [distPointGood distancePoints(k)];
    else
        distPointNotGood = [distPointNotGood distancePoints(k)];
    end
end
percIn = numel(distPointGood)/numel(distancePoints) * 100;
percOut = numel(distPointNotGood)/numel(distancePoints) * 100;
% figure
% plot(distPointGood,'g*')
% hold on
% plot(distPointNotGood,'r*')
% yline(20e3,'r')
% ylabel('Distance')
% xlabel('iteration')

% figure; clf
% b = bar([percIn percOut], 'stacked');
% 

% 
% ylim([0 100]);
% ylabel('[%]');
% title('Values inside/outside the range');

figure; clf
y = [percIn; percOut];              % 2 valori separati

x = categorical({'success','failure'});
x = reordercats(x,{'success','failure'});   % opzionale, per l'ordine

b = bar(x, y);
b.FaceColor = 'flat';          % abilita colori per-segmento
b.CData(1,:) = [0 0.8 0];      % primo pezzo (percIn)  -> verde
b.CData(2,:) = [0.8 0 0];      % secondo pezzo (percOut) -> rosso
ylabel('[%]');
title('Success vs Failure');



%% Propagation of the other uncertainties

sizeMC = 5;


for i = 1:launcher(1)
    configuration.stage{i}.engine = mission.engines{launcher(1+i)};
end

[mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);

meanStructMass1 = staging{1}.mStruct;
meanStructMass2 = staging{2}.mStruct;
meanPropMass1 = staging{1}.mProp;
meanPropMass2 = staging{2}.mProp;
meanIsp1 = configuration.stage{1}.engine.ispZero;
meanIsp2 = configuration.stage{2}.engine.ispVac;

structMass1Distrib = meanStructMass1.*(1 + 0.01*randn(sizeMC,1));
structMass2Distrib = meanStructMass2.*(1 + 0.01*randn(sizeMC,1));
propMass1Distrib = meanPropMass1.*(1 + 0.01*randn(sizeMC,1));
propMass2Distrib = meanPropMass2.*(1 + 0.01*randn(sizeMC,1));
isp1Distrib = meanIsp1.*(1 + 0.001*randn(sizeMC,1));
isp2Distrib = meanIsp2.*(1 + 0.001*randn(sizeMC,1));

mStage1Distrib = structMass1Distrib + propMass1Distrib;
mStage2Distrib = structMass2Distrib + propMass2Distrib;

[A,B,C,D,E,F] = ndgrid(structMass1Distrib, structMass2Distrib,propMass1Distrib,propMass2Distrib,isp1Distrib,isp2Distrib);
Matr =[A(:),B(:),C(:),D(:),E(:),F(:)];

% Shuffle of the Matrix elements
for i=1:size(Matr,2)
    shuffledMatrix(:,i) = Matr(randperm(size(Matr,1)),i);
end

distanceFromTargetControlled = zeros(size(shuffledMatrix,1),1);

% Initialization of the variables
sizeMC = 1;
hVec = 0:100;

% Computation of the Wind Profiles

[meanWind, varWind] = GRAM07_HWM07_annual(hVec);
WindVelocityMag = meanWind ;
windAngVel = WindVelocityMag ./ (mission.environment.rEarth + hVec);
lonInitial = mission.launchBase.lonInitial;
montecarlo.vxWind = - windAngVel .* (mission.environment.rEarth + hVec) .* sin(lonInitial) ;
montecarlo.vyWind = windAngVel .* (mission.environment.rEarth + hVec) .* cos(lonInitial) ;

% Functions for wind profile on ECI (rotated inside the dynamics)
windVelXFun = griddedInterpolant(hVec,montecarlo.vxWind(1,:),'linear','linear');
windVelYFun = griddedInterpolant(hVec,montecarlo.vyWind(1,:),'linear','linear');
settings.gaControl = optimoptions("ga", ...
                        "Display","iter", ...
                        "MaxGenerations",20, ...
                        "PopulationSize",100,...
                        "UseParallel",true,...
                        "FunctionTolerance", 1e-4,...
                        "PlotFcn",{'gaplotbestf','gaplotbestindiv'}...
                        );

% set the GA
nVarsGA = launcher(1) * 6;
fun = @(gainGA) objGAGainsMonte(gainGA,launcher,configuration,mission,settings,windVelXFun,windVelYFun,stateCollocationRef,timeCollocationRef,10,thrustData);
lb = zeros(nVarsGA,1);
ub = inf * ones(nVarsGA,1);

figure
EarthPlot(mission.environment.rEarth)
hold on
plot3(mission.target.initialPointECI(1),mission.target.initialPointECI(2),mission.target.initialPointECI(3),'bo')
[timeCollocationRef, stateCollocationRef] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustDataVecFMC);



for i=1:size(shuffledMatrix,1)

    configuration.stage{1}.mProp = shuffledMatrix(i,3);
    configuration.stage{2}.mProp = shuffledMatrix(i,4);
    configuration.stage{1}.mStage = shuffledMatrix(i,1) + shuffledMatrix(i,3);
    configuration.stage{2}.mStage = shuffledMatrix(i,2) + shuffledMatrix(i,4);
    configuration.stage{1}.engine.ispZero = shuffledMatrix(i,5);
    configuration.stage{2}.engine.ispVac = shuffledMatrix(i,6);

    % [timeCollocationRef, stateCollocationRef] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustDataVecFMC);

    % Tuning the gains

    % if i==1
    %     % set the GA
    % nVarsGA = launcher(1) * 6;
    % fun = @(gainGA) objGAGainsMonte(gainGA,launcher,configuration,mission,settings,windVelXFun,windVelYFun,stateCollocationRef,timeCollocationRef,10,thrustData);
    % lb = zeros(nVarsGA,1);
    % ub = inf * ones(nVarsGA,1);
    % 
    % % Run GA
    % [xga,fval,exitFlag,output,population,scores] = ga(fun,nVarsGA,[],[],[],[],lb,ub,[],[],settings.gaControl);
    % 
    % % extrapolate gains
    % gainGA = xga;
    % end

    % Integration of the trajectory
    [timeCollocation, stateCollocation] = totalTrajectoryControlled(launcher,configuration,mission,settings,windVelXFun,windVelYFun,stateCollocationRef,timeCollocationRef,10,thrustData,gainGA);

    % Error
    distanceFromTargetControlled(i) = norm(stateCollocation(1:3,end,end)-mission.target.initialPointECI);
    
    plot3(stateCollocation(1,end,end),stateCollocation(2,end,end),stateCollocation(3,end,end),'ro')
    drawnow
end