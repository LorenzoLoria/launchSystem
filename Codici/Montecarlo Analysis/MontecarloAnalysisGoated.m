
% Initialization

clear all
clc
close all

addpath(genpath("..\..\"))

[mission,settings] = dataStructGlobal;

launcher = [2,2,3,4,0.47,0.5,0.7];

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

% Nominal Trajectory
[timeCollocationRef, stateCollocationRef] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustData);
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

maxThrust = [mission.engines{1}.thrustZero*4 ; mission.engines{4}.thrustVacum];

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
fun = @(gainGA) objGAGainsMonte(gainGA,launcher,configuration,mission,settings,windVelXFun,windVelYFun,stateCollocationRef,timeCollocationRef,maxThrust,10,thrustData);
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
    [timeCollocation, stateCollocation] = totalTrajectoryControlled(launcher,configuration,mission,settings,windVelXFun,windVelYFun,stateCollocationRef,timeCollocationRef,maxThrust,10,thrustData,gainGA);

    % Error
    distanceFromTargetControlled(parforiter) = norm(stateCollocation(1:3,end,end)-mission.target.initialPointECI);

end

% Computation of the cumulative mean

k = 0;

for j = 1:1:length(distanceFromTargetControlled)
    k = k+1;
    cumulativeMeanControlled(k) = mean(distanceFromTargetControlled(1:j));
end

%%

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
