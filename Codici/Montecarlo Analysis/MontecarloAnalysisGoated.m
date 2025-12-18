clear all
clc
close all

% Path Directory

addpath(genpath("..\..\"))

% Upload Mission Struct
[mission,settings] = dataStructGlobal;
% Initialization

launcher = [2.0	2.0	3.0	3.0	0.5718839104841471	0.8679383834568386	0.6027644098891288];
for i = 1:launcher(1)
    configuration.stage{i}.engine = mission.engines{launcher(1+i)};
end

[mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);

 thrustDataVecFMC = [   0.997027764065679
   0.942654752527805
   0.951101595174460
   0.955848071611001
   0.958803211573321
   0.266317283506581
  18.826711728279765
  32.515170831440209
  33.514016831344129
  46.141997578498405
   0.999087083722307
   0.853687176825034
   0.697414992914764
   0.646557363727949
   0.912180915015670
  55.283673653565842
  86.115064742122868
  81.590553439030884
  89.388813826003727
  98.329715236102246];

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
%%

close all
% extrapolate gains
gainGA = xga;

% Number of elements

sizeMC = 100;

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

distMeanCon = mean(distanceFromTargetControlled);
distStdCon = std(distanceFromTargetControlled);

distMeanUnc = mean(distanceFromTargetUncontrolled);
distStdUnc = std(distanceFromTargetUncontrolled);


%% ============================== PLOTS ===================================


NominalTrajectory = figure(1);
EarthPlot(mission.environment.rEarth)
hold on
plot3(stateCollocationRef(1,:,1),stateCollocationRef(2,:,1),stateCollocationRef(3,:,1))
plot3(stateCollocationRef(1,:,2),stateCollocationRef(2,:,2),stateCollocationRef(3,:,2))
plot3(stateCollocationRef(1,:,3),stateCollocationRef(2,:,3),stateCollocationRef(3,:,3))
plot3(mission.target.initialPointECI(1),mission.target.initialPointECI(2),mission.target.initialPointECI(3),'ro')
xlabel('X_{ECI}')
ylabel('Y_{ECI}')
zlabel('Z_{ECI}')
setPlotSettings(title('Reference Trajectory'))

exportStandardizedFigure(NominalTrajectory,'NominalTrajectory',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')


% plot of the cumulative mean

ComparisonMean=figure(2);
subplot(2,1,1)
plot(cumulativeMeanControlled,'Color',settings.color.blu)
setPlotSettings(title('Cumulative Mean for Controlled Trajectory'))
subplot(2,1,2)
plot(cumulativeMeanUncontrolled,'Color',settings.color.blu)
setPlotSettings(title('Cumulative Mean for Uncontrolled Trajectory'))

exportStandardizedFigure(ComparisonMean,'ComparisonMean',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')

% Distribution of the 'distance' population

% figure
% distMeanControlled  = mean(distanceFromTargetControlled);
% distStdControlled   = std(distanceFromTargetControlled);
% distGaussControlled = linspace(min(distanceFromTargetControlled),max(distanceFromTargetControlled));
% distPdfControlled   = normpdf(distGaussControlled, distMeanControlled, distStdControlled);
% histogram(distanceFromTargetControlled,'Normalization','pdf','NumBins',2*ceil(log(length(distanceFromTargetControlled))/log(2)+1))
% hold on
% plot(distGaussControlled, distPdfControlled, 'r', 'LineWidth', 2)
% hold off
% title('Distribution of the d_{taget} for controlled trajectory')
% legend('Real distribution', 'Estimated Gaussian Curve')
% 
% figure
% distMeanUncontrolled  = mean(distanceFromTargetUncontrolled);
% distStdUncontrolled   = std(distanceFromTargetUncontrolled);
% distGaussUncontrolled = linspace(min(distanceFromTargetUncontrolled),max(distanceFromTargetUncontrolled));
% distPdfUncontrolled   = normpdf(distGaussUncontrolled, distMeanUncontrolled, distStdUncontrolled);
% histogram(distanceFromTargetUncontrolled,'Normalization','pdf','NumBins',2*ceil(log(length(distanceFromTargetUncontrolled))/log(2)+1))
% hold on
% plot(distGaussUncontrolled, distPdfUncontrolled, 'r', 'LineWidth', 2)
% hold off
% title('Distribution of the d_{taget} for uncontrolled trajectory')
% legend('Real distribution', 'Estimated Gaussian Curve')

Points = figure(3);
distancePointsControlled = distanceFromTargetControlled ;
plot(distancePointsControlled,'k*')
yline(20e3,'r','LineWidth',2)
ylabel('Distance')
xlabel('iteration')


distPointGood =[]; %NaN * ones(length(distancePoints));
distPointNotGood =[]; % NaN * ones(length(distancePoints));
for k = 1:length(distancePointsControlled)
    if distancePointsControlled(k) <= 20e3
        distPointGood = [distPointGood distancePointsControlled(k)];
    else
        distPointNotGood = [distPointNotGood distancePointsControlled(k)];
    end
end
percIn = numel(distPointGood)/numel(distancePointsControlled) * 100;
percOut = numel(distPointNotGood)/numel(distancePointsControlled) * 100;
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

succesvFailControl=figure(4); clf
y = [percIn; percOut];              % 2 valori separati

x = categorical({'success','failure'});
x = reordercats(x,{'success','failure'});   % opzionale, per l'ordine

b = bar(x, y);
b.FaceColor = 'flat';          % abilita colori per-segmento
b.CData(1,:) = [0 0.8 0];      % primo pezzo (percIn)  -> verde
b.CData(2,:) = [0.8 0 0];      % secondo pezzo (percOut) -> rosso
ylabel('[%]');

setPlotSettings(title('Success vs Failure'))

exportStandardizedFigure(succesvFailControl,'succesvFailControl',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')


PointsUnc = figure(5);
distancePointsUnontrolled = distanceFromTargetUncontrolled ;
plot(distancePointsUnontrolled,'k*')
yline(20e3,'r','LineWidth',2)
ylabel('Distance')
xlabel('iteration')

distPointGood =[]; %NaN * ones(length(distancePoints));
distPointNotGood =[]; % NaN * ones(length(distancePoints));
for k = 1:length(distancePointsUnontrolled)
    if distancePointsUnontrolled(k) <= 20e3
        distPointGood = [distPointGood distancePointsUnontrolled(k)];
    else
        distPointNotGood = [distPointNotGood distancePointsUnontrolled(k)];
    end
end
percIn = numel(distPointGood)/numel(distancePointsUnontrolled) * 100;
percOut = numel(distPointNotGood)/numel(distancePointsUnontrolled) * 100;

succesvFailUncontrol=figure(6); clf
y = [percIn; percOut];              % 2 valori separati

x = categorical({'success','failure'});
x = reordercats(x,{'success','failure'});   % opzionale, per l'ordine

b = bar(x, y);
b.FaceColor = 'flat';          % abilita colori per-segmento
b.CData(1,:) = [0 0.8 0];      % primo pezzo (percIn)  -> verde
b.CData(2,:) = [0.8 0 0];      % secondo pezzo (percOut) -> rosso
ylabel('[%]');

setPlotSettings(title('Success vs Failure'))

exportStandardizedFigure(succesvFailUncontrol,'succesvFailUncontrol',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')


%% Propagation of the other uncertainties (stage 1)
clear all
clc
close all

% Path Directory

addpath(genpath("..\..\"))

% Upload Mission Struct
[mission,settings] = dataStructGlobal;
% Initialization

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

close all

sizeMC = 30;


for i = 1:launcher(1)
    configuration.stage{i}.engine = mission.engines{launcher(1+i)};
end

[mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);

meanStructMass1 = staging{1}.mStruct;
%meanStructMass2 = staging{2}.mStruct;
meanPropMass1 = staging{1}.mProp;
%meanPropMass2 = staging{2}.mProp;
meanIsp1 = configuration.stage{1}.engine.ispZero;
%meanIsp2 = configuration.stage{2}.engine.ispVac;

structMass1Distrib = meanStructMass1.*(1 + 0.01*randn(sizeMC,1));
%structMass2Distrib = meanStructMass2.*(1 + 0.01*randn(sizeMC,1));
propMass1Distrib = meanPropMass1.*(1 + 0.01*randn(sizeMC,1));
%propMass2Distrib = meanPropMass2.*(1 + 0.01*randn(sizeMC,1));
isp1Distrib = meanIsp1.*(1 + 0.001*randn(sizeMC,1));
%isp2Distrib = meanIsp2.*(1 + 0.001*randn(sizeMC,1));

mStage1Distrib = structMass1Distrib + propMass1Distrib;
%mStage2Distrib = structMass2Distrib + propMass2Distrib;

[A,B,C] = ndgrid(structMass1Distrib,propMass1Distrib,isp1Distrib);
Matr =[A(:),B(:),C(:)];

% Shuffle of the Matrix elements

shuffledMatrix = zeros(size(Matr));

for i=1:size(Matr,2)
    shuffledMatrix(:,i) = Matr(randperm(size(Matr,1)),i);
end

positionError = zeros(3,size(shuffledMatrix,1));
velocityError = zeros(3,size(shuffledMatrix,1));

% Initialization of the variables
sizeMC = 1;
hVec = 0:100;


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

for i=1:size(shuffledMatrix,1)

    configuration.stage{1}.mProp = shuffledMatrix(i,2);
    %configuration.stage{2}.mProp = shuffledMatrix(i,4);
    configuration.stage{1}.mStage = shuffledMatrix(i,1) + shuffledMatrix(i,2);
    %configuration.stage{2}.mStage = shuffledMatrix(i,2) + shuffledMatrix(i,4);
    configuration.stage{1}.engine.ispZero = shuffledMatrix(i,3);
    %configuration.stage{2}.engine.ispVac = shuffledMatrix(i,6);

    

    % %Tuning the gains
    % 
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
    [timeCollocation, stateCollocation] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustDataVecFMC);

    % Error
    positionError(:,i) = stateCollocation(1:3,end,end-1)-stateCollocationRef(1:3,end,end-1);
    velocityError(:,i) = stateCollocation(4:6,end,end-1)-stateCollocationRef(4:6,end,end-1);
end


meanErrPosStage1 = mean(positionError);
stdErrPosStage1 = std(positionError);

meanErrVelStage1 = mean(velocityError);
stdErrVelStage1 = std(velocityError);

histPosx1 = figure(1);
histogram(positionError(1,:),'FaceColor',settings.color.blu,'EdgeColor','k','NumBins',30)
%xlim([-45000,45000])
setPlotSettings(title('Position error along $x_{ECI}$ [m]'))
exportStandardizedFigure(histPosx1,'histPosx1',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')

histPosy1 = figure(2);
histogram(positionError(2,:),'FaceColor',settings.color.blu,'EdgeColor','k','NumBins',30)
%xlim([-45000,45000])
setPlotSettings(title('Position error along $y_{ECI}$ [m]'))
exportStandardizedFigure(histPosy1,'histPosy1',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')

histPosz1 = figure(3);
histogram(positionError(3,:),'FaceColor',settings.color.blu,'EdgeColor','k','NumBins',30)
%xlim([-45000,45000])
setPlotSettings(title('Position error along $z_{ECI}$ [m]'))
exportStandardizedFigure(histPosz1,'histPosz1',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')

histVelx1 = figure(4);
histogram(velocityError(1,:),'FaceColor',settings.color.blu,'EdgeColor','k','NumBins',30)
%xlim([-250,250])
setPlotSettings(title('Velocity error along $x_{ECI}$ [m/s]'))
exportStandardizedFigure(histVelx1,'histVelx1',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')

histVely1 = figure(5);
histogram(velocityError(2,:),'FaceColor',settings.color.blu,'EdgeColor','k','NumBins',30)
%xlim([-250,250])
setPlotSettings(title('Velocity error along $y_{ECI}$ [m/s]'))
exportStandardizedFigure(histVely1,'histVely1',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')

histVelz1 = figure(6);
histogram(velocityError(3,:),'FaceColor',settings.color.blu,'EdgeColor','k','NumBins',30)
%xlim([-250,250])
setPlotSettings(title('Velocity error along $z_{ECI}$ [m/s]'))
exportStandardizedFigure(histVelz1,'histVelz1',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')

% %% Computation of the cumulative mean
% 
% k = 0;
% 
% cumulativeMeanPosition = zeros(length(positionError),1);
% cumulativeMeanVelocity = zeros(length(velocityError),1);
% 
% for j = 1:1:length(positionError)
%     k = k+1;
%     cumulativeMeanPosition(j) = mean(positionError(1:j));
%     cumulativeMeanVelocity(j) = mean(velocityError(1:j));
% end
% positionErrorPlotStage1 = figure(1);
% plot(cumulativeMeanPosition,'Color',settings.color.blu)
% setPlotSettings(title('Cumulative Mean for the position error'))
% exportStandardizedFigure(positionErrorPlotStage1,'positionErrorPlotStage1',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')
% 
% velocityErrorPlotStage1 = figure(2);
% plot(cumulativeMeanVelocity,'Color',settings.color.blu)
% setPlotSettings(title('Cumulative Mean for the velocity error'))
% exportStandardizedFigure(velocityErrorPlotStage1,'velocityErrorPlotStage1',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')
% 

% Propagation of the other uncertainties (stage 2)
close all
sizeMC = 30;


for i = 1:launcher(1)
    configuration.stage{i}.engine = mission.engines{launcher(1+i)};
end

[mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);

%meanStructMass1 = staging{1}.mStruct;
meanStructMass2 = staging{2}.mStruct;
%meanPropMass1 = staging{1}.mProp;
meanPropMass2 = staging{2}.mProp;
%meanIsp1 = configuration.stage{1}.engine.ispZero;
meanIsp2 = configuration.stage{2}.engine.ispVac;

%structMass1Distrib = meanStructMass1.*(1 + 0.01*randn(sizeMC,1));
structMass2Distrib = meanStructMass2.*(1 + 0.01*randn(sizeMC,1));
%propMass1Distrib = meanPropMass1.*(1 + 0.01*randn(sizeMC,1));
propMass2Distrib = meanPropMass2.*(1 + 0.01*randn(sizeMC,1));
%isp1Distrib = meanIsp1.*(1 + 0.001*randn(sizeMC,1));
isp2Distrib = meanIsp2.*(1 + 0.001*randn(sizeMC,1));

%mStage1Distrib = structMass1Distrib + propMass1Distrib;
mStage2Distrib = structMass2Distrib + propMass2Distrib;

[A,B,C] = ndgrid(structMass2Distrib,propMass2Distrib,isp2Distrib);
Matr =[A(:),B(:),C(:)];

% Shuffle of the Matrix elements

shuffledMatrix = zeros(size(Matr));

for i=1:size(Matr,2)
    shuffledMatrix(:,i) = Matr(randperm(size(Matr,1)),i);
end

positionError = zeros(3,size(shuffledMatrix,1),1);
velocityError = zeros(3,size(shuffledMatrix,1),1);

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

[timeCollocationRef, stateCollocationRef] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustDataVecFMC);

for i=1:size(shuffledMatrix,1)

    %configuration.stage{1}.mProp = shuffledMatrix(i,2);
    configuration.stage{2}.mProp = shuffledMatrix(i,2);
    %configuration.stage{1}.mStage = shuffledMatrix(i,1) + shuffledMatrix(i,2);
    configuration.stage{2}.mStage = shuffledMatrix(i,1) + shuffledMatrix(i,2);
    %configuration.stage{1}.engine.ispZero = shuffledMatrix(i,3);
    configuration.stage{2}.engine.ispVac = shuffledMatrix(i,3);

    

    % %Tuning the gains
    % 
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
    [timeCollocation, stateCollocation] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustDataVecFMC);

    % Error
    positionError(:,i) = stateCollocation(1:3,end,end-1)-stateCollocationRef(1:3,end,end-1);
    velocityError(:,i) = stateCollocation(4:6,end,end-1)-stateCollocationRef(4:6,end,end-1);


end

meanErrPosStage2 = mean(positionError);
stdErrPosStage2 = std(positionError);

meanErrVelStage2 = mean(velocityError);
stdErrVelStage2 = std(velocityError);


histPosx2 = figure(1);
histogram(positionError(1,:),'FaceColor',settings.color.blu,'EdgeColor','k','NumBins',30)
%xlim([-1.5e4,1.5e4])
setPlotSettings(title('Position error along $x_{ECI} $[m]'))
exportStandardizedFigure(histPosx2,'histPosx2',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')

histPosy2 = figure(2);
histogram(positionError(2,:),'FaceColor',settings.color.blu,'EdgeColor','k','NumBins',30)
%xlim([-6e4,6e4])
setPlotSettings(title('Position error along $y_{ECI} $[m]'))
exportStandardizedFigure(histPosy2,'histPosy2',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')

histPosz2 = figure(3);
histogram(positionError(3,:),'FaceColor',settings.color.blu,'EdgeColor','k','NumBins',30)
%xlim([-1.5e4,1.5e4])
setPlotSettings(title('Position error along $z_{ECI} $ [m]'))
exportStandardizedFigure(histPosz2,'histPosz2',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')

histVelx2 = figure(4);
histogram(velocityError(1,:),'FaceColor',settings.color.blu,'EdgeColor','k','NumBins',30)
%xlim([-60,60])
setPlotSettings(title('Velocity error along $x_{ECI} $[m/s]'))
exportStandardizedFigure(histVelx2,'histVelx2',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')

histVely2 = figure(5);
histogram(velocityError(2,:),'FaceColor',settings.color.blu,'EdgeColor','k','NumBins',30)
%xlim([-200,200])
setPlotSettings(title('Velocity error along $y_{ECI}$ [m/s]'))
exportStandardizedFigure(histVely2,'histVely2',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')

histVelz2 = figure(6);
histogram(velocityError(3,:),'FaceColor',settings.color.blu,'EdgeColor','k','NumBins',30)
%xlim([-60,60])
setPlotSettings(title('Velocity error along $z_{ECI}$ [m/s]'))
exportStandardizedFigure(histVelz2,'histVelz2',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')

%% Computation of the cumulative mean

k = 0;

cumulativeMeanPosition = zeros(length(positionError),1);
cumulativeMeanVelocity = zeros(length(velocityError),1);

for j = 1:1:length(positionError)
    k = k+1;
    cumulativeMeanPosition(j) = mean(positionError(1:j));
    cumulativeMeanVelocity(j) = mean(velocityError(1:j));
end
positionErrorPlotStage2 = figure(1);
plot(cumulativeMeanPosition,'Color',settings.color.blu)
setPlotSettings(title('Cumulative Mean for the position error'))
exportStandardizedFigure(positionErrorPlotStage2,'positionErrorPlotStage2',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')

velocityErrorPlotStage2 = figure(2);
plot(cumulativeMeanVelocity,'Color',settings.color.blu)
setPlotSettings(title('Cumulative Mean for the velocity error'))
exportStandardizedFigure(velocityErrorPlotStage2,'velocityErrorPlotStage2',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')


%% Influence of the errors in position and velocity wrt trajectory height
clear all
clc
close all

addpath(genpath("..\..\"))

[mission,settings] = dataStructGlobal;

launcher = [2.0	2.0	3.0	3.0	0.459952176990556	0.7543705311589043	0.6347957418855592;...
2.0	4.0	3.0	3.0	0.546875	0.8812500000000001	0.43125;...
3.0	2.0	3.0	3.0	0.5428936903202396	0.3086421045633821	0.7967415421952768;...
3.0	4.0	3.0	3.0	0.5892783340106565	0.7402654405122346	0.613707074894946;...
3.0	2.0	3.0	3.0	0.49766283303833814	0.6632650657874963	0.8773977802975306;...
2.0	2.0	3.0	3.0	0.5367409770035041	0.8290068052727061	0.8663007317366302;...
3.0	2.0	3.0	3.0	0.4973304292947056	0.40468340649728135	0.8814431147582379;...
3.0	2.0	3.0	3.0	0.5428936903202396	0.3086421045633821	0.7967415421952768;...
2.0	2.0	3.0	3.0	0.5708877449876226	0.8001322138002916	0.7375746684496541;...
2.0	2.0	3.0	3.0	0.5367409770035041	0.8290068052727061	0.8663007317366302;...
2.0	2.0	3.0	3.0	0.40545722914246285	0.8652055699518671	0.7882838965431354;...
2.0	4.0	3.0	3.0	0.6011668795260986	0.7819774421784801	0.8158827144968238;...
3.0	2.0	3.0	3.0	0.5400291127513209	0.4942159883194645	0.847831302218087;...
3.0	2.0	3.0	3.0	0.5994066768849554	0.3700755287295609	0.7429660545310444;...
2.0	2.0	3.0	3.0	0.513465673550142	0.7850011366362578	0.8682316772357497;...
2.0	2.0	3.0	3.0	0.5332321562661497	0.862063900613571	0.8760583948911163;...
2.0	2.0	3.0	3.0	0.5225772656447571	0.8141900136416045	0.7688846153271951;...
2.0	2.0	3.0	3.0	0.5343961934728245	0.7764568902527008	0.8680663839021251;...
2.0	2.0	3.0	3.0	0.5107214190869522	0.7611369258076272	0.6763840990319829;...
2.0	2.0	3.0	3.0	0.5225772656447571	0.8141900136416045	0.7688846153271951;...
2.0	2.0	3.0	3.0	0.5387270219549453	0.8250583722152607	0.8761683568795317;...
2.0	2.0	3.0	3.0	0.5107214190869522	0.7611369258076272	0.6763840990319829;...
2.0	2.0	3.0	3.0	0.5996262334097416	0.7803193698270777	0.8851982317179204;...
3.0	4.0	3.0	3.0	0.5159912736203377	0.3588443181659504	0.8832452760621861;...
2.0	2.0	3.0	3.0	0.5260509367911556	0.8741076026398458	0.6762532795731154;...
2.0	2.0	3.0	3.0	0.5550786187243343	0.8652248740324553	0.6749645278820562;...
2.0	2.0	3.0	3.0	0.5260509367911556	0.8741076026398458	0.6762532795731154;...
2.0	2.0	3.0	3.0	0.5327556637915469	0.8839945121998461	0.4281023970511908;...
2.0	2.0	3.0	3.0	0.41132525337528275	0.8306853204956812	0.8149820366325383;...
2.0	2.0	3.0	3.0	0.5107214198613309	0.7611369265796956	0.6763840969393151;...
2.0	2.0	3.0	3.0	0.43444917856433163	0.7862352390623101	0.7075833080800635;...
2.0	2.0	3.0	3.0	0.5255718332549055	0.7532192099047599	0.6762709355836022;...
2.0	2.0	3.0	3.0	0.5075069937717044	0.81908503500597	0.780319721323177;...
2.0	2.0	3.0	3.0	0.43437499999999996	0.8062500000000001	0.35625;...
2.0	4.0	3.0	3.0	0.546875	0.8812500000000001	0.43125;...
2.0	4.0	3.0	3.0	0.546875	0.8812500000000001	0.43125;...
2.0	2.0	3.0	3.0	0.5422168405107969	0.8270366462834824	0.41434291552917174;...
2.0	4.0	3.0	3.0	0.5468750020953832	0.8812499881896585	0.43124999733314867;...
2.0	2.0	3.0	3.0	0.5502792477885152	0.8455405439139619	0.40775613838869174;...
2.0	2.0	3.0	3.0	0.5502792477885152	0.8455405439139619	0.40775613838869174;...
3.0	2.0	3.0	3.0	0.528125	0.39375	0.54375;...
2.0	2.0	3.0	3.0	0.5722473975667838	0.8575377808140997	0.384017819782447;...
3.0	2.0	3.0	3.0	0.5207113576743595	0.4525869676416344	0.7640586641234791;...
2.0	2.0	3.0	3.0	0.5722473975667838	0.8575377808140997	0.384017819782447;...
3.0	2.0	3.0	3.0	0.5114620374796529	0.48919456475634054	0.8289300766789469;...
3.0	4.0	3.0	3.0	0.565625	0.46875	0.61875;...
3.0	4.0	3.0	3.0	0.5660495703033552	0.46790085939328957	0.6201652343445174;...
2.0	2.0	3.0	3.0	0.5718839104841471	0.8679383834568386	0.6027644098891288;...
3.0	4.0	3.0	3.0	0.5660495703033552	0.46790085939328957	0.6201652343445174;...
2.0	2.0	3.0	3.0	0.5718839104841471	0.8679383834568386	0.6027644098891288;...
2.0	4.0	3.0	3.0	0.5379794858812149	0.7614574415209188	0.6432887903432223;...
3.0	4.0	3.0	3.0	0.5296141538502213	0.42330089425881895	0.6964242939412005;...
2.0	2.0	3.0	3.0	0.5718839104841471	0.8679383834568386	0.6027644098891288];


thrustData = [0.901784296794966	0.9557332507188437	0.959077988241489	0.9854897637057565	0.9879060788152333	0.9982340890548043	0.9787668403172282	0.9965081863412928	0.9487423275062546	0.9487404730639827	0.9997344613670269	0.9558152857779411	0.97834377707324	0.9973193938418461	0.948742327073947	0.9996601515722676	0.9487423275062546	0.9487423270566042	0.9487423275062546	0.982819810725315	0.9487408474459559	0.9865663062865458	0.9487423275063142	0.9871276701051527	0.9487423275062546	0.9487423275062546	0.9487423275064402	0.9487421713615554	0.9998778016315012	0.9651301937335665	0.9825669710191074	0.9959102227151216	0.9392402217668435	0.9785572366925552	0.9487423275062546	0.9487423100509781	0.9997785385962358	0.9487423275062546	0.9487423274736503	0.9487423148080698	0.9788794242513349	0.9803644629318801	0.9692922585082769	0.9487423275062536	0.9616910579682059	1.0	0.9407367152187078	0.9921898476149137	0.9557891864821773	0.9785032458706234	0.9488259826251969	0.999154916538961	0.997027764065679;...
0.9999674580286427	0.9425779013850505	0.9467276994411381	0.9677575471628994	1.0	0.9902761116692427	0.9444315574973697	0.9720385014969348	0.9277908817456332	0.9902762900468526	0.9999999882289988	0.9988639212957003	0.9349392018420066	0.9833485773638354	0.9902763197011585	0.9810802881918627	0.9902763197161293	0.9988146170535571	0.9902763197161293	0.9829512069471825	0.9902762640999196	0.9902763197161293	0.9546760241798989	0.9186681470945937	0.9902763197161293	0.9902763197161293	0.9902763197161364	0.9902763090056965	0.9000048432929755	0.9503371097410461	0.9606030333168889	0.975404664403355	0.9892788657084552	0.9052226684032352	0.9905204603411293	0.9106661450143029	0.9998881780125194	0.990295288894329	0.9902763196961125	0.9902763183057867	0.9376421022483573	0.9902763197190582	0.9961274035773754	0.9282627603580484	0.9469561357472137	1.0	0.9136668331165754	0.9902763170033021	0.9410631909374766	0.990276316360522	0.9422249550089791	0.9698043666331716	0.9426547525278047;...
0.9000025616008186	0.9000009685755028	0.9227563282294452	0.9	0.9772384307393943	0.9898055931576069	0.9	0.9567900628692967	0.9516021300924457	0.9818894091811315	0.9998582199416037	0.9	0.9825118375549773	0.9674887988644236	0.9888985075718449	0.9669978063652296	0.9930852199117732	0.9	0.9559792307206758	0.9785564407180731	0.9743009631523325	0.9942933668881295	0.9000000596046448	0.9658633300466226	0.9937577015881497	0.9849717773069264	0.9	0.9613842300644535	0.9856785276009491	0.9917003251467494	0.9674058552375464	0.9474464617771119	0.900000001654897	0.9991994862503942	0.9699332420707425	0.9	0.9336737611353957	0.9	0.980500024176705	0.996330385941665	0.9704880179300459	0.9801431200365865	0.9364834189177433	0.9009841918945313	0.9645448934266782	1.0	0.9082655076433312	0.973695585147831	0.9997480061314544	0.9000000149011612	0.9	0.9416841583899985	0.9511015951744598;...
0.9000000000000067	0.9715966623931802	0.9858852729263722	0.9115859798541472	1.0	0.9	0.9401547376495807	0.929911689503517	0.9543181008247372	0.9	0.9998874920689337	0.9663446894990211	0.9549702717687348	0.9597755008626173	0.9000000596046448	0.9813672878743042	0.9886798980295173	0.9009765624998192	0.9897679852472069	0.9000190734864733	0.9317891416286326	0.9128070122014705	0.9268103598207652	0.9601055214613236	0.9866592497658307	0.9427456212224672	0.9731916543051815	0.9602559703827449	0.989456834409506	0.9472032279100113	0.9867936900339188	0.9873588146201034	0.9000000016943238	0.9803736690338773	0.9000038146972656	0.9	0.9934864008644843	0.9	0.9993605604668871	0.9953864302638927	0.9835244195379267	0.9523567814984778	0.9240727721005719	0.9	0.9249861399828703	0.9843723374957625	0.9248230540047486	0.9935184832549608	0.9037164738535813	0.9405544152665621	0.9	0.9998884738124039	0.9558480716110014;...
0.9000000000000058	0.9	0.9015480885071626	0.956302907487519	0.9554888429254983	0.9835132620989395	0.9032919080926147	0.9062344655817036	0.90006103515625	0.9	0.9998950215400455	0.9	0.9104545495035388	0.9423803333705523	0.9000010132789612	0.9533355306312766	0.9867663147777839	0.9993230413928806	0.915625	0.9679956030637371	0.9	0.9747077583855126	0.9	0.9659164680713493	0.9	0.9002443790435791	0.9	0.9	0.9392689661186178	0.9406146483295297	0.9000000000000002	0.9508406711690679	0.9000000009528917	0.9585170208726972	0.9479207287246911	0.9000038146972686	0.9000786152982689	0.9000038146972656	0.9000152587893196	0.9859550574438751	0.9556761064802887	0.900061035156244	0.9101160821536949	0.9	0.9106932576091417	0.9000153183937291	0.9459782599471976	0.9002441406247069	0.934281782489719	0.900000968575478	0.9000002384185791	0.9	0.9588032115733206;...
1.323764277511055	0.32432030219594743	2.7371981701971015	4.758340863299063	1.803432408441833	0.32588556896398846	2.3867597917354444	8.452492382202099	2.847641245801038	2.357101270971357	8.653518686058149	0.32584897222772385	1.3089365679589624	4.698553876958179	0.3258489752839863	1.325848448238901	3.450848972227724	1.3258489719577307	2.325848972227724	2.4508480412726805	2.3414647172223693	1.3258489722277238	0.38068810644403295	2.863176992655458	3.590253269102724	2.8363866072893518	1.3258489722266926	0.3102242703007216	0.07126103805056193	1.5134118441050786	3.3348666795058737	0.31307109269960814	2.3757316724891266	5.373510460825523	1.3375677222277238	0.3258491154976331	5.8563159554086415	2.327435886290224	1.325848972162331	0.82584896962388	3.107302545263983	1.0757269267428422	9.312576807570844	0.32584897222772075	7.075893482461403	0.016125063796186324	1.315066285794553	0.32951109187103683	0.1726176772233692	1.3258486919411399	1.82598797025934	8.409056435311664	0.2663172835065809;...
22.472181011082554	17.312578493352326	6.638086521922955	5.748307774338563	2.266457825225162	18.81467012785494	2.206572471695661	13.682242644439036	16.814653545694025	17.81465903346146	12.527992198053415	18.814653545694025	12.623537823440309	15.729672408432364	17.816034465474598	19.81467412161959	22.814653545694025	17.81758323169993	18.814653545694025	17.814650117947057	20.814648448849255	18.814653545694025	18.85391482718436	8.914599079344129	22.666636842377564	17.782915264444025	19.8146535456945	20.31734046718772	19.69121043981574	24.048318283531067	16.845403215650318	17.78377730138846	19.408464627862845	20.307099438438957	16.30464704262647	17.814409582689144	23.54389251784447	18.814653545694025	19.814653546893116	18.815889514304896	17.980982531659357	17.814653472823352	5.666866890735455	17.81465354569407	6.516495148298802	9.770177219223832	14.298018736101678	20.816851080247236	14.597949116701836	19.84199691445766	17.814653545694025	17.299762499019945	18.826711728279765;...
51.90257696694407	20.122428037152865	11.164135791125396	8.362606781184036	20.407339341904297	33.253460824437624	24.249034834466862	9.535610716011993	48.96904507606302	44.49048228817597	44.08116603143353	26.378921489402693	25.149879032694972	41.277749543803424	42.638713867585075	38.74523332451936	25.11158381013032	37.31709196357019	41.50682293222509	45.500690145663015	34.27649558236843	49.15407472087337	28.087250496348496	22.249072464232817	23.11939631013032	33.98680390707431	37.32283084617314	38.44468403579112	41.1361664894414	52.52432897235279	47.788476780864016	35.11474538906305	44.01765208334889	48.483343898318786	20.14283381013032	19.119396567995697	25.210900125699748	20.11939631013032	40.52550986635032	39.23832774083607	40.25792125587822	32.57062420075704	35.56562961226804	18.119396310130295	16.927672853602807	28.048634737412804	23.499094282276133	41.72003520154378	30.954560937707374	29.55074705680309	20.11939631013032	17.74619196453377	32.51517083144021;...
58.30430707703795	34.20164735616873	26.020986445039952	30.9881044484178	14.973993730593111	40.889792566807365	27.866225443607156	39.11276990559619	49.044349614779435	41.20311280131707	64.50739841969052	33.203134423210564	19.955266429442556	47.698523048906424	52.237903041755935	45.20308711152945	48.841697961037774	53.930347731988576	46.09899826088369	43.20313245194959	40.35148262669646	46.82634342253949	55.27026892869371	37.55622219930927	54.072140429975896	55.55849605511651	38.20313442322013	45.94621070586433	66.80753571370222	46.2234690133822	61.74124841524301	57.872202521171346	49.073707528894886	62.151641374788255	36.203134423210564	33.20313479276058	43.50891747469531	33.703134423210564	38.20313441427486	38.19922783379375	66.13150199964664	42.870154845919934	23.294609526967584	53.625818976081426	15.376910236272538	47.278444404825485	34.821397419985935	38.203485565946046	38.16691089668368	37.437622049768365	33.203134423210564	39.942376003524565	33.51401683134413;...
54.08373492923317	32.798257806026875	33.95087943606947	34.06809616762181	28.64278806439104	50.58395508482959	22.35776690488648	29.450091522030945	51.92662330915284	52.230519902507744	69.02615752647039	44.14998633814755	19.596106171621763	61.097545275094994	52.19944739268946	39.79901331436349	52.56124011641032	49.148960604401694	58.97053958300241	44.05004847941088	49.19268713960143	44.798972919526534	59.77312887595662	37.0360043357699	45.963890396096566	50.88461489181215	54.17867863681169	42.29897361206178	65.87031128359031	46.85771554883104	60.66884104495723	52.42198510206184	45.822808628521315	56.89304278545685	32.814597919526534	33.29515994663236	46.25604601705202	32.798972919526534	37.79897292405112	49.869785554286146	58.84845389849032	52.15968239820038	23.36461702555383	48.18784397862178	28.97131068345025	42.702802384020295	35.891984633380055	43.42418384913963	39.88761248029097	47.11473609641013	38.045043056258926	51.96041247935146	46.141997578498405;...
0.40088481839913526	0.7847172163097644	0.6843829935977722	0.6601804709713119	0.8613596827230859	0.9999785499452274	0.5932056334612723	0.6524277843926491	0.7837667500227848	0.8931197921088871	0.9998464085848418	0.8423605000227848	0.7128963105561301	0.9810877188630631	0.9717712230958685	0.7837667287499421	0.9654887920949269	0.7837667500683358	0.7837667500227848	0.7838118463297662	0.8326779829734435	0.8462667500227848	0.783766750022755	0.9797867665493438	0.7837667500227848	0.7837667500227848	0.8915335813807078	0.7837643304520796	0.9988160888167089	0.6814824548091383	0.652186357351644	0.6363762046120591	0.816845954569191	0.9435708449799385	0.5576492706988962	0.4543320640778322	0.40204761285125495	0.7214635122673458	0.9856319247633168	0.7858263584309177	0.88700133026755	0.8800807876810915	0.937312195226425	0.7857198750227845	0.4865095748306403	0.999963087898482	0.6731815771045399	0.785719743508168	0.9121341137886515	0.783766710459388	0.7837667500227848	0.9992584740162559	0.999087083722307;...
0.966674994308654	0.9650411017539193	0.7046105830406719	0.9370805072996553	0.9370805073650564	0.9650369657860582	0.9800975588628543	0.47935280580176937	0.9650392415566812	0.965038541605365	0.9999094983436786	0.9650392415566812	0.8066264757529964	0.984840047418555	0.9668817404422126	0.9650398059585864	0.9650392415566812	0.965039241456082	0.9650392415566812	0.965039205080018	0.9650391305262953	0.9650392415566812	0.9650392415566812	0.9370837951915788	0.9650392415566812	0.9650392415566812	0.9519847107054251	0.9650391099412065	0.9980532102198286	0.6751225626698045	0.8124189126878401	0.7805410433823297	0.6138893299616942	0.9721440814378841	0.4767579915566812	0.52914364132359	0.5776564395005958	0.9650392415566812	0.9650392412344923	0.9664430519731553	0.8811527466061156	0.9650392415749075	0.9879320372252417	0.9650392415566891	0.49918285537991003	0.9999814604564892	0.5435095363354792	0.965039232665146	0.4217860291806856	0.9650392506838394	0.4964512581545834	0.9982229183227128	0.8536871768250345;...
0.9740499056761404	0.48062673195259886	0.7369711982294874	0.4050795760692105	0.9118255791744145	0.7462207303650679	0.997905382945969	0.41797830415179743	0.9853404886140084	0.6961557008848572	0.9998833979312102	0.5753589017350644	0.8112300429958045	0.9318551211091535	0.8749999164215948	0.6952024427837669	0.6961746324893574	0.696174623706799	0.7142902893015167	0.9958283762683359	0.9655496353268858	0.920306002187712	0.696174632489505	0.9318468876769724	0.9291368787980127	0.9461746324893574	0.8957321838415419	0.9384126143273346	0.9994106152849603	0.8906798289149146	0.9919590051814201	0.9840235222564497	0.6795792644246044	0.9854932122540956	0.6961746324893574	0.5406701856548034	0.8426897683705645	0.47989274283322336	0.9734773351534389	0.6961744895563734	0.9883880742841034	0.6961745622301775	0.953671954701032	0.963055482261194	0.706633404108062	0.9997711181671044	0.7715834936992843	0.9536094715544027	0.6335431058928185	0.9881986026090699	0.5211837878693717	0.5745260251988032	0.6974149929147638;...
0.9920546041066133	0.6403671993670956	0.5015418621204234	0.5170991851126415	0.9465474979139664	0.7845392695933369	0.8144519851268289	0.4000152587890625	0.7689329411621377	0.9727805278423551	0.9997582571263477	0.6403037582308928	0.6941934270625469	0.7770867391809231	0.6866171588734238	0.9900083361231276	0.8348324480627392	0.6960094836776212	0.6403037582308928	0.8099346748989855	0.9268542334789035	0.8903037582308928	0.6442100082304508	0.5404616495427896	0.9440173434228111	0.9701478325202411	0.9335078663497498	0.7398120333048809	0.998118176743241	0.9292716614714346	0.9588176042864929	0.7441854212866852	0.9674798540021489	0.8932595709839778	0.5361056980789805	0.6403034100488512	0.934308899051914	0.6403037582308928	0.8479057955465835	0.87866220810611	0.8069934261815652	0.8063208911023365	0.9504126703270382	0.6403037582308895	0.5066307899031983	1.0	0.41750060819123397	0.6798161732978671	0.40311398960807493	0.9157847446054607	0.6403294925362746	0.4014876867394041	0.6465573637279493;...
0.9928886953836925	0.9792229831433253	0.8567816645834364	0.8925181343648972	0.7586793099793234	0.47922681268773315	0.756350409373052	0.7920935473906224	0.4792229564172124	0.9792229315629707	0.9999517763916196	0.4818322093468999	0.6960373345056036	0.9486648317000157	0.4792229563056761	0.9792229563870831	0.9914392874109299	0.4792229566643972	0.4792229564172124	0.889379105992736	0.47922222214332155	0.9792229564172124	0.47922295641722973	0.9644278873663527	0.9572247753078602	0.4792229564172124	0.9039352886800426	0.9726984968619357	0.9939803443191189	0.9807008187603684	0.9795654695867414	0.7546668531247764	0.9993798430576166	0.9688681891636535	0.9987542064172124	0.47922293232130586	0.9996387189741439	0.4792229564172124	0.479222956778282	0.6491603397371017	0.989671757931972	0.7449236670272005	0.7645541734417021	0.7493227820020922	0.4726188074661589	0.6892164569212522	0.49346499656188764	0.8222674180436702	0.611889127121886	0.47922292687211115	0.4792229564172124	0.9997029017888067	0.9121809150156703;...
64.84160163971507	54.70200417842421	16.862318222098043	39.60926390230574	8.966045934068038	58.202917342701596	17.7334299788167	37.40083777655211	55.20290626149756	58.20291086808604	70.69813978068673	55.2055063949166	45.10969679173504	44.85722658020757	59.09056240634012	61.944644276186835	61.23415626149756	56.403211249491285	59.22829688649756	58.20290603180887	57.20289644742026	59.26150001149756	55.20309419494272	28.689700512691058	64.65515703541259	56.20290626149756	59.20290626149693	56.202907034912194	67.35242721502651	59.41876948904724	59.4587060194255	55.296023402860556	60.030270045985176	66.55293135177239	49.96403397739792	54.20290625743094	59.04566977235487	54.70290626149756	54.202960357362734	56.807839512844765	63.93670190411977	55.20290626167012	52.65077548375323	56.06294667586458	12.303047972642172	40.89660336355452	41.03909784921581	55.20429467455379	43.00179666035485	56.171656663562466	53.20290626149756	38.10485933152	55.28367365356584;...
78.75719074873983	62.76257707618675	46.13738831097915	37.64012235382129	28.872980299540053	78.3054202187485	38.689810334677524	31.98227285495313	62.2622084872238	68.76220625124425	75.9709683667148	63.7622084872238	24.594753350442804	78.34847575176052	62.762208480045004	66.7621989983101	78.70181498730796	62.762208481488635	74.82802961207386	65.76220813376565	76.9057185574814	72.8208022372238	63.762208487225074	47.48492883840865	67.2622084872238	62.7622084872238	63.77319481535295	72.60682971663593	79.1886173476507	73.67661579759249	76.3343442161265	80.39169055108317	69.0533856526884	77.32768719790892	63.7694220797531	62.762229225581656	81.37098603984975	61.2622084872238	79.95589522676566	65.76291027126182	77.72027207057691	77.95433265603559	18.58995211685283	65.15016925595467	18.336933186305792	35.122839907401016	63.81324356808722	66.7622082001694	42.9955353078672	77.9707039104868	63.0122084872238	36.921948935888615	86.11506474212287;...
90.63940693149004	75.76879194771003	35.309232315041456	30.698611609832522	18.60500342825406	79.51709157215127	25.13614289834661	55.214925820393994	68.08118245151417	79.5170693419977	83.5993331337388	73.51730693809108	30.063917065407082	58.93310316256527	76.51706279186128	77.5170627995727	91.54511375676445	82.81070755919346	82.32527406230946	87.52937291464036	78.01704195286003	88.64449697274851	75.5287815474676	38.62713947091339	81.51706279746608	76.51706279746608	77.52048076621463	79.51706553365038	87.39437198351155	90.00867796230204	88.89456597349805	78.65418825461086	79.58831304979968	93.73957636866668	74.97910189434818	76.51707833588598	86.43403992494	75.51706279746608	86.21508742771269	79.51723057276354	64.87635824475092	76.51706277674374	49.773866484206856	76.02074445928108	51.84013135493913	68.67704038857427	73.66858483517724	79.51706257144473	45.48902680650367	78.51706194072285	75.51749004355983	55.96986954412272	81.59055343903088;...
87.64894515511924	69.85815042986864	39.70416919625959	43.27072867164678	25.791991168657255	91.63273857058812	53.46599639871748	44.047278035348114	85.03546720034903	79.69845114334294	89.16277928610043	64.09593589269782	51.56267396110732	68.46178413808902	98.43454149497572	84.86870941723397	77.85796050840257	90.3628557070452	89.72229009661758	90.19089555952227	94.87444114760018	87.68597984326587	73.90226832170461	47.18333323615762	74.85796050840257	69.85796050840257	91.5805636674368	76.85796323359425	86.65198270966921	94.96622209439877	89.52484768077107	74.14796789165194	87.81929640166112	91.12781872896116	69.88530425840257	69.85796097402175	92.23863326709451	69.85796050840257	94.7859227484457	90.83323450571496	70.52392718016607	72.85796050640279	23.785893453846175	88.33864145896527	33.978514540196684	70.37136201506891	42.57511819018768	82.41517754158326	59.30611437989076	86.67451591095491	81.21256919849687	38.160913369475935	89.38881382600373;...
98.5947513930343	76.89104505103779	38.21832956466325	43.965486238443326	26.455569721221075	90.75888015525469	49.59787917263711	59.27741906478796	90.2608190727917	93.26085674508178	96.4523597284124	65.94105916455881	50.238124185368626	76.258070826078	91.26081909739518	92.26081907884905	94.7608190727917	91.26081907178009	93.38893186576045	92.26081920881079	92.26081297148963	94.2842565727917	90.2686315727921	56.531326986684554	95.2608190727917	91.2608190727917	92.2608190727913	93.26082088020937	96.60479177958659	96.13934490729385	97.44340522607767	91.38947745803536	99.69163100315706	96.55141206239341	75.96102623079003	90.76081899799345	99.68810252240341	85.54790052020411	92.26069700646012	91.26451166605318	78.03781714585348	93.26081907277107	35.05565546872951	91.2608190727917	51.187654575434294	64.89952679069974	65.12365396167144	91.26081894906802	64.72174260871745	91.26081871562755	90.2608190727917	58.69522844356848	98.32971523610225;...
0.0	0.0	0.5815044606075295	0.44022094494829855	0.6044628948967983	0.0	0.7638397760330673	0.5657934788645814	0.0	0.0	0.0	0.0	0.7799120788457772	0.985455323710195	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.9510949499754489	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.7682888777391816	0.0	0.4653326034441252	0.0	0.7893973720828832	0.8025187977160486	0.41817945481747737	0.0	0.9074242036041307	0.0	0.0	0.6914329997664286	0.0;...
0.0	0.0	0.49640247394828985	0.4328168319102153	0.5317267814314268	0.0	0.8074690846666429	0.829214111475369	0.0	0.0	0.0	0.0	0.798629267031362	0.8021979074129104	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.9483441745315234	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.9326778607586724	0.0	0.5159941901833737	0.0	0.886989989973124	0.9081894840791611	0.5891020655160554	0.0	0.6098612187333146	0.0	0.0	0.5699148099980222	0.0;...
0.0	0.0	0.556898443914759	0.557440768554419	0.7579057737232748	0.0	0.8966468009488173	1.0	0.0	0.0	0.0	0.0	0.5236536897415314	0.678176979768926	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.9211679637549681	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.996895059002205	0.0	0.4588815483166649	0.0	0.7840956407180928	0.772885604387448	0.4269351137980568	0.0	0.7944634411001499	0.0	0.0	0.6250459833162858	0.0;...
0.0	0.0	0.43780265795102485	0.4038991015804802	0.6618761584059021	0.0	0.948833069562776	0.5707271674464229	0.0	0.0	0.0	0.0	0.5609218681913851	0.938165395862527	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.7018571479214883	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.9475156890581817	0.0	0.667327837170729	0.0	0.8814475711018122	0.9413663425962853	0.5956777221882347	0.0	0.8045826012546872	0.0	0.0	0.6094935432277647	0.0;...
0.0	0.0	0.684457869491361	0.5681682911231181	0.9000057472659275	0.0	0.9177817460107172	0.41360260394922915	0.0	0.0	0.0	0.0	0.40570778288089326	0.9655238272051597	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.9411105101960724	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.9589600194872068	0.0	0.8905998061890874	0.0	0.6423259169449117	0.9847251314522988	0.6319412723049007	0.0	0.7970534569917503	0.0	0.0	0.7655889524837929	0.0;...
0.0	0.0	53.69874513545629	54.797667558316206	50.20243126768813	0.0	58.02715938237844	45.39777469137701	0.0	0.0	0.0	0.0	39.26434879380165	63.601341347789116	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	60.84509464241868	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	81.95125571712617	0.0	53.69019855122936	0.0	73.69511937821851	66.35629014598047	84.78717137105501	0.0	79.75361768377708	0.0	0.0	52.153684456148305	0.0;...
0.0	0.0	54.8378356485984	50.94378537949284	52.5157042117117	0.0	68.14752957784918	52.5809995582279	0.0	0.0	0.0	0.0	55.82108068649572	71.20277208713911	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	55.53845348740266	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	100.09183567392371	0.0	42.36011686046707	0.0	61.298749347597315	70.53561269706108	99.25532431583515	0.0	78.31571603832488	0.0	0.0	81.59479420271748	0.0;...
0.0	0.0	54.94304021292862	77.55916568634457	84.21909464642098	0.0	75.45621668198028	59.232713961511166	0.0	0.0	0.0	0.0	62.13069445687329	84.26647557205393	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	62.51174739463667	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	85.50382830162084	0.0	43.57390429177145	0.0	80.33445484796742	68.97592238493432	82.08647066706529	0.0	93.04225178051396	0.0	0.0	93.73383036033987	0.0;...
0.0	0.0	87.1876580782111	85.89865471343944	76.5959104479488	0.0	60.60542378482036	76.42685282365197	0.0	0.0	0.0	0.0	70.30050509508308	83.06191370009489	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	69.68120996857373	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	91.0669712395448	0.0	87.63795369833892	0.0	69.68121004682354	65.76707464622241	85.00159646945322	0.0	75.25693805485702	0.0	0.0	100.78776322362415	0.0;...
0.0	0.0	85.77904538478508	87.03856320896531	73.80438837020688	0.0	68.57077836182191	75.80106678085497	0.0	0.0	0.0	0.0	85.0463757089571	96.85296000636563	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	70.40570391042567	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	89.48251169066603	0.0	87.04672861958485	0.0	85.0463757089571	68.15351725783142	103.17558659229651	0.0	89.46691197908268	0.0	0.0	91.91097917511448	0.0];



hMax = [];
x0Caps = [];
errorMean = [];

for i = 1:size(launcher,1)
    for j = 1:launcher(i,1)
        configuration.stage{j}.engine = mission.engines{launcher(i,1+j)};
    end

    [mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher(i,:));
    
    if launcher(i,1) == 2
        thrustDataVec = reshape(thrustData(1:20,i),settings.nOptPointsTraj,2,2);
        [timeCollocation, stateCollocation] = totalTrajectoryGlobalGA(launcher(i,:),configuration,mission,settings,thrustDataVec);
    else
        thrustDataVec = reshape(thrustData(1:30,i),settings.nOptPointsTraj,2,3);
        [timeCollocation, stateCollocation] = totalTrajectoryGlobalGA(launcher(i,:),configuration,mission,settings,thrustDataVec);
    end

    if norm(stateCollocation(1:3,end,end)-mission.target.initialPointECI) < 5e5
        hMax = [hMax,max(vecnorm(stateCollocation(1:3,:,end),2,1)-mission.environment.rEarth)];
        x0Caps = [x0Caps,stateCollocation(1:6,1,end)];

        uncertaintyx0 = 0.01.*randn(6,100);
        x0CapsUnc = x0Caps(:,end) + uncertaintyx0.*x0Caps(:,end);
        for k = 1:100
            [tt,xx] = ballisticTrajectory(x0CapsUnc(:,k),mission,[1;0;0],0,100);
            error(k) = norm(xx(1:3,end)-stateCollocation(1:3,end,end));
        end
        errorMean = [errorMean mean(error)];
        

    end
    x0CapsUnc = [];
    
    
end
%%
hVec = linspace(2e5,1e7,1000);
coeff = polyfit(hMax,1./errorMean,2);

p = polyval(coeff,hVec);

regression = figure(1);
plot(hMax,errorMean,'Color',settings.color.orange,'Marker','x','LineStyle','none');
hold on
plot(hVec,1./p,'Color',settings.color.blu);
coeff(1)=0;
p = polyval(coeff,hVec);
plot(hVec,1./p,'Color',settings.color.terracotta)
xlabel('$h_{max}$ [m]')
ylabel('Mean Error')
legend('points','quadratic regression','linear regression','Location','northeast','Orientation','vertical')
setPlotSettings(title('Mean error vs Max altitude'))
exportStandardizedFigure(regression,'regression',0.55,1.5,'ChangeColors',false,'AddMarkers',false,'overwriteFigure',true,'exportFIG',true,'exportPDF',false,'figurePath','images')