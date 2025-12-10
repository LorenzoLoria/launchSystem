% -------------------------------------------------------------------------
% ---------------------- MONTECARLO ANALYSIS ------------------------------
% -------------------------------------------------------------------------

% Initialization

clear all
clc
close all

addpath(genpath("..\..\"))

[mission,settings] = dataStructGlobal;

launcher = [2,1,4,4,0.4056,0.4016,0.7];

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

[timeCollocation, stateCollocation] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustData);
figure
EarthPlot(mission.environment.rEarth)
hold on
plot3(stateCollocation(1,:,1),stateCollocation(2,:,1),stateCollocation(3,:,1))
plot3(stateCollocation(1,:,2),stateCollocation(2,:,2),stateCollocation(3,:,2))
plot3(stateCollocation(1,:,3),stateCollocation(2,:,3),stateCollocation(3,:,3))
plot3(mission.target.initialPointECI(1),mission.target.initialPointECI(2),mission.target.initialPointECI(3),'ro')

% Number of elements

sizeMC = 10000;

% Initialization of the variables

distanceFromTarget = zeros(1,sizeMC);
cumulativeMean = zeros(length(distanceFromTarget),1);
hVec = 0:100;

% Computation of the Wind Profiles

[meanWind, varWind] = GRAM07_HWM07_annual(hVec);
windUncertainty = sqrt(varWind) .* randn(sizeMC,1);
WindVelocityMag = meanWind + windUncertainty;
windAngVel = WindVelocityMag ./ (mission.environment.rEarth + hVec);
lonInitial = mission.launchBase.lonInitial;
montecarlo.vxWind = - windAngVel .* (mission.environment.rEarth + hVec) .* sin(lonInitial) ;
montecarlo.vyWind = windAngVel .* (mission.environment.rEarth + hVec) .* cos(lonInitial) ;

hold on
parfor parforiter = 1:sizeMC
    
    % Functions for wind profile on ECI (rotated inside the dynamics)
    windVelXFun = griddedInterpolant(hVec,montecarlo.vxWind(parforiter,:),'linear','linear');
    windVelYFun = griddedInterpolant(hVec,montecarlo.vyWind(parforiter,:),'linear','linear');

    % Integration of the trajectory
    [timeCollocation, stateCollocation] = totalTrajectoryMontecarlo(launcher,configuration,mission,settings,windVelXFun,windVelYFun,thrustData,parforiter);
    
    % Error
    distanceFromTarget(parforiter) = norm(stateCollocation(1:3,end,end)-mission.target.initialPointECI);

end

% Computation of the cumulative mean

k = 0;

for j = 1:1:length(distanceFromTarget)
    k = k+1;
    cumulativeMean(k) = mean(distanceFromTarget(1:j));
end

figure
plot(cumulativeMean)

%%
% % %%
% % % structuralMassUncertainty = 0.0005 * randn(sizeMC,1);
% % % propellantMassUncertainty = 0.0005 * randn(sizeMC,1);
% % % specificImpulseUncertainty = 0.00 * randn(sizeMC,1);
% % 
% % % assembly of all possible combination
% % Matr = [];
% % for i = 1 :  optimisation.nStages
% %     structuralMass(:,i) = optimisation.stage{i}.structuralMass + structuralMassUncertainty .* optimisation.stage{i}.structuralMass;
% %     propellantMass(:,i) = optimisation.stage{i}.mProp + propellantMassUncertainty .* optimisation.stage{i}.mProp;
% %     specificImpulse(:,i) = optimisation.stage{i}.Isp + specificImpulseUncertainty .* optimisation.stage{i}.Isp;
% %     [T, A, N] = ndgrid(structuralMass(:,i), propellantMass(:,i), specificImpulse(:,i));
% %     Matr =[Matr ,T(:),A(:),N(:)];
% %     optimisation.stage{i}.structuralMass = 0;
% %     optimisation.stage{i}.mProp = 0;
% %     optimisation.stage{i}.Isp = 0;
% % end
% % 
% % % Shuffle of the Matrix elements
% % for i=1:size(Matr,2)
% %     shuffledMatrix(:,i) = Matr(randperm(size(Matr,1)),i);
% % end
% % 
% % % Start of Montecarlo analysis
% % 
% % 
% % nStages = optimisation.nStages;
% % 
% % Ms = zeros(1,nStages);
% % Mp = zeros(1,nStages);
% % Isp = zeros(1,nStages);
% % 
% % for stage = 1:nStages
% %     optimisation.montecarlo.massStruct(:,stage) = shuffledMatrix(:,1 + (stage-1)*3);
% %     optimisation.montecarlo.massProp(:,stage) = shuffledMatrix(:,2 + (stage-1)*3);
% %     optimisation.montecarlo.impSpec(:,stage) = shuffledMatrix(:,3 + (stage-1)*3);
% %     optimisation.montecarlo.StageMass(:,stage) = optimisation.montecarlo.massStruct(:,stage) + optimisation.montecarlo.massProp(:,stage);
% % end
% % 
% % optimisation.montecarlo.initialTime = timeCollocation(1,1:end);
% % optimisation.montecarlo.finalTime = timeCollocation(end,1:end-1);
% % 
% % distanceFromTarget = zeros(sizeMC^NvarsUnc,1);
% % finalPos = zeros(3,size(shuffledMatrix,1));
% % 
% % parfor parforiter = 1:size(shuffledMatrix,1)
% % 
% %     [timeCollocation, stateCollocation] = totalTrajectoryMontecarlo(mission,optimisation,thrustData,parforiter);
% %     latInitial = mission.target.latInitial ;
% %     latFinal = latInitial ;
% %     lonInitial = mission.target.lonInitial ;
% %     omega = mission.target.omega ;
% %     lonFinal = lonInitial + omega * timeCollocation(end,end) ;
% %     targetFinalPos = 6371000*[cos(latFinal)*cos(lonFinal); cos(latFinal)*sin(lonFinal); sin(latFinal) ];
% % 
% %     distanceFromTarget(parforiter) = norm(stateCollocation(1:3,end,end) - targetFinalPos);
% %     finalPos(:,parforiter) = stateCollocation(1:3,end,end);
% % 
% % end
% % 
% % k = 0;
% % for j = 1:100:length(distanceFromTarget)
% %     k = k+1;
% %     cumulativeMean(k) = mean(distanceFromTarget(1:j));
% % end
% % 
% % plot(cumulativeMean)
% % %%
% % figure
% % EarthPlot(mission.environment.rEarth)
% % hold on
% % %scatter3(targetFinalPos(1),targetFinalPos(2),targetFinalPos(3),'rx')
% % scatter3(finalPos(1,:),finalPos(2,:),finalPos(3,:),'bx')
% % 
% % 
% % %% PROVA PORCODIO 2 VISTO CHE LA PRIMA FA CAGARE
% % clear all
% % clc
% % close all
% % addpath(genpath("..\..\"))
% % 
% % sizeMC = 30;
% % NvarsUnc = 3;
% % 
% % [mission,optimisation] = dataStruct;
% % thrustData(:,:,1) =     [0.9956   -2.7415   13.0559; ...
% %     0.9243   -2.3127    4.7079; ...
% %     0.9720   -2.4832   27.0284; ...
% %     0.9607   13.3110   22.9276; ...
% %     0.9868  -33.8890   83.4509];
% % thrustData(:,:,2) =    [0.7495   26.9887   29.4877; ...
% %     0.7732  -22.4865   93.7944; ...
% %     0.8558   42.2771   77.4467; ...
% %     0.8221   30.4242   77.8149; ...
% %     0.4019   35.0922   87.7839];
% % 
% % 
% % 
% % [timeCollocation, stateCollocation] = totalTrajectory(mission,optimisation,thrustData);
% % 
% % % nominalCapsuleMass = mission.capsule.weigth;
% % % capsuleMassUncertainty = ( 8.8 * 6 + 10 * 6 ) / mission.capsule.weigth + 0.1;
% % % stdDeviationMass = capsuleMassUncertainty * randn(30,1);
% % nominalCapsulePosition = stateCollocation(1:3,end,end-1);
% % capsulePositionUncertainty = max(norm(nominalCapsulePosition) * 1e-6,1e-2);
% % stdDeviationPosition = capsulePositionUncertainty * randn(30,3);
% % nominalCapsuleVelocity = stateCollocation(4:6,end,end-1);
% % capsuleVelocityUncertainty = max(norm(nominalCapsuleVelocity) * 1e-6,1e-2);
% % stdDeviationVelocity = capsuleVelocityUncertainty * randn(30,3);
% % 
% % capsulePosition = nominalCapsulePosition' + stdDeviationPosition;
% % %capsuleVelocity = nominalCapsuleVelocity' + stdDeviationVelocity; % columns are vx,vy,vz
% % 
% % [A,B,C] = ndgrid(capsulePosition(:,1),capsulePosition(:,2),capsulePosition(:,3));
% % Matr =[A(:),B(:),C(:)];
% % 
% % % Shuffle of the Matrix elements
% % for i=1:size(Matr,2)
% %     shuffledMatrix(:,i) = Matr(randperm(size(Matr,1)),i);
% % end
% % 
% % distanceFromTarget = zeros(sizeMC^NvarsUnc,1);
% % nDeval = 100;
% % windDirection = [1 0 0]';
% % parfor parforiter = 1:size(shuffledMatrix,1)
% % 
% %     x0 = [shuffledMatrix(parforiter,:)';stateCollocation(4:6,end,end-1)];
% %     [tt,xx] = ballisticTrajectory(x0,mission,windDirection,timeCollocation(end,end-1),nDeval);
% %     latInitial = mission.target.latInitial ;
% %     latFinal = latInitial ;
% %     lonInitial = mission.target.lonInitial ;
% %     omega = mission.target.omega ;
% %     lonFinal = lonInitial + omega * tt(end) ;
% %     targetFinalPos = 6371000*[cos(latFinal)*cos(lonFinal); cos(latFinal)*sin(lonFinal); sin(latFinal) ];
% %     distanceFromTarget(parforiter) = norm(xx(1:3,end)-targetFinalPos);
% % 
% % end
% % 
% % k = 0;
% % for j = 1:100:length(distanceFromTarget)
% %     k = k+1;
% %     cumulativeMean(k) = mean(distanceFromTarget(1:j));
% % end
% % 
% % plot(cumulativeMean)
