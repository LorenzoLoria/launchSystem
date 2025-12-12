function [outputOBJ,outputNLC] = launcherSimulation(launcher,mission,settings,nlconFlag)

persistent previousLauncher previousTof previousLauncherMass

if isempty(previousLauncher)
    previousLauncher = zeros(length(launcher),1);
end

if launcher([1,2:launcher(1)+1,5:4+launcher(1)]) == previousLauncher([1,2:launcher(1)+1,5:4+launcher(1)])

    totalMass = previousLauncherMass ;
    tof = previousTof ;
    
else

    % retrive launcher engines
    for i = 1:launcher(1)
        configuration.stage{i}.engine = mission.engines{launcher(1+i)};
    end


    [mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);


    if launcher(1) == 1

        localGAsettings = settings.gaTrajOptions ;
        localGAsettings.InitialPopulationMatrix = settings.initialPopulationGATraj.oneStage' ;

    elseif launcher(1) == 2
        localGAsettings = settings.gaTrajOptions ;
        localGAsettings.InitialPopulationMatrix = settings.initialPopulationGATraj.twoStages' ;
    else
        localGAsettings = settings.gaTrajOptions ;
        localGAsettings.InitialPopulationMatrix = settings.initialPopulationGATraj.threeStages' ;
    end

    localLowerBoundsGA = settings.lowerBoundsGA(1:launcher(1)*2*settings.nOptPointsTraj) ;
    localUpperBoundsGA = settings.upperBoundsGA(1:launcher(1)*2*settings.nOptPointsTraj) ;

    [xGATraj, fvalGATraj] = ga( @(x) objFunGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings), ...
        launcher(1)*2*settings.nOptPointsTraj,...
        [],[],[],[],localLowerBoundsGA,localUpperBoundsGA, ...
        @(x) nlconGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings),localGAsettings);

if isnan(fvalGATraj) || ~isreal(fvalGATraj) || isinf(fvalGATraj)
    keyboard
end
    
if fvalGATraj > 10000
        outputNLC = 100;
        outputOBJ.launcherMass = 1e9 ;           
        outputOBJ.tof = 1e9 ;
previousLauncherMass = 1e9 ;
previousTof = 1e9 ;
        return
end

    
    currentStructuralMass = 0 ;
    for i = 1:launcher(1)
        currentStructuralMass = currentStructuralMass + mer.stage{i}.tankMassFuel + mer.stage{i}.tankMassOx + mer.stage{i}.interStage ;
    end
    

    maxMassErr = 0.05 ; 
    error = maxMassErr + 1 ;

    xGATrajRS = reshape(xGATraj,settings.nOptPointsTraj,2,launcher(1));
    localLowerBoundsFMC = settings.lowerBoundsFMC(:,:,1:launcher(1)) ;
    localUpperBoundsFMC = settings.upperBoundsFMC(:,:,1:launcher(1)) ;

    count = 1 ;
    maxIter = 3 ;
    while error > maxMassErr && count <= maxIter
        
        count = count + 1 ;

        [thrustDataVecFMC,fvalFMCTraj,~,checkViolation] = fmincon ( @(x)objFunFMCTraj(x,launcher,configuration,mission,settings),...
            xGATrajRS,[],[],[],[],...
            localLowerBoundsFMC-eps,localUpperBoundsFMC+eps,...
            @(x) nlconFMCTraj(x,launcher,configuration,mission,settings),...
            settings.fminconTrajOptions);
        
         if isempty(checkViolation.bestfeasible)
            outputNLC = checkViolation.constrviolation;
            outputOBJ.launcherMass = 1e9 ;           
            outputOBJ.tof = 1e9 ;
            previousLauncherMass = 1e9 ;
            previousTof = 1e9 ;
            return
         end

        [timeCollocation,stateCollocation] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustDataVecFMC);
        [maxQData] = externalLoads(timeCollocation, stateCollocation, mission, configuration, launcher, mer, 0) ;
        [internalActions] = loadsFinder(mission, launcher, configuration, maxQData) ; 
        [updatedStructuralMass] = thicknessFunction(mission, launcher, configuration, maxQData, internalActions) ; 

        configuration.totalMass = configuration.totalMass - currentStructuralMass + updatedStructuralMass ;
        error = abs(updatedStructuralMass - currentStructuralMass) / currentStructuralMass ;
        currentStructuralMass = updatedStructuralMass ;   
        
    end

    totalMass = configuration.totalMass ;
    tof = checkViolation.bestfeasible.fval ;
    xConfig = checkViolation.bestfeasible.x;
    
    xConfigSave = zeros(30,1);
    xConfigSave(1:numel(xConfig)) = xConfig(:);
    filename = 'configMAT.mat';
    if isfile(filename)
        result = load(filename,"xConfigSave");
        xConfigSave = [result.xConfigSave,xConfigSave];
    end
        save(filename, 'xConfigSave');
  
    filenameLauncher = 'LauncherMAT.mat';
    if isfile(filenameLauncher)
        result = load(filenameLauncher,"launcher");
        launcherSave = [result.launcherSave;launcher];
    end
    save(filenameLauncher, 'launcherSave');

    previousLauncherMass = totalMass ;
    previousTof = tof ;
    previousLauncher = launcher;
end

    outputNLC = 0;
    outputOBJ.launcherMass = totalMass;
    outputOBJ.tof = tof;

end
