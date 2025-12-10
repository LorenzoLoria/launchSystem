function [output] = launcherSimulation(launcher,mission,settings,nlconFlag)

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


    [~,~,configuration] = initialMassEstimation(mission,configuration,settings,launcher);


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
    if fvalGATraj > 6000

        if nlconFlag
            output = 1e9;
        else
            output.launcherMass = 1e9 ;
            output.tof = 1e9 ;
        end

        return

    end



    error = 101;
    xGATrajRS = reshape(xGATraj,settings.nOptPointsTraj,2,launcher(1));

    localLowerBoundsFMC = settings.lowerBoundsFMC(:,:,1:launcher(1)) ;
    localUpperBoundsFMC = settings.upperBoundsFMC(:,:,1:launcher(1)) ;

    while error > 100

        [~,fvalFMCTraj] = fmincon ( @(x)objFunFMCTraj(x,launcher,configuration,mission,settings),...
            xGATrajRS,[],[],[],[],...
            localLowerBoundsFMC-eps,localUpperBoundsFMC+eps,...
            @(x) nlconFMCTraj(x,launcher,configuration,mission,settings),...
            settings.fminconTrajOptions);
        error = 99;

    end

    totalMass = configuration.totalMass ;
    tof = fvalFMCTraj ;

    previousLauncherMass = totalMass ;
    previousTof = tof ;

end

if nlconFlag
    output = [];
else

    output.launcherMass = totalMass;
    output.tof = tof;

end

end