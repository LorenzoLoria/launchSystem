function [output] = launcherSimulation(launcher,mission,settings,nlconFlag)

persistent previousLauncher 

if isempty(previousLauncher)
    previousLauncher = zeros(length(launcher),1);
end

if launcher([1,2:launcher(1)+1,5:4+launcher(1)]) == previousLauncher([1,2:launcher(1)+1,5:4+launcher(1)])

else

% retrive launcher engines 
for i = 1:launcher(1)
    configuration.stage{i}.engine = mission.engines{launcher(1+i)};
end


[mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);


if launcher(1) == 1

elseif launcher(1) == 2
    localGAsettings = settings.gaTrajOptions ;
    localGAsettings.InitialPopulationMatrix = settings.initialPopulationGATraj.twoStages' ;
else

end


[xGATraj, fvalGATraj] = ga( @(x) objFunGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings), ...
                        launcher(1)*2*settings.nOptPointsTraj,...
                        [],[],[],[],settings.lowerBoundsGA,settings.upperBoundsGA, ...
                        @(x) nlconGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings),localGAsettings);


if fvalGATraj > 6000

    if nlconFlag
        output = 1e9;
    else
        output.launcheMass = 1e9 ;
        output.tof = 1e9 ;
    end

    return

end



error = 101;
xGATrajRS = reshape(xGATraj,settings.nOptPointsTraj,2,launcher(1));
while error > 100

    [~,fvalFMCTraj] = fmincon ( @(x)objFunFMCTraj(x,launcher,configuration,mission,settings),...
        xGATrajRS,[],[],[],[],...
        settings.lowerBoundsFMC-eps,settings.upperBoundsFMC+eps,...
        @(x) nlconFMCTraj(x,launcher,configuration,mission,settings),...
        settings.fminconTrajOptions);
    error = 99;

end

if nlconFlag
    output = [];
else
    
    output.launcherMass = configuration.totalMass;
    output.tof = fvalFMCTraj; 

end

end