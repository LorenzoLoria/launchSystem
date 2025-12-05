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


[mer,~,configuration] = initialMassEstimation(mission,configuration,settings,launcher);


[xGATraj, fvalGATraj] = ga( @(x) objFunGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings), ...
                        launcher(1)*2*settings.nOptPointsTraj,...
                        [],[],[],[],settings.lowerBoundsGA,settings.upperBoundsGA, ...
                        @(x) nlconGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings),settings.gaTrajOptions);


if fvalGATraj > 6000  
    output = 3000000;
return

end

error = 101;

while error > 100

    [~,fvalFMCTraj] = fmincon ( @(x)objFunFMCTraj(x,launcher,configuration,mission),...
        xGATraj,[],[],[],[],...
        settings.lowerBoundsFMC-eps,settings.upperBoundsFMC+eps,...
        @(x) nlconFMCTraj(x,launcher,configuration,mission,settings.trajectoryOption2D),...
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