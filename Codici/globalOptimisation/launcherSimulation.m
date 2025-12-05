function [output] = launcherSimulation(launcher,mission,opt,settings,option2D,nlconFlag)

persistent previousLauncher 

if isempty(previousLauncher)
    previousLauncher = zeros(length(launcher),1);
end

if launcher([1,2:launcher(1)+1,5:4+launcher(1)]) == previousLauncher([1,2:launcher(1)+1,5:4+launcher(1)])

else


[mer,staging] = initialMassEstimation(mission,opt,settings,launcher);

initialMassGuess= mission.capsule.weigth;
opt.m0Tot = initialMassGuess;

        for i = 1:launcher(1)
            initialMassGuess = initialMassGuess +staging{i}.mStage;
        end    



[xGATraj, fvalGATraj] = ga( @(x) objFunGATraj( reshape(x,settings.TrajOptimisationPoints,2,launcher(1)),launcher, mission,opt,option2D), ...
                        launcher(1)*2*settings.TrajOptimisationPoints,...
                        [],[],[],[],settings.lowerBoundsGA,settings.upperBoundsGA, ...
                        @(x) nlconGATraj( reshape(x,settings.TrajOptimisationPoints,2,launcher(1)), mission,opt,launcher,option2D),settings.gaTrajOptions);


if fvalGATraj > 6000  
    output = 3000000;
return

end

error = 101;

while error > 100

    [~,fvalFMCTraj] = fmincon ( @(x)objFunFMCTraj(x,launcher,mission,opt,1),...
        xGATraj,[],[],[],[],...
        settings.lowerBoundsFMC-eps,settings.upperBoundsFMC+eps,...
        @(x) nlconFMCTraj(x,launcher,mission,opt,1),...
        settings.fminconTrajOptions);
    error = 99;

end





if nlconFlag
    output = [];
else
    totalMass = mission.capsule.weigth;

        for i = 1:launcher(1)
            totalMass = totalMass +staging{i}.mStage;
        end    
    output.launcherMass = totalMass;
    output.tof = fvalFMCTraj; 

end



end