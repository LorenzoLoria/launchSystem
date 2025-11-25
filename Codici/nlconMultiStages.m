function[cin,ceq] = nlconMultiStages(x,mission,opt)

thrustDataVec = x;

[timeCollocation,stateCollocation] = totalTrajectory(mission,opt,thrustDataVec);

for i = 1:opt.nStages
    time      = linspace(timeCollocation(1,i),timeCollocation(end,i),length(timeCollocation)-1);
    v         = stateCollocation(4:6,:,i);
    vNorm     = sqrt( stateCollocation(4,:,i).^2 + stateCollocation(5,:,i).^2 +stateCollocation(6,:,i).^2 );
    acc       = diff(vNorm)./(time(2)-time(1));
    accMaxStage(i) = max(acc);
end

accMax = max(accMaxStage);
cin = [accMax-6*mission.environment.g0];
ceq = norm(stateCollocation(1:3,end,end) - mission.target);

end