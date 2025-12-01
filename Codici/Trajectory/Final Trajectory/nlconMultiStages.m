function[cin,ceq] = nlconMultiStages(x,mission,opt,option2D)

thrustDataVec = x;

[timeCollocation,stateCollocation] = totalTrajectory(mission,opt,thrustDataVec,option2D);

latInitial = mission.target.latInitial ;
latFinal = latInitial ;
lonInitial = mission.target.lonInitial ;
if option2D == 1
    omega = 0;
else
omega = mission.target.omega ;
end
lonFinal = lonInitial + omega * timeCollocation(end,end) ;
targetFinalPos = 6371000*[cos(latFinal)*cos(lonFinal); cos(latFinal)*sin(lonFinal); sin(latFinal) ];


for i = 1:opt.nStages
    time      = linspace(timeCollocation(1,i),timeCollocation(end,i),length(timeCollocation)-1);
    v         = stateCollocation(4:6,:,i);
    vNorm     = sqrt( stateCollocation(4,:,i).^2 + stateCollocation(5,:,i).^2 +stateCollocation(6,:,i).^2 );
    acc       = diff(vNorm)./(time(2)-time(1));
    accMaxStage(i) = max(acc);
end

accMax = max(accMaxStage);
%cin = [(accMax-6*mission.environment.g0)/mission.environment.g0 ; (norm(stateCollocation(1:3,end,end) - targetFinalPos) - 5e3)/1e3];
%ceq = norm(stateCollocation(1:3,end,end) - targetFinalPos);
cin = [(accMax-6*mission.environment.g0)/mission.environment.g0 ;];
ceq = [ norm(stateCollocation(1:3,end,end) - targetFinalPos) ];
end