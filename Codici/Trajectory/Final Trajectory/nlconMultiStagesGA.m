function[cin,ceq] = nlconMultiStagesGA(x,mission,opt,option2D)

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

[angleWrtVelMax,accMax] = findAccAngle(mission,opt,stateCollocation,timeCollocation,thrustDataVec,option2D);

%cin = [(accMax-6*mission.environment.g0)/mission.environment.g0 ; (norm(stateCollocation(1:3,end,end) - targetFinalPos) - 5e3)/1e3];
%ceq = norm(stateCollocation(1:3,end,end) - targetFinalPos);
cin = [(accMax-6*mission.environment.g0)/mission.environment.g0]; %; abs(angleWrtVelMax)-20];
ceq = [];
end