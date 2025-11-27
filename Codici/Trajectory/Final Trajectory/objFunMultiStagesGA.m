function[objective] = objFunMultiStagesGA(x,mission,opt)

thrustDataVec = x;
%thrustDataVec is the variable containing informations about the
%optimisation paramets. As of right now, the trajectory optimisation
%variables are throttoling and angle. The best way to divide them is by
%using a optVarNum x 2 x nStages matrix, since it will help in keeping the
%code organised correctly.

[timeCollocation,stateCollocation] = totalTrajectory(mission,opt,thrustDataVec);

latInitial = mission.target.latInitial ;
latFinal = latInitial ;
lonInitial = mission.target.lonInitial ;
omega = mission.target.omega ;
lonFinal = lonInitial + omega * timeCollocation(end,end) ;
targetFinalPos = 6371000*[cos(latFinal)*cos(lonFinal); cos(latFinal)*sin(lonFinal); sin(latFinal) ];

objective =   0.8* norm (stateCollocation(1:3,end,end)-targetFinalPos)/1000 + 0.2*timeCollocation(end,end)/3000;

end