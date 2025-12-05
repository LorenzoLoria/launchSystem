function[objective] = objFunFMCTraj(x,launcher,opt,mission,settings)

thrustDataVec = x;
%thrustDataVec is the variable containing informations about the
%optimisation paramets. As of right now, the trajectory optimisation
%variables are throttoling and angle. The best way to divide them is by
%using a optVarNum x 2 x nStages matrix, since it will help in keeping the
%code organised correctly.

[timeCollocation,stateCollocation] = totalTrajectoryGlobalGA(launcher,opt,mission,settings,thrustDataVec);

% latInitial = mission.target.latInitial ;
% latFinal = latInitial ;
% lonInitial = mission.target.lonInitial ;
% omega = mission.target.omega ;
% lonFinal = lonInitial + omega * timeCollocation(end,end) ;
% targetFinalPos = 6371000*[cos(latFinal)*cos(lonFinal); cos(latFinal)*sin(lonFinal); sin(latFinal) ];

objective = timeCollocation(end,end) ;

objective = objective;
end