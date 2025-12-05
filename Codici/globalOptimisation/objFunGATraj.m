function [objective] = objFunGATraj(x,launcher, opt,mission,settings)

thrustDataVec = x;
%thrustDataVec is the variable containing informations about the
%optimisation paramets. As of right now, the trajectory optimisation
%variables are throttoling and angle. The best way to divide them is by
%using a optVarNum x 2 x nStages matrix, since it will help in keeping the
%code organised correctly.

[timeCollocation,stateCollocation] = totalTrajectoryGlobalGA(launcher,opt,mission,settings,thrustDataVec);

latInitial = mission.target.latInitial ;
latFinal = latInitial ;
lonInitial = mission.target.lonInitial ;
if settings.trajectoryOption2D == 1
    omega = 0;
else
omega = mission.target.omega ;
end
lonFinal = lonInitial + omega * timeCollocation(end,end) ;
targetFinalPos = 6371000*[cos(latFinal)*cos(lonFinal); cos(latFinal)*sin(lonFinal); sin(latFinal) ];

objective =   norm(stateCollocation(1:3,end,end)-mission.target.initialPointECI)/1000;

%thetaOscillation = sqrt(sum((diff([x(:,2,1) ; x(:,2,2)])).^2)) ;
%phiOscillation = sqrt(sum((diff([x(:,3,1) ; x(:,3,2)])).^2)) ;

objective = objective; %+ 0.1 * thetaOscillation / 30 + 0.1 * phiOscillation / 60;
end