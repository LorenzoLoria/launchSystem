function [cin,ceq] = nlconFMCTraj(x,launcher,opt,mission,settings)


option2D = settings.trajectoryOption2D ;
thrustDataVec = x;

[timeCollocation,stateCollocation] = totalTrajectoryGlobalGA(launcher,opt,mission,settings,thrustDataVec);

% latInitial = mission.target.latInitial ;
% latFinal = latInitial ;
% lonInitial = mission.target.lonInitial ;
% if option2D == 1
%     omega = 0;
% 
% else
%     omega = mission.target.omega ;
% 
% end
% lonFinal = lonInitial + omega * timeCollocation(end,end) ;
%targetFinalPos = 6371000*[cos(latFinal)*cos(lonFinal); cos(latFinal)*sin(lonFinal); sin(latFinal) ];

[angleWrtVelMax,accMax] = findAccAngle(launcher,opt,mission,stateCollocation,timeCollocation,thrustDataVec,settings) ; 





%cin = [abs(angleWrtVelMax)-20;(accMax-6*mission.environment.g0)/mission.environment.g0 ; (norm(stateCollocation(1:3,end,end) - mission.target.initialPointECI))/1e3 - 2] ; 
cin = [(accMax-6*mission.environment.g0)/mission.environment.g0 ; (norm(stateCollocation(1:3,end,end) - mission.target.initialPointECI))/1e3 - 2] ; 
ceq = [] ;












end