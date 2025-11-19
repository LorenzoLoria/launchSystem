function[objective] = objFun(x,mission,rtarg)

thrustDataVec = x;
tSpan = [0 3*24*3600];
tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));
%thrustData = @(t) [interp1(tVec,thrustDataVec(:,1),t);interp1(tVec,thrustDataVec(:,2),t);interp1(tVec,thrustDataVec(:,3),t)];

F_thrust = griddedInterpolant(tVec, thrustDataVec, 'linear', 'none');  
thrustData = @(t) F_thrust(t).';  

x0 = [0; 0; mission.environment.rEarth; 0; 0; 0; mission.launcher.engines{1}.m0];
[ttL,xxL] = launcherTrajectory(x0,mission, thrustData,0);

x0C = xxL(1:6,end);

windDirection = -1+2*rand(3,1);
windDirection = windDirection/norm(windDirection);

[ttC,xxC] = ballisticTrajectory(x0C,mission,windDirection,ttL(end));

objective = ttC(end)/1000 ;

end