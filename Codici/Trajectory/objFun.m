function[objective] = objFun(x,mission)

thrustDataVec = x;
tSpan = [0 mission.launcher.engines{1}.isp];
tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));
%thrustData = @(t) [interp1(tVec,thrustDataVec(:,1),t);interp1(tVec,thrustDataVec(:,2),t);interp1(tVec,thrustDataVec(:,3),t)];

F_thrust = griddedInterpolant(tVec, thrustDataVec, 'linear', 'none');  
thrustData = @(t) F_thrust(t).';  

x0 = [mission.initialPoint'; 0; 0; 0; mission.launcher.engines{1}.m0];
[ttL,xxL] = launcherTrajectory(x0,mission, thrustData,0,thrustDataVec);

x0C = xxL(1:6,end);

windDirection = -1;
windDirection = windDirection/norm(windDirection);

[ttC,xxC] = ballisticTrajectory(x0C,mission,windDirection,ttL(end));

objective = ttC(end) ;

end