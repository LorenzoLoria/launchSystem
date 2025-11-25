function[objective] = objFun(x,mission)

% thrustDataVec is a 5x3x2 matrix 

thrustDataVec = x;
tSpan = [0 3*24*3600];
tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));

% 1st stage
F_thrust1 = griddedInterpolant(tVec, thrustDataVec(:,:,1), 'linear', 'none'); 
thrustData1 = @(t) F_thrust1(t).'; 
x0 = [0; 0; mission.environment.rEarth; 0; 0; 0; mission.launcher.engines{1}.m0];
[ttL1,xxL1] = launcherTrajectory(x0,mission, thrustData1,mission.launcher.engines{1}.mPropellant1,0);

% 2nd stage
F_thrust2 = griddedInterpolant(tVec, thrustDataVec(:,:,2), 'linear', 'none'); 
thrustData2 = @(t) F_thrust2(t).';
x0L2 = [xxL1(1:6,end);xxL1(7,end)-mission.launcher.engines{1}.ms1];
[ttL2,xxL2] = launcherTrajectory(x0L2,mission, thrustData2,mission.launcher.engines{1}.mPropellant2,ttL1(end));

windDirection = -1+2*rand(3,1);
windDirection = windDirection/norm(windDirection);

x0C = xxL2(1:6,end);
[ttC,~] = ballisticTrajectory(x0C,mission,windDirection,ttL2(end));

objective = ttC(end) ;

end