function[cin, ceq] = nlcon(x,mission,rtarg)

thrustDataVec = x;
tSpan = [0 3*24*3600];
tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));
%thrustData = @(t) [interp1(tVec,thrustDataVec(:,1),t);interp1(tVec,thrustDataVec(:,2),t);interp1(tVec,thrustDataVec(:,3),t)];
F_thrust = griddedInterpolant(tVec, thrustDataVec, 'linear', 'none');  
thrustData = @(t) F_thrust(t).';  

x0 = [0; 0; mission.environment.rEarth; 0; 0; 0; mission.launcher.engines{1}.m0];

[ttL,xxL] = launcherTrajectory(x0,mission, thrustData);

x0C = xxL(1:6,end);

windDirection = -1+2*rand(3,1);
windDirection = windDirection/norm(windDirection);

[ttC,xxC] = ballisticTrajectory(x0C,mission,windDirection,ttL(end));

Tmax = 845e3*10;
ceq = [norm(xxC(1:3,end) - rtarg)] ;
cin=[];
for i = 1:size(x,1)
    Ti = norm(x(i,:));
    cin(i) = norm(Ti)-Tmax;

end
