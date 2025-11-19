function[cin, ceq] = nlconGA(x,mission)

thrustDataVec = x;
GM  = mission.environment.GM;
A   = mission.capsule.Area;
Cd  = mission.capsule.supersonicCD;
g0  = mission.environment.g0;
tSpan = [0 3*24*3600];
tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));
%thrustData = @(t) [interp1(tVec,thrustDataVec(:,1),t);interp1(tVec,thrustDataVec(:,2),t);interp1(tVec,thrustDataVec(:,3),t)];
F_thrust = griddedInterpolant(tVec, thrustDataVec, 'linear', 'none');  
thrustData = @(t) F_thrust(t).';  

x0 = [0; 0; mission.environment.rEarth; 0; 0; 0; mission.launcher.engines{1}.m0];

[ttL,xxL] = launcherTrajectory(x0,mission, thrustData);

Tmax = 845e3*4;
cin=[];
for i = 1:length(ttL)
    r(i) = norm(xxL(1:3,i));
    h = r(i) - mission.environment.rEarth;
    rho = mission.environment.gridInterp(h);
    v = xxL(4:6,i);
    T(:,i) = thrustData(ttL(i));
    Tnorm(i) = norm(T(:,i));
    D(:,i) = - 0.5 .* rho .* norm(v) .* A .* Cd .* v;
    G(:,i) = - GM * r(i) /norm(r(i))^3;
    m(i) = xxL(7,i);
    acc(i) = norm((T(:,i) + D(:,i))/m(i) + G(:,i));
end
cin = [(max(Tnorm)-Tmax)/Tmax ;max(acc)-10*g0;(mission.environment.rEarth-min(r))];
x0C = xxL(1:6,end);

windDirection = -1+2*rand(3,1);
windDirection = windDirection/norm(windDirection);

[ttC,xxC] = ballisticTrajectory(x0C,mission,windDirection,ttL(end));

ceq = [];


end