
clearvars
clc
close all

addpath(genpath("Capsule"))
addpath(genpath("LauncherTrajectory"))
mission = dataStruct;
xx0 = [0*ones(6,1) , 4*845e3*cos(deg2rad(63)) .* ones(6,1), 4*844e3*sin(deg2rad(63)) .* ones(6,1)];
rtarg = [0;mission.environment.rEarth * cos(deg2rad(87.5));mission.environment.rEarth * sin(deg2rad(87.5))];
opts = optimoptions("fmincon","Display","iter","MaxIterations",200,'MaxFunctionEvaluations',10000,'StepTolerance',1e-16,'OptimalityTolerance',1e-8,'FunctionTolerance',1e-19,'ConstraintTolerance',1e-10);

lb = 0*ones(size(xx0,1));
ub = 4*845e3*ones(size(xx0,1),size(xx0,2));

disp('xx0:');
disp(xx0);
disp('LB:');
disp(lb);
disp('UB:');
disp(ub);

[X,FVAL,EXITFLAG,OUTPUT] = fmincon(@(x) objFun(x,mission),xx0,[],[],[],[],-1*ones(size(xx0,1),size(xx0,2)),4*845e3*ones(size(xx0,1),size(xx0,2)),@(x) nlcon(x,mission,rtarg),opts)

%%
x0 = [0; 0; mission.environment.rEarth; 0; 0; 0; mission.launcher.engines{1}.m0];
tSpan = [0 3*24*3600];
thrustDataVec = xx0;
tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));
thrustData = @(t) [interp1(tVec,thrustDataVec(:,1),t);interp1(tVec,thrustDataVec(:,2),t);interp1(tVec,thrustDataVec(:,3),t)];

[ttL,xxL] = launcherTrajectory(x0,mission, thrustData);
plot3(xxL(1,:),xxL(2,:),xxL(3,:),'r')
hold on
x0C = xxL(1:6,end);
windDirection = -1+2*rand(3,1);
windDirection = windDirection/norm(windDirection);
[ttC,xxC] = ballisticTrajectory(x0C,mission,windDirection,ttL(end));
plot3(xxC(1,:),xxC(2,:),xxC(3,:),'r')
plot3(x0(1),x0(2),x0(3),'bo')
plot3(rtarg(1),rtarg(2),rtarg(3),'ro')
err = norm(xxC(1:3,end)-rtarg)
axis equal
%
x0 = [0; 0; mission.environment.rEarth; 0; 0; 0; mission.launcher.engines{1}.m0];
tSpan = [0 3*24*3600];
thrustDataVec = X;
tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));
thrustData = @(t) [interp1(tVec,thrustDataVec(:,1),t);interp1(tVec,thrustDataVec(:,2),t);interp1(tVec,thrustDataVec(:,3),t)];

[ttL,xxL] = launcherTrajectory(x0,mission, thrustData);
plot3(xxL(1,:),xxL(2,:),xxL(3,:),'g')
hold on
x0C = xxL(1:6,end);
windDirection = -1+2*rand(3,1);
windDirection = windDirection/norm(windDirection);

[ttC,xxC] = ballisticTrajectory(x0C,mission,windDirection,ttL(end));
plot3(xxC(1,:),xxC(2,:),xxC(3,:),'g')
plot3(x0(1),x0(2),x0(3),'bo')
plot3(rtarg(1),rtarg(2),rtarg(3),'ro')
err = norm(xxC(1:3,end)-rtarg)
%axis equal
theta = rad2deg(atan(X(:,end)./X(:,2)));
figure
plot(theta)
