clearvars
clc
close all

%% Path Directory

addpath(genpath("..\..\"))

%% Upload Mission Struct

mission = dataStruct;

%% Optimisation Code

% Initial fMinCon Guess
xx0 = [0*ones(6,1) , 4*845e3*0.5*cos(deg2rad(70)) .* ones(6,1), 4*844e3*0.5*sin(deg2rad(70)) .* ones(6,1)];

% Target Point
rtarg = [0;mission.environment.rEarth * cos(deg2rad(75.5));mission.environment.rEarth * sin(deg2rad(75.5))];

% Setting Boundaries
lb = 0*ones(size(xx0,1));
ub = 4*845e3*ones(size(xx0,1),size(xx0,2));

disp('xx0:');
disp(xx0);
disp('LB:');
disp(lb);
disp('UB:');
disp(ub);

% Optimisation with fMinCon
[X,FVAL,EXITFLAG,OUTPUT] = fmincon(@(x) objFun(x,mission),xx0,[],[],[],[],-1*ones(size(xx0,1),size(xx0,2)),4*845e3*ones(size(xx0,1),size(xx0,2)),@(x) nlcon(x,mission,rtarg),mission.options.fmincon)

%% Initial Guess Plot

x0 = [0; 0; mission.environment.rEarth; 0; 0; 0; mission.launcher.engines{1}.m0];
tSpan = [0 1*24*3600];
thrustDataVec = xx0;
tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));
thrustData = @(t) [interp1(tVec,thrustDataVec(:,1),t);interp1(tVec,thrustDataVec(:,2),t);interp1(tVec,thrustDataVec(:,3),t)];

[ttL,xxL] = launcherTrajectory(x0,mission, thrustData);
x0C = xxL(1:6,end);
windDirection = -1+2*rand(3,1);
windDirection = windDirection/norm(windDirection);
[ttC,xxC] = ballisticTrajectory(x0C,mission,windDirection,ttL(end));

figure(1)
plot3(xxL(1,:)/1000,xxL(2,:)/1000,xxL(3,:)/1000,'r.')
hold on
%EarthPlot(6371)
plot3(xxC(1,:)/1000,xxC(2,:)/1000,xxC(3,:)/1000,'r')
plot3(x0(1)/1000,x0(2)/1000,x0(3)/1000,'bo')
plot3(rtarg(1)/1000,rtarg(2)/1000,rtarg(3)/1000,'ro')

err = norm(xxC(1:3,end)-rtarg)

axis equal
%% Optimised Solution Plot


x0 = [0; 0; mission.environment.rEarth; 0; 0; 0; mission.launcher.engines{1}.m0];
tSpan = [0 3*24*3600];
thrustDataVec = X;
tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));
thrustData = @(t) [interp1(tVec,thrustDataVec(:,1),t);interp1(tVec,thrustDataVec(:,2),t);interp1(tVec,thrustDataVec(:,3),t)];

[ttL,xxL] = launcherTrajectory(x0,mission, thrustData);
x0C = xxL(1:6,end);
windDirection = -1+2*rand(3,1);
windDirection = windDirection/norm(windDirection);
[ttC,xxC] = ballisticTrajectory(x0C,mission,windDirection,ttL(end));

plot3(xxL(1,:)/1000,xxL(2,:)/1000,xxL(3,:)/1000,'g.')
plot3(xxC(1,:)/1000,xxC(2,:)/1000,xxC(3,:)/1000,'g')
hold off
err = norm(xxC(1:3,end)-rtarg)
axis equal
legend([ "" "Initial Guess Trajectory " "Initial Point" "Final Point" "Optimised Trajectory"])

theta = rad2deg(atan(X(:,end)./X(:,2)));
figure(2)
plot(theta)



