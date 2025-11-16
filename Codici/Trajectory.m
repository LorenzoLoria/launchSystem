clearvars
clc
close all

addpath(genpath("Capsule"))
addpath(genpath("LauncherTrajectory"))
mission = dataStruct;

x0 = [0; 0; mission.environment.rEarth; 0; 0; 0; mission.launcher.engines{1}.m0];
tSpan = [0 3*24*3600];
thrustDataVec = [zeros(5,1) , 4*845e3*cos(deg2rad(60)) .* ones(5,1),4*845e3*sin(deg2rad(60)) .* ones(5,1)];
tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));
thrustData = @(t) [interp1(tVec,thrustDataVec(:,1),t);interp1(tVec,thrustDataVec(:,2),t);interp1(tVec,thrustDataVec(:,3),t)];


[ttL,xxL] = launcherTrajectory(x0,mission, thrustData);
plot3(xxL(1,:),xxL(2,:),xxL(3,:))
hold on
x0C = xxL(1:6,end);
windDirection = -1+2*rand(3,1);
windDirection = windDirection/norm(windDirection);

[ttC,xxC] = ballisticTrajectory(x0C,mission,windDirection);
plot3(xxC(1,:),xxC(2,:),xxC(3,:))
%EarthPlot(mission.environment.rEarth)
axis equal
