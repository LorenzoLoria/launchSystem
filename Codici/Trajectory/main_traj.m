clear all
clc
close all

%% Path Directory

addpath(genpath("..\..\"))

%% Upload Mission Struct

mission = dataStruct;

%% Optimisation Code

%% thetaGimball varia max da -pi/2 a pi/2 mentre gammaGImball da 0 a pi

% Initial Guess using GA
obj_ga = @(x) objFunGA( reshape(x,mission.optimisation.GA.variables,2), mission);
nonlcon_ga = @(x) nlconGA( reshape(x,mission.optimisation.GA.variables,2), mission );

lbGA = [0.1*ones(mission.optimisation.GA.variables,1);0*ones(mission.optimisation.GA.variables,1)];
ubGA = [ones(mission.optimisation.GA.variables,1);90*ones(mission.optimisation.GA.variables,1)];
%% GA initialisation

options_ga = optimoptions("ga", ...
    "Display","iter", ...
    "MaxGenerations",20, ...
    "PopulationSize",100,...
    "UseParallel",true,"HybridFcn","fmincon"); 

[x_ga, fval_ga] = ga(obj_ga,2*mission.optimisation.GA.variables,[],[],[],[],lbGA,ubGA,nonlcon_ga,options_ga);
%% FminCon final

T0 = reshape(x_ga,mission.optimisation.GA.variables,2);
% Optimisation with fMinCon
[X,FVAL,EXITFLAG,OUTPUT] = fmincon(@(x) objFun(x,mission),T0,[],[],[],[],reshape(lbGA,mission.optimisation.GA.variables,2),reshape(ubGA,mission.optimisation.GA.variables,2),@(x) nlcon(x,mission),mission.options.fmincon);







%% Initial Guess Plot

x0 = [mission.initialPoint'; 0; 0; 0; mission.launcher.engines{1}.m0];
tSpan = [0 1*24*3600];
thrustDataVec = reshape(x_ga,mission.optimisation.GA.variables,2);
%T0 = [0.9*ones(5,1),(-20)*ones(5,1)];
%thrustDataVec = T0;
tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));
thrustData = @(t) [interp1(tVec,thrustDataVec(:,1),t);interp1(tVec,thrustDataVec(:,2),t)];

[ttL,xxL] = launcherTrajectory(x0,mission, thrustData,0,thrustDataVec);
x0C = xxL(1:6,end);
windDirection = -1+2*rand(3,1);
windDirection = windDirection/norm(windDirection);
[ttC,xxC] = ballisticTrajectory(x0C,mission,windDirection,ttL(end));

figure(1)
plot3(xxL(1,:)/1000,xxL(2,:)/1000,xxL(3,:)/1000,'r.')
hold on
EarthPlot(mission.environment.rEarth/1000)
plot3(xxC(1,:)/1000,xxC(2,:)/1000,xxC(3,:)/1000,'r')
plot3(x0(1)/1000,x0(2)/1000,x0(3)/1000,'bo')
plot3(mission.target(1)/1000,mission.target(2)/1000,mission.target(3)/1000,'ro')

err = norm(xxC(1:3,end)-mission.target)

axis equal

%% Optimised Solution Plot


x0 = [mission.initialPoint'; 0; 0; 0; mission.launcher.engines{1}.m0];
tSpan = [0 3*24*3600];
thrustDataVec = X;
tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));
thrustData = @(t) [interp1(tVec,thrustDataVec(:,1),t);interp1(tVec,thrustDataVec(:,2),t)];

[ttL,xxL] = launcherTrajectory(x0,mission, thrustData,0,thrustDataVec);
x0C = xxL(1:6,end);
windDirection = -1+2*rand(3,1);
windDirection = windDirection/norm(windDirection);
[ttC,xxC] = ballisticTrajectory(x0C,mission,windDirection,ttL(end));

plot3(xxL(1,:)/1000,xxL(2,:)/1000,xxL(3,:)/1000,'g.')
plot3(xxC(1,:)/1000,xxC(2,:)/1000,xxC(3,:)/1000,'g')

hold off
err = norm(xxC(1:3,end)-mission.target)
axis equal
legend([ "" "Initial Guess Trajectory " "Initial Point" "Final Point" "Optimised Trajectory"])

theta = rad2deg(atan(X(:,end)./X(:,2)));
figure(2)
plot(X(:,1))
figure(3)
plot(X(:,2))

