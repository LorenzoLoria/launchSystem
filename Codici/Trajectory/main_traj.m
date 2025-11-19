clearvars
clc
close all

%% Path Directory

addpath(genpath("..\..\"))

%% Upload Mission Struct

mission = dataStruct;

%% Optimisation Code

%% thetaGimball varia max da -pi/2 a pi/2 mentre gammaGImball da 0 a pi


% Target Point
rtarg = [0;mission.environment.rEarth * cos(deg2rad(75));mission.environment.rEarth * sin(deg2rad(75))];

% Initial Guess using GA
obj_ga = @(x) objFunGA( reshape(x,5,3), mission, rtarg);
nonlcon_ga = @(x) nlconGA( reshape(x,5,3), mission );

%lbGA = 845e3*ones(length(Tvec),1);
%ubGA = 4*845e3*ones(length(Tvec),1);

lbGA = [0.001*ones(5,1);-45*ones(5,1);-45*ones(5,1)];
ubGA = [ones(5,1);45*ones(5,1);45*ones(5,1)];
%%
options_ga = optimoptions("ga", ...
    "Display","iter", ...
    "MaxGenerations",20, ...
    "PopulationSize",50); 

[x_ga, fval_ga] = ga(obj_ga,15,[],[],[],[],lbGA,ubGA,nonlcon_ga,options_ga)
%%
T0 = [0.9*ones(5,1),deg2rad(0)*ones(5,1),deg2rad(-35).*ones(5,1)] %reshape(x_ga,5,3)
%%


% Optimisation with fMinCon
[X,FVAL,EXITFLAG,OUTPUT] = fmincon(@(x) objFun(x,mission,rtarg),T0,[],[],[],[],reshape(lbGA,5,3),reshape(ubGA,5,3),@(x) nlcon(x,mission,rtarg),mission.options.fmincon)

%% Initial Guess Plot

x0 = [0; 0; mission.environment.rEarth; 0; 0; 0; mission.launcher.engines{1}.m0];
tSpan = [0 1*24*3600];
thrustDataVec = T0;
tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));
thrustData = @(t) [interp1(tVec,thrustDataVec(:,1),t);interp1(tVec,thrustDataVec(:,2),t);interp1(tVec,thrustDataVec(:,3),t)];

[ttL,xxL] = launcherTrajectory(x0,mission, thrustData,0);
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
plot3(rtarg(1)/1000,rtarg(2)/1000,rtarg(3)/1000,'ro')

err = norm(xxC(1:3,end)-rtarg)

axis equal

%% Optimised Solution Plot


x0 = [0; 0; mission.environment.rEarth; 0; 0; 0; mission.launcher.engines{1}.m0];
tSpan = [0 3*24*3600];
thrustDataVec = X;
tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));
thrustData = @(t) [interp1(tVec,thrustDataVec(:,1),t);interp1(tVec,thrustDataVec(:,2),t);interp1(tVec,thrustDataVec(:,3),t)];

[ttL,xxL] = launcherTrajectory(x0,mission, thrustData,0);
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
plot(X(:,1))
figure(3)
plot(X(:,2))
hold on
plot(X(:,3))

hold off


