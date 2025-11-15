
clearvars
clc
close all

addpath(genpath("Capsule/"))

[mission] = dataStruct;
tSpan = [0 3*24*3600];
x0=[0;0;mission.environment.rEarth+100e3;5000;0;1000];
windDirection = -1 + 2*rand(3,1) ;
windDirection = windDirection / norm(windDirection) ; 

[tt,xx] = ballisticTrajectory(x0,mission,windDirection);
plot3(xx(1,:),xx(2,:),xx(3,:),"r", "LineWidth",1)
hold on

EarthPlot(mission.environment.rEarth)

<<<<<<< Updated upstream
=======
%%

hold on
plot3(x0(1),x0(2),x0(3), 'ob')
>>>>>>> Stashed changes
