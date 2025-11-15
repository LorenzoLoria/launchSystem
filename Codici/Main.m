
clearvars
clc
close all

[mission] = dataStruct;
tSpan = [0 3*24*3600];
x0=[0;0;mission.envirnoment.rEarth+100e3;8000;0;0];
[tt,xx] = ballisticTrajectory(x0,mission);
plot3(xx(1,:),xx(2,:),xx(3,:),"r", "LineWidth",1)
hold on

EarthPlot(mission.envirnoment.rEarth)

