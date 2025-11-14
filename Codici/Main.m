
clearvars
clc
close all

%% Data
launcher.M0 = 210.4e3;

launcher.stage1.Mp=150.5e3;
launcher.stage1.Ms=12e3;
launcher.stage1.F=
launcher.stage1.mDot=

launcher.stage2.Mp=37.6e3;
launcher.stage2.Ms=3e3;
launcher.stage2.F=
launcher.stage2.mDot=

capsule.A=3.7^2*pi;
capsule.M=8600;
capsule.cd=1.23; %https://space.stackexchange.com/questions/55071/dragon-re-entry-flight-profile?utm_source=chatgpt.com

function [dsdt] = capsuledyn(x)

dsdt(1:2)=x(3:4);
vx = x(3);
vy = x(4);

v = sqrt(vx^2+vy^2);

