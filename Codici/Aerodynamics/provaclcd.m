clear all
clc
close all

%% Path Directory

addpath(genpath("..\"))

%% Upload Mission Struct
tSpan = [0 500];

[mission,opt] = dataStruct;


[CL,CD] = CLCDcomputation(2,deg2rad(3),283220,1,mission)