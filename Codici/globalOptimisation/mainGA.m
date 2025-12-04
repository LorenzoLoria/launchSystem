clear all
clc
close all

%% Path Directory

addpath(genpath("..\..\"))

%% Upload Mission Struct
[mission,opt] = dataStruct;

%% GA

%opt.stage{1}.percentage = 0.5;
%opt.stage{2}.percentage = 0.4;


%1 --> 3
%2 --> 5
%3 --> 7



obj_ga = @(x) objFunGlobalGA( reshape(x,mission.optimisation.GA.variables,3,2), mission,opt);
nonlcon_ga = @(x) nlconGlobalGA( reshape(x,mission.optimisation.GA.variables,3,2), mission,opt );

lbGlobalGA = [1,1,1,1,0.2,0.2,0.2];
ubGlobalGA = [3,4,4,4,1,1,1];

[mer,staging] = initialMassEstimation(mission,opt);


