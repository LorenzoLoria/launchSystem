clear all
clc
close all

%% Path Directory

addpath(genpath("..\..\"))

%% Upload Mission Struct
[mission,opt,settings] = dataStruct;

%% GA

%opt.stage{1}.percentage = 0.5;
%opt.stage{2}.percentage = 0.4;


%1 --> 3
%2 --> 5
%3 --> 7


% obj_ga = @(x) objFunGlobalGA( reshape(x,mission.optimisation.GA.variables,3,2), mission,opt);
% nonlcon_ga = @(x) nlconGlobalGA( reshape(x,mission.optimisation.GA.variables,3,2), mission,opt );
 
lbGlobalGA = [1,1,1,1,0.2,0.2,0.2];
ubGlobalGA = [3,4,4,4,1,1,1];
 
intcon = [1 2 3 4];


option2D = 1;
options_ga = optimoptions("ga", ...
    "Display","iter", ...
    "MaxGenerations",20, ...
    "PopulationSize",50,...
    "UseParallel",true,...
    "FunctionTolerance", 1e-4);


[x_ga, fval_ga] = ga(@(x) objFunGlobalGA(x,mission,opt,settings,option2D) ,7,[],[],[],[],lbGlobalGA,ubGlobalGA, @(x) nlconGlobalGA(x,mission,opt,settings,option2D),intcon,options_ga);
