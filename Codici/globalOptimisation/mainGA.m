clear all
clc
close all

%% Path Directory

addpath(genpath("..\..\"))

%% Upload Mission Struct
[mission,settings] = dataStructGlobal;

%% GA


[x_ga, fval_ga] = ga(@(x) objFunGlobalGA(x,mission,settings),...
                    settings.globalGAoptVariables,[],[],[],[],...
                    settings.lowerBoundsGlobalGA,settings.upperBoundsGlobalGA,...
                    @(x) nlconGlobalGA(x,mission,settings),...
                    settings.intconGlobalGA,settings.globalGAOptions);


%% PLOTYS Sfonv9uwn'enu0n erqni0











