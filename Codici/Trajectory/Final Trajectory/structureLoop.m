%% Loop for Structure convergence

clear all
clc
close all

%% Path Directory

addpath(genpath("..\..\"))

%% Upload Mission Struct


[mission,opt] = dataStruct;


%% Initial Guess using GA
obj_ga = @(x) objFunMultiStagesGA( reshape(x,mission.optimisation.GA.variables,3,2), mission,opt);
nonlcon_ga = @(x) nlconMultiStagesGA( reshape(x,mission.optimisation.GA.variables,3,2), mission,opt );

lbFmincon(:,:,1) = [0.9*ones(mission.optimisation.GA.variables,1);-90*ones(mission.optimisation.GA.variables,1);0*ones(mission.optimisation.GA.variables,1)];
ubFmincon(:,:,1) = [ones(mission.optimisation.GA.variables,1);90*ones(mission.optimisation.GA.variables,1);100*ones(mission.optimisation.GA.variables,1)];

lbFmincon(:,:,2) = [0.3*ones(mission.optimisation.GA.variables,1);-90*ones(mission.optimisation.GA.variables,1);0*ones(mission.optimisation.GA.variables,1)];
ubFmincon(:,:,2) = [ones(mission.optimisation.GA.variables,1);90*ones(mission.optimisation.GA.variables,1);100*ones(mission.optimisation.GA.variables,1)];


lbGA = lbFmincon(:);
ubGA = ubFmincon(:);

options_ga = optimoptions("ga", ...
    "Display","iter", ...
    "MaxGenerations",20, ...
    "PopulationSize",100,...
    "UseParallel",true,...
    "FunctionTolerance", 1e-9); 


%% Loop 


err = 10000; %initialize error

whileStruct = struct();

whileStruct.stage{1} = opt.stage{1};
whileStruct.stage{2} = opt.stage{2};
whileStruct.stage{1}.engine = opt.stage{1}.engine ;
whileStruct.stage{2}.engine = opt.stage{2}.engine ;

mOld = opt.m0Tot;


%% retrive initial diameter and length

lOverD = 4 ;  

OF = whileStruct.stage{1}.engine.OF;
mOx = whileStruct.stage{1}.mProp * OF/(OF+1);
mFuel = whileStruct.stage{1}.mProp * 1/(OF+1);  
rhoFuel = whileStruct.stage{1}.engine.fuelDens ; 
rhoOx = whileStruct.stage{1}.engine.oxDens ;

fun1 =@(x) [ mFuel/rhoFuel/0.97-4/3*pi*x(1)^3 - pi*x(1)^2*x(2); mOx/rhoOx-4/3*pi*x(1)^3 - pi*x(1)^2*x(3); 0.8*lOverD- (4*x(1) + x(2)+ x(3))/(2*x(1)) ];
x01 = [2;10;10];
sol1 = fsolve(fun1,x01);

h1 = 4*sol1(1) + sol1(2) + sol1(3);


lOverD = 6;
OF = whileStruct.stage{2}.engine.OF;
mOx = whileStruct.stage{2}.mProp * OF/(OF+1);
mFuel = whileStruct.stage{2}.mProp * 1/(OF+1);  
rhoFuel = whileStruct.stage{2}.engine.fuelDens ; 
rhoOx = whileStruct.stage{2}.engine.oxDens ;

fun2 =@(x) [ mFuel/rhoFuel/0.97-4/3*pi*x(1)^3 - pi*x(1)^2*x(2); mOx/rhoOx-4/3*pi*x(1)^3 - pi*x(1)^2*x(3); 0.8*lOverD- (4*x(1) + x(2)+ x(3))/(2*x(1)) ];
x02 = [2;10;10];
sol2 = fsolve(fun2,x02);

h2 = 4*sol2(1) + sol2(2) + sol2(3);

%%
%while abs(err)>100

[x_ga, fval_ga] = ga(obj_ga,3*2*mission.optimisation.GA.variables,[],[],[],[],lbGA,ubGA,nonlcon_ga,options_ga);

T0 = reshape(x_ga,mission.optimisation.GA.variables,3,2);

[X,FVAL,EXITFLAG,OUTPUT] = fmincon(@(x) objFunMultiStages(x,mission,opt),T0,[],[],[],[],lbFmincon-eps,ubFmincon+eps,@(x) nlconMultiStages(x,mission,opt),mission.options.fmincon);


%% Retrive load conditions

thrustData = X;
[timeCollocation, stateCollocation] = totalTrajectory(mission,opt,thrustData);

[q,an,at,T,D,angle,gamma] = externalLoads(timeCollocation,stateCollocation,mission,opt,thrustData);

% .... codice di strutture ... 

%% Draw launcher