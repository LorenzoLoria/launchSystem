clear all
clc
close all

%% Path Directory

addpath(genpath("..\..\"))

%% Upload Mission Struct
tSpan = [0 500];

[mission,opt] = dataStruct;
%%

% Generate optimisation variables for general stages

thrustDataVec1 = [[1; 1 ;1; 1; 1] , [0; 0; 0; 0; 0] ];
thrustDataVec2 = [[1; 1 ; 1; 1; 1] , [0; 60; 70; 80; 90] ];

thrustData(:,:,1) =thrustDataVec1;
thrustData(:,:,2) =thrustDataVec2;

[timeCollocation, stateCollocation] = totalTrajectory(mission,opt,thrustData);

%%

% Initial Guess using GA
obj_ga = @(x) objFunMultiStagesGA( reshape(x,mission.optimisation.GA.variables,2,2), mission,opt);
nonlcon_ga = @(x) nlconMultiStagesGA( reshape(x,mission.optimisation.GA.variables,2,2), mission,opt );

lbFmincon(:,:,1) = [0.9*ones(mission.optimisation.GA.variables,1);0*ones(mission.optimisation.GA.variables,1)];
ubFmincon(:,:,1) = [ones(mission.optimisation.GA.variables,1);150*ones(mission.optimisation.GA.variables,1)];

lbFmincon(:,:,2) = [0.3*ones(mission.optimisation.GA.variables,1);0*ones(mission.optimisation.GA.variables,1)];
ubFmincon(:,:,2) = [ones(mission.optimisation.GA.variables,1);150*ones(mission.optimisation.GA.variables,1)];


lbGA = lbFmincon(:);
ubGA = ubFmincon(:);


%% GA initialisation

options_ga = optimoptions("ga", ...
    "Display","iter", ...
    "MaxGenerations",20, ...
    "PopulationSize",100,...
    "UseParallel",true); 

[x_ga, fval_ga] = ga(obj_ga,2*2*mission.optimisation.GA.variables,[],[],[],[],lbGA,ubGA,nonlcon_ga,options_ga);
%%
T0 = reshape(x_ga,mission.optimisation.GA.variables,2,2);
% Optimisation with fMinCon
[X,FVAL,EXITFLAG,OUTPUT] = fmincon(@(x) objFunMultiStages(x,mission,opt),T0,[],[],[],[],lbFmincon-eps,ubFmincon+eps,@(x) nlconMultiStages(x,mission,opt),mission.options.fmincon);


%%

thrustData = X;

[timeCollocation, stateCollocation] = totalTrajectory(mission,opt,thrustData);
figure(1)

EarthPlot(mission.environment.rEarth)
hold on
plot3(stateCollocation(1,:,1),stateCollocation(2,:,1),stateCollocation(3,:,1),'r')
plot3(stateCollocation(1,:,2),stateCollocation(2,:,2),stateCollocation(3,:,2),'y')
plot3(stateCollocation(1,:,3),stateCollocation(2,:,3),stateCollocation(3,:,3),'g')
plot3(mission.target(1),mission.target(2),mission.target(3),'bo')
title("Trajectory")
hold off


%%
figure
plot(timeCollocation(:),[stateCollocation(7,:,1),stateCollocation(7,:,2) ,stateCollocation(7,:,3)]/1000)
title("Mass")


figure
plot(diff([stateCollocation(7,:,1),stateCollocation(7,:,2) ,stateCollocation(7,:,3)]))
title("Mass flow rate")


figure
for i = 1:opt.nStages
    time      = linspace(timeCollocation(1,i),timeCollocation(end,i),length(timeCollocation)-1);
    v         = stateCollocation(4:6,:,i);
    vNorm     = sqrt( stateCollocation(4,:,i).^2 + stateCollocation(5,:,i).^2 +stateCollocation(6,:,i).^2 );
    acc       = diff(vNorm)./(time(2)-time(1));
    plot(time,acc./mission.environment.g0)
    hold on
end
title("Accelerations")


figure
plot([X(:,1,1);X(:,1,2)])
title("Throttling")



figure
plot([X(:,2,1);X(:,2,2)])
title("Angle")

% dist = sqrt ( stateCollocation(1,:,end).*stateCollocation(1,:,end) + stateCollocation(2,:,end).*stateCollocation(2,:,end)+stateCollocation(3,:,end).*stateCollocation(3,:,end));
% h = dist-6371000;
% 
% figure
% plot(h)