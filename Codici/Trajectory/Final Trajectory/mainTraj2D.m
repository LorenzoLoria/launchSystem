clear all
clc
close all

%% Path Directory

addpath(genpath("..\..\"))

%% Upload Mission Struct
tSpan = [0 500];

[mission,opt] = dataStruct;
%%

% %Generate optimisation variables for general stages

% thrustDataVec1 = [[1; 1 ;1; 1; 1] , [0; 0; 0; 0; 0] ];
% thrustDataVec2 = [[1; 1 ; 1; 1; 1] , [0; 30; 60; 70; 90] ];
% 
% thrustData(:,:,1) =thrustDataVec1;
% thrustData(:,:,2) =thrustDataVec2;
% 
% [timeCollocation, stateCollocation] = totalTrajectory(mission,opt,thrustData,1);

%%


% Initial Guess using GA
obj_ga = @(x) objFunMultiStagesGA( reshape(x,mission.optimisation.GA.variables,2,2), mission,opt,1);
nonlcon_ga = @(x) nlconMultiStagesGA( reshape(x,mission.optimisation.GA.variables,2,2), mission,opt,1);

lbFmincon(:,:,1) = [0.9*ones(mission.optimisation.GA.variables,1);0*ones(mission.optimisation.GA.variables,1)];
ubFmincon(:,:,1) = [ones(mission.optimisation.GA.variables,1);90*ones(mission.optimisation.GA.variables,1)];

lbFmincon(:,:,2) = [0.4*ones(mission.optimisation.GA.variables,1);0*ones(mission.optimisation.GA.variables,1)];
ubFmincon(:,:,2) = [ones(mission.optimisation.GA.variables,1);120*ones(mission.optimisation.GA.variables,1)];

lbGA = lbFmincon(:);
ubGA = ubFmincon(:);


options_ga = optimoptions("ga", ...
    "Display","iter", ...
    "MaxGenerations",20, ...
    "PopulationSize",100,...
    "UseParallel",true,...
    "FunctionTolerance", 1e-4,...
    'EliteCount',  6,...
    "CrossoverFraction", 0.9);

% n = 2*2*mission.optimisation.GA.variables;   % dimensione del vettore delle variabili
% 
% % blocchi di indici in cui vuoi x crescente
% blocks = [6 10;
%           15 20];    % righe: [inizio fine]
% 
% % numero di vincoli: (10-6) + (20-15) = 4 + 5 = 9
% m = sum(diff(blocks,1,2));
% 
% Aineq = zeros(m, n);
% row   = 0;
% 
% for b = 1:size(blocks,1)
%     iStart = blocks(b,1);
%     iEnd   = blocks(b,2);
%     for k = iStart:iEnd-1
%         row = row + 1;
%         Aineq(row, k)   = 1;
%         Aineq(row, k+1) = -1;
%     end
% end
% 
% bineq = zeros(m,1);


[x_ga, fval_ga] = ga(obj_ga,2*2*mission.optimisation.GA.variables,[],[],[],[],lbGA,ubGA,nonlcon_ga,options_ga);
%%

T0 = reshape(x_ga,mission.optimisation.GA.variables,2,2);

% Optimisation with fMinCon
[X,FVAL,EXITFLAG,OUTPUT] = fmincon(@(x) objFunMultiStages(x,mission,opt,1),T0,[],[],[],[],lbFmincon-eps,ubFmincon+eps,@(x) nlconMultiStages(x,mission,opt,1),mission.options.fmincon);


%% 
thrustData = X ;

[timeCollocation, stateCollocation] = totalTrajectory(mission,opt,thrustData,1);
figure(1)

EarthPlot(mission.environment.rEarth)
hold on
plot3(stateCollocation(1,:,1),stateCollocation(2,:,1),stateCollocation(3,:,1),'r')
plot3(stateCollocation(1,:,2),stateCollocation(2,:,2),stateCollocation(3,:,2),'y')
plot3(stateCollocation(1,:,3),stateCollocation(2,:,3),stateCollocation(3,:,3),'g')
plot3(mission.target.initialPointECI(1),mission.target.initialPointECI(2),mission.target.initialPointECI(3),'bo')
targetFinalLat = mission.target.latInitial ; 
targetFinalLon = mission.target.lonInitial + mission.target.omega * timeCollocation(end,end); 
targetFinalPosECI = 6371000*[cos(targetFinalLat)*cos(targetFinalLon); cos(targetFinalLat)*sin(targetFinalLon); sin(targetFinalLat) ];
plot3(targetFinalPosECI(1),targetFinalPosECI(2),targetFinalPosECI(3), 'ob')
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
title("Angle1")

X = [X(:,1,1) ; X(:,2,1) ; X(:,1,2) ;  X(:,2,2)] ;




if norm(stateCollocation(1:3,end,end) - mission.target.initialPointECI) < 2000
    
    filename = 'ThrustData.mat';
    if isfile(filename)
    result = load(filename,"X");
    X = [result.X,X];
    end
    % Risalvo tutto
    save(filename, 'X');
    
end
