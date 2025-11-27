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

% thrustDataVec1 = [[1; 1 ;1; 1; 1] , [0; 0; 0; 0; 0] ];
% thrustDataVec2 = [[1; 1 ; 1; 1; 1] , [0; 60; 70; 80; 90] ];
% 
% thrustData(:,:,1) =thrustDataVec1;
% thrustData(:,:,2) =thrustDataVec2;
% 
% [timeCollocation, stateCollocation] = totalTrajectory(mission,opt,thrustData);

%%

% Initial Guess using GA
obj_ga = @(x) objFunMultiStagesGA( reshape(x,mission.optimisation.GA.variables,3,2), mission,opt);
nonlcon_ga = @(x) nlconMultiStagesGA( reshape(x,mission.optimisation.GA.variables,3,2), mission,opt );

lbFmincon(:,:,1) = [0.9*ones(mission.optimisation.GA.variables,1);-90*ones(mission.optimisation.GA.variables,1);0*ones(mission.optimisation.GA.variables,1)];
ubFmincon(:,:,1) = [ones(mission.optimisation.GA.variables,1);90*ones(mission.optimisation.GA.variables,1);100*ones(mission.optimisation.GA.variables,1)];

lbFmincon(:,:,2) = [0.2*ones(mission.optimisation.GA.variables,1);-90*ones(mission.optimisation.GA.variables,1);0*ones(mission.optimisation.GA.variables,1)];
ubFmincon(:,:,2) = [ones(mission.optimisation.GA.variables,1);90*ones(mission.optimisation.GA.variables,1);100*ones(mission.optimisation.GA.variables,1)];


lbGA = lbFmincon(:);
ubGA = ubFmincon(:);

ind1(:,:,1) = [0.9880   -2.7485    2.3756...
    1.0000   -1.7611   21.2746...
    0.9832   -1.0682    5.9076...
    0.9110   26.3551   44.3236...
    0.9252    5.9692   34.9518] ; 
ind1(:,:,2) = [0.3156 0.2109 50.9289...
    0.8615   -6.1012   88.9859...
    0.9844    6.0644   72.1261...
    0.9876    9.2176   85.1640...
    0.5500   81.1500   77.5762] ; 
ind1 = ind1(:) ;
ind2(:,:,1) = [1.0000    0.9680    8.2502...
    0.9926   -1.1427    9.2594...
    0.9189   -3.4637   18.1407...
    0.9000   -8.2757   48.5878...
    0.9000   -5.5477   65.7544];
ind2(:,:,2) = [0.3824   -1.9564   38.6008...
    0.8894  -51.3270   67.9861...
    0.9827  -69.4701   55.9126...
    1.0000  -60.7321   83.4617...
    0.3000  -52.5123   54.4304];
ind2 = ind2(:) ;

ind3(:,:,1) = [0.9851   -0.6762    0.0946...
    0.9981   -1.2033   33.7250...
    0.9484   -1.5424    8.6267...
    0.9704    2.3461   16.3068...
    0.9554    2.9494   92.4029];

ind3(:,:,2) = [0.9267  -39.5366   83.2831...
    0.9896    5.2593   76.5156...
    0.9358   11.8028   74.9192...
    0.8802  -72.6178   99.9643...
    0.4762    8.0744   83.6577];
ind3 = ind3(:) ;

ind4(:,:,1) = [0.9584   -0.7362    3.0108...
    0.9617   -7.1524   13.0258...
    0.9526    2.4299    9.6667...
    0.9510    5.8943   35.5056...
    0.9523   -7.2826   65.8391];

ind4(:,:,2) = [0.5684   -5.7325   44.2025...
    0.6112  -30.2430   77.1205...
    0.6345   30.8746   74.5043...
    0.6585  -58.0053   80.6452...
    0.3684   30.8006   82.8541];
ind4 = ind4(:) ;

ind5(:,:,1) = [0.9852  -17.5905   10.0802...
    0.9218   23.7661    6.2794...
    0.9708  -22.2833   13.0833...
    0.9906   16.1841   34.2804...
    0.9927    9.4728   59.8558];

ind5(:,:,2) = [0.8634   24.8162   71.7692...
    0.9959   14.4336   72.0577...
    0.8019   16.0032   82.4152...
    0.9518   31.0518   87.0766...
    0.2050   27.2781   67.9287];
ind5 = ind5(:) ;

ind6(:,:,1) = [0.9956   -2.7727   13.0488...
    0.9243   -2.3127    4.7023...
    0.9720   -2.4859   27.0284...
    0.9607   13.3108   22.9273...
    0.9868  -33.8890   83.4509];

ind6(:,:,2) = [0.7494   26.9874   29.4877...
    0.7732  -22.4865   93.7943...
    0.8558   41.5260   77.4466...
    0.8221   30.4242   77.8149...
    0.4019   35.0883   87.7839];
ind6 = ind6(:) ;



initialPopulation = [ind1' ; ind2' ; ind3' ; ind4' ; ind5'; ind6'];


%% GA initialisation

options_ga = optimoptions("ga", ...
    "Display","iter", ...
    "MaxGenerations",30, ...
    "PopulationSize",200,...
    "UseParallel",true,...
    "FunctionTolerance", 1e-9,...
    'EliteCount',  6,...
    "InitialPopulationMatrix", initialPopulation);
    % 'SelectionFcn', {@selectiontournament, 4},...
    %     'MutationFcn', {@mutationpower},...
    %     'CrossoverFcn', {@crossoverintermediate, 0.8},...
    %     "InitialPopulationMatrix", initialPopulation); 

[x_ga, fval_ga] = ga(obj_ga,3*2*mission.optimisation.GA.variables,[],[],[],[],lbGA,ubGA,nonlcon_ga,options_ga);
%%
T0 = reshape(x_ga,mission.optimisation.GA.variables,3,2);
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

figure
plot([X(:,3,1);X(:,3,2)])
title("Angle2")

% dist = sqrt ( stateCollocation(1,:,end).*stateCollocation(1,:,end) + stateCollocation(2,:,end).*stateCollocation(2,:,end)+stateCollocation(3,:,end).*stateCollocation(3,:,end));
% h = dist-6371000;
% 
% figure
% plot(h)