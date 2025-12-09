clc
close all
clear all

% Optimal Solution form an old GA
launcher = [2,1,4,4,0.4056,0.4016,0.7];
for i = 1:launcher(1)
    configuration.stage{i}.engine = mission.engines{launcher(1+i)};
end

[mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);

[xGATraj, fvalGATraj] = ga( @(x) objFunGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings), ...
                        launcher(1)*2*settings.nOptPointsTraj,...
                        [],[],[],[],settings.lowerBoundsGA,settings.upperBoundsGA, ...
                        @(x) nlconGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings),settings.gaTrajOptions);

%% Plots
thrustData = reshape(xGATraj,settings.nOptPointsTraj,2,2);

% Upload state and time from Traj2D
load('stateCollocation.mat') ;
load('timeCollocation.mat') ;

% Initial Condition
x0 = stateCollocation(:,1,1) ;

% Upload mission data
[mission,opt] = dataStruct;
tSpan = [0 500];

Nstage = 3 ;

idx = ceil(linspace(1, 100 , 25)) ;
% Point for the contorl law of the Guidance
guidancePoints = zeros(7,length(idx), Nstage) ; 
guidanceTime = zeros(length(idx), Nstage) ;


for ii = 1 : Nstage
 guidancePoints(:,:,ii) = stateCollocation(:, idx, ii) ;
 guidanceTime(:,ii)   = timeCollocation(idx, ii) ;
end

% prova di gains
gains = [1, 1, 1, 1] ;


%%
