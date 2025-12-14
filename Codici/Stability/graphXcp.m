clc
close all
clear all

% Path Directory

addpath(genpath("..\..\"))

launcher = [2,2,3,3,0.45,0.77, 0.5] ;

% Path Directory

addpath(genpath("..\..\"))

% Upload Mission Struct
[mission,settings] = dataStructGlobal;

for i = 1:launcher(1)
    configuration.stage{i}.engine = mission.engines{launcher(1+i)};
end

[mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);

%%
[xGATraj, fvalGATraj] = ga( @(x) objFunGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings), ...
                        launcher(1)*2*settings.nOptPointsTraj,...
                        [],[],[],[],settings.lowerBoundsGA,settings.upperBoundsGA, ...
                        @(x) nlconGATraj( reshape(x,settings.nOptPointsTraj,2,launcher(1)),launcher,configuration, mission,settings),settings.gaTrajOptions);


thrustData = reshape(xGATraj,settings.nOptPointsTraj,2,2);


[timeCollocation, stateCollocation] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustData);

%% PARAMETERS

t = timeCollocation(:);
m = [stateCollocation(end, :, 1) stateCollocation(end, :, 2) stateCollocation(end, :, 3)];

% Flare data

Xcg_curve = zeros(size(t));
Xcp_curve = zeros(size(t));
Xcp_flare_curve = zeros(size(t));
Xcp_finflare_curve = zeros(size(t));
Xcp_fins_curve = zeros(size(t));
d_curve = zeros(size(t));

%% Computations
for i = 1:length(t)
    
    if i <= 100
        stageNumber = 1;
        collocationIndex = i;
          
        
    elseif i <= 200
        stageNumber = 2;
        collocationIndex = i-100;
        
     

    end
  
    stateVector = stateCollocation(1:6, collocationIndex, stageNumber) ;
    d = 2*configuration.geometry.stage{stageNumber}.radius ;
    d_curve(i) = d;

    Xcg = computeXcgGlobal(mission, configuration, launcher, mer, m(i), stageNumber);
    Xcg_curve(i) = Xcg;
   
    Xcp = computeXcpGlobal(mission, configuration,launcher, 0, stageNumber);
    Xcp_curve(i) = Xcp;

    Xcp_flare = computeFlareXcp( mission, configuration, d);
    Xcp_flare_curve(i) = Xcp_flare;

    Xcp_finflare = computeFinFlareXcp(mission,configuration, d, stageNumber);
    Xcp_finflare_curve(i) = Xcp_finflare;

 
    Xcp_fins = computetotalFinXcp(configuration, mission, stageNumber, stateVector) ;
    Xcp_fins_curve(i) = Xcp_fins;

end

% Find the indices where t <= timeCollocation(end,2)
valid_indices = t <= timeCollocation(end,2);

% Truncate t and Xcg_curve
t_truncated = t(valid_indices);
Xcg_curve_truncated = Xcg_curve(valid_indices);
Xcp_curve_truncated = Xcp_curve(valid_indices);
Xcp_flare_curve_truncated = Xcp_flare_curve(valid_indices);
Xcp_finflare_curve_truncated = Xcp_finflare_curve(valid_indices);
Xcp_fins_curve_truncated = Xcp_fins_curve(valid_indices);
d_curve_truncated = d_curve(valid_indices);


%% PLOT CURVES 
figure('Position', [100, 100, 800, 600]);

% Plot all curves on the same graph
hold on;
plot(t_truncated, Xcg_curve_truncated, 'LineWidth', 2, 'DisplayName', 'Xcg');
plot(t_truncated, Xcp_curve_truncated, 'LineWidth', 2, 'DisplayName', 'Xcp');
plot(t_truncated, Xcp_flare_curve_truncated, 'LineWidth', 2, 'DisplayName', 'Xcp with Flare');
plot(t_truncated, Xcp_finflare_curve_truncated, 'LineWidth', 2, 'DisplayName', 'Xcp with Fins and Flare');
plot(t_truncated, Xcp_fins_curve_truncated, 'LineWidth', 2, 'DisplayName', 'Xcp with Fins');
hold off;

% Add labels and title
xlabel('Time (s)');
ylabel('Position (m)');
title('Center of Gravity and Center of Pressure vs. Time (Truncated)');
grid on;

legend('show', 'Location', 'best');

%% PLOT STABILITY MARGIN 

stabilityMargin = (Xcp_curve_truncated - Xcg_curve_truncated) ./ d_curve_truncated;
stabilityMargin_fins = (Xcp_fins_curve_truncated - Xcg_curve_truncated) ./ d_curve_truncated;
stabilityMargin_flare = (Xcp_flare_curve_truncated - Xcg_curve_truncated) ./ d_curve_truncated;
stabilityMargin_finflare = (Xcp_finflare_curve_truncated - Xcg_curve_truncated) ./ d_curve_truncated;

figure;

% Subplot 2: Stability margin vs. Time
subplot(4, 1, 1);
plot(t_truncated, stabilityMargin, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('(Xcp-Xcg)/d');
title('Stability margin vs. Time');
grid on;

% Subplot 2: Stability margin with fins vs. Time
subplot(4, 1, 2);
plot(t_truncated, stabilityMargin_fins, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('(Xcp-Xcg)/d');
title('Stability margin with fins vs. Time');
grid on;

% Subplot 3: Stability margin with flare vs. Time
subplot(4, 1, 3);
plot(t_truncated, stabilityMargin_flare, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('(Xcp-Xcg)/d');
title('Stability margin with flare vs. Time');
grid on;

% Subplot 4: Stability margin with fins and flare vs. Time
subplot(4, 1, 4);
plot(t_truncated, stabilityMargin_finflare, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('(Xcp-Xcg)/d');
title('Stability margin with fins and flare vs. Time');
grid on;




function Xcg = computeXcgGlobal(mission, configuration, launcher, mer, m, stageNumber)
% Computes the center of gravity Xcg of a three-stage launcher with a nose cone
% Assumption: Uniform Weight Distribution 
%
% Inputs:
%   N        : number of stages 
%   lc1, lc2 : length of stage 1 and 2
%   lco      : length of conical nose
%   m        : mass of the launcher at time t 
%   mc2      : wet mass of stage 2 
%   mco      : wet mass of payload 
%
% Output:
%   Xcg      : center of gravity measured from the top of the launcher

N = launcher(1) - (stageNumber - 1);

% Payload Data
% m = massMaxQ;
lco = mission.capsule.height;
mco = mission.capsule.weight;

lc = zeros(1,N);
li = zeros(1,N);
mStage = zeros(1,N);
x = lco;
massVec = mco;
totalMassIter = mission.capsule.weight ;

for i = (launcher(1)-stageNumber+1) : -1 : 1
  
    totalMassIter = totalMassIter + configuration.stage{i}.mStage ;
    lc(i) = configuration.geometry.stage{i}.tanksLength;
    li(i) = configuration.geometry.stage{i}.interstage.length;

    if i == 1
        mStage(i) = configuration.stage{i}.mStage - (totalMassIter - m) ;
    else
        mStage(i) = configuration.stage{i}.mStage;

    end
    if i == launcher(1)
        rBottom = configuration.geometry.stage{i}.radius; % bottom = lower stage
        rTop = mission.capsule.radius;
    else
        rBottom = configuration.geometry.stage{i}.radius; % bottom = lower stage
        rTop = configuration.geometry.stage{i+1}.radius;
    end

    mInterstage(i) = mer.stage{i}.interStage;
    x = [x li(i) lc(i)];
    massVec = [massVec mInterstage(i) mStage(i)];
end

xi = zeros(1,length(x));


for i = 1:length(x)
    if i==1
        xi(i) = sum(x(1:i))-x(i)/4;
    elseif isinteger(i/2)

        xi(i) = sum(x(1:i)) - x(i) * (rBottom^2 + 2*rBottom*rTop + 3*rTop^2) / (4 * (rBottom^2 + rBottom*rTop + rTop^2));
    else
        xi(i) = sum(x(1:i))-x(i)/2;
    end
end

Xcg = (massVec*xi') / m;

end