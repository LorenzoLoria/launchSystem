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

for i = (launcher(1)-stageNumber+1) : -1 : 1
    lc(i) = configuration.geometry.stage{i}.tanksLength;
    li(i) = configuration.geometry.stage{i}.interstage.length;
    if i == 1
        mStage(i) = configuration.stage{i}.mStage - (configuration.totalMass - m) ;
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