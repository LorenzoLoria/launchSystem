function [timeCollocation, stateCollocation] = totalTrajectory(mission,opt,thrustData,thrustDataVec)

nStages = opt.nStages;

m0Tot = opt.stage{1}.mStage +opt.stage{2}.mStage + mission.capsule.weigth; 
nDeval = 100;
stateCollocation = zeros(7,nDeval,nStages+1);
timeCollocation = zeros(nDeval,nStages+1);

for i = 1:nStages
    
if i == 1
    m0 = m0Tot;
    x0 = [mission.initialPoint'; 0 ;0 ; 0;  m0];
    t0 = 0;
else
    m0 = m0-opt.stage{i-1}.mStage;
    x0 = [stateCollocation(1:3,end,i-1); stateCollocation(4:6,end,i-1); m0];
    t0 = timeCollocation(end,i-1);
end

[tt,xx] = launcherTrajectory(x0,mission,thrustData,t0,thrustDataVec,nDeval);

stateCollocation(:,:,i) = xx;
timeCollocation(:,i) = tt;
end

x0Capsule = [stateCollocation(1:3,end,nStages); stateCollocation(4:6,end,nStages)];
t0Capsule = timeCollocation(end,end-1);
windDirection = [1 0 0]';

[ttCaps,xxCaps] = ballisticTrajectory(x0Capsule,mission,windDirection,t0Capsule,nDeval);

stateCollocation(:,:,end) = [xxCaps;mission.capsule.weigth*ones(1,nDeval)];
timeCollocation(:,end) = ttCaps;

end





