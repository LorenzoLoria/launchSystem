function [timeCollocation, stateCollocation] = totalTrajectory(mission,opt,thrustDataVec,option2D)

persistent thrustDataVecPrev timeCollocationPrev stateCollocationPrev

if isempty(thrustDataVecPrev)
    thrustDataVecPrev = zeros(size(thrustDataVec)) ;
end


if thrustDataVec == thrustDataVecPrev

timeCollocation = timeCollocationPrev ;
stateCollocation = stateCollocationPrev ;

else

thrustDataVecPrev = thrustDataVec ; 

nStages = opt.nStages;
m0Tot = opt.m0Tot;
 
nDeval = 100;
stateCollocation = zeros(7,nDeval,nStages+1);
timeCollocation = zeros(nDeval,nStages+1);

latInitial = mission.launchBase.latInitial;
lonInitial = mission.launchBase.lonInitial;

if option2D == 1
    omega =0;
else
    omega = mission.target.omega ;
end
vxInitial = - omega * mission.environment.rEarth * cos(latInitial) * sin(lonInitial) ;
vyInitial = omega * mission.environment.rEarth * cos(latInitial) * cos(lonInitial) ;


for i = 1:nStages

        if i == 1
            m0 = m0Tot;
            x0 = [mission.launchBase.initialPointECI'; vxInitial; vyInitial; 0;  m0];
            t0 = 0;
            tf = opt.stage{i}.mProp * (length(thrustDataVec(:,1,i))-1)* opt.stage{i}.engine.ispZero * mission.environment.g0 * 2 ;

        else
            m0 = m0-opt.stage{i-1}.mStage;
            x0 = [stateCollocation(1:3,end,i-1); stateCollocation(4:6,end,i-1); m0];
            t0 = timeCollocation(end,i-1);
             tf = opt.stage{i}.mProp * (length(thrustDataVec(:,1,i))-1)* opt.stage{i}.engine.ispVac * mission.environment.g0 * 2 ;

        end


     tf = tf / (opt.stage{i}.nEngines *opt.stage{i}.engine.thrust) / (2*sum(thrustDataVec(:,1,i)) - thrustDataVec(1,1,i) - thrustDataVec(end,1,i)) ;
     tSpan = [t0 t0+tf]; %da rivedere nel caso i tempi non vadano bene
    
    tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));    
    fThrust = griddedInterpolant(tVec, thrustDataVec(:,:,i), 'linear', 'none'); 
    thrustData= @(t) fThrust(t).';

    opt.m0Tot = m0;
    
    [tt,xx] = launcherTrajectory(x0,mission,thrustData,tSpan,nDeval,i,opt,option2D);
    
    stateCollocation(:,:,i) = xx;
    timeCollocation(:,i) = tt;
end

x0Capsule = [stateCollocation(1:3,end,nStages); stateCollocation(4:6,end,nStages)];
t0Capsule = timeCollocation(end,end-1);
windDirection = [1 0 0]';

[ttCaps,xxCaps] = ballisticTrajectory(x0Capsule,mission,windDirection,t0Capsule,nDeval);

stateCollocation(:,:,end) = [xxCaps;mission.capsule.weigth*ones(1,nDeval)];
timeCollocation(:,end) = ttCaps;

timeCollocationPrev = timeCollocation ;
stateCollocationPrev = stateCollocation ; 

end



end




