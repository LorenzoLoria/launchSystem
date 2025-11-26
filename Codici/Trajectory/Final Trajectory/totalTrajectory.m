function [timeCollocation, stateCollocationGlobal] = totalTrajectory(mission,opt,thrustDataVec)

nStages = opt.nStages;
m0Tot = opt.m0Tot;
 
nDeval = 100;
stateCollocationGlobal = zeros(9,nDeval,nStages+1);
stateCollocation = zeros(7,nDeval,nStages);
timeCollocation = zeros(nDeval,nStages+1);
positionGlobalIRF = zeros(3,nDeval,nStages+1);
velocityGlobalIRF = zeros(3,nDeval,nStages+1);


for i = 1:nStages

        if i == 1
            m0 = m0Tot;
            x0 = [mission.environment.rEarth; 0; 0; 0; 0; 0; m0];
            t0 = 0;
        else
            m0 = m0-opt.stage{i-1}.mStage;
            x0 = [stateCollocation(1:6,end,i-1); m0];
            t0 = timeCollocation(end,i-1);
        end



   tf = opt.stage{i}.mProp * (length(thrustDataVec(:,1,i))-1)* opt.stage{i}.Isp * mission.environment.g0 * 2 ;
   tf = tf / (opt.stage{i}.nEngines *opt.stage{i}.engine.thrust) / (2*sum(thrustDataVec(:,1,i)) - thrustDataVec(1,1,i) - thrustDataVec(end,1,i)) ;
   tSpan = [t0 t0+tf]; %da rivedere nel caso i tempi non vadano bene
    
    tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));    
    fThrust = griddedInterpolant(tVec, thrustDataVec(:,:,i), 'linear', 'none'); 
    thrustData= @(t) fThrust(t).';

    opt.m0Tot = m0;
    
    [tt,xx] = launcherTrajectory(x0,mission,thrustData,tSpan,nDeval,i,opt);
    
    stateCollocation(:,:,i) = xx;
    timeCollocation(:,i) = tt;
end

localIRFtoGlobalIRF = [mission.initialPoint'/norm(mission.initialPoint) (cross(mission.environment.normalToTrajectoryPlane,mission.initialPoint'/norm(mission.initialPoint)))' mission.environment.normalToTrajectoryPlane' ]' ;
localIRFtoGlobalIRF(:,:,1:nStages) = localIRFtoGlobalIRF.*ones(size(localIRFtoGlobalIRF,1),size(localIRFtoGlobalIRF,2),nStages);
positionGlobalIRF(:,:,1:nStages) = pagemtimes(localIRFtoGlobalIRF , [stateCollocation(1:2,:,1:end);zeros(1,size(stateCollocation,2),nStages)]);
velocityGlobalIRF(:,:,1:nStages) = pagemtimes(localIRFtoGlobalIRF , [stateCollocation(3:4,:,1:end);zeros(1,size(stateCollocation,2),nStages)]);

x0Capsule = [positionGlobalIRF(:,end,nStages); velocityGlobalIRF(:,end,nStages)];
t0Capsule = timeCollocation(end,end-1);
windDirection = [1 0 0]';

[ttCaps,xxCaps] = ballisticTrajectory(x0Capsule,mission,windDirection,t0Capsule,nDeval);

positionGlobalIRF(:,:,end) = xxCaps(1:3,:);
velocityGlobalIRF(:,:,end) = xxCaps(4:6,:);
timeCollocation(:,end) = ttCaps;
stateCollocationGlobal(1:6,:,:) = [positionGlobalIRF;velocityGlobalIRF];
stateCollocationGlobal(7:8,:,1:nStages) = stateCollocation(5:6,:,:);
stateCollocationGlobal(9,:,1:nStages) = stateCollocation(end,:,:);
stateCollocationGlobal(9,:,end) = mission.capsule.weigth;
end





