function [timeCollocation, stateCollocation] = totalTrajectoryControlled(launcher,opt,mission,settings,windVelXFun,windVelYFun,stateCollocationRef,timeCollocationRef,maxGimball,thrustData,gainGA)

 

nStages = launcher(1);

% Retrieve the variables from the shuffled matrix 

% for j = 1:nStages
% 
%     opt.stage{j}.mProp = opt.montecarlo.massProp(parforiter,j);
%     opt.stage{j}.structuralMass = opt.montecarlo.massStruct(parforiter,j);
%     opt.stage{j}.Isp = opt.montecarlo.impSpec(parforiter,j);
%     opt.stage{j}.mStage = opt.montecarlo.StageMass(parforiter,j);
% 
% end


totalMass = opt.totalMass;

nDeval = settings.nEvalPointsTraj;
stateCollocation = zeros(7,nDeval,nStages+1);
timeCollocation = zeros(nDeval,nStages+1);

latInitial = mission.launchBase.latInitial;
lonInitial = mission.launchBase.lonInitial;

if settings.trajectoryOption2D == 1
    omega =0;
else
    omega = mission.target.omega ;
end
vxInitial = - omega * mission.environment.rEarth * cos(latInitial) * sin(lonInitial) ;
vyInitial = omega * mission.environment.rEarth * cos(latInitial) * cos(lonInitial) ;


for i = 1:nStages

    if i == 1
        m0 = totalMass;
        x0 = [mission.launchBase.initialPointECI'; vxInitial; vyInitial; 0;  m0];
        t0 = 0;
        tf = timeCollocationRef(end,i);
    else
        m0 = m0-opt.stage{i-1}.mStage;
        x0 = [stateCollocation(1:3,end,i-1); stateCollocation(4:6,end,i-1); m0];
        t0 = timeCollocation(end,i-1);
        tf = timeCollocationRef(end,i);
    end

    tSpan = [t0 t0+tf]; %da rivedere nel caso i tempi non vadano bene

    opt.totalMass = m0;
    stateCollocationRefIdx = stateCollocationRef(:,:,i);
    timeCollocationRefIdx = timeCollocationRef(:,i);
    refStateFun = griddedInterpolant(timeCollocationRefIdx,stateCollocationRefIdx(1:6,:)','linear','linear');
    tvec = linspace(timeCollocationRefIdx(1),timeCollocationRefIdx(end),size(thrustData,1));
    thrustDataFun = griddedInterpolant(tvec,thrustData,"linear","none");
    [tt,xx] = launcherTrajectoryControlled(x0,mission,tSpan,nDeval,i,opt,1,windVelXFun,windVelYFun,refStateFun,maxGimball,thrustDataFun,gainGA);

    stateCollocation(:,:,i) = xx;
    timeCollocation(:,i) = tt;
end

x0Capsule = [stateCollocation(1:3,end,nStages); stateCollocation(4:6,end,nStages)];
t0Capsule = timeCollocation(end,end-1);
windDirection = [1 0 0]';

[ttCaps,xxCaps] = ballisticTrajectory(x0Capsule,mission,windDirection,t0Capsule,nDeval);

stateCollocation(:,:,end) = [xxCaps;mission.capsule.weight*ones(1,nDeval)];
timeCollocation(:,end) = ttCaps;



end
