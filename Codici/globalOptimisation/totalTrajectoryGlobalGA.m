function [timeCollocation, stateCollocation] = totalTrajectoryGlobalGA(launcher,opt,mission,settings,thrustDataVec)

persistent thrustDataVecPrev timeCollocationPrev stateCollocationPrev

if isempty(thrustDataVecPrev)
    thrustDataVecPrev = zeros(size(thrustDataVec)) ;
end

% if length(size(thrustDataVecPrev)) ~= length(size(thrustDataVec))
%     thrustDataVecPrev = zeros(size(thrustDataVec)) ; 
% end



if thrustDataVec == thrustDataVecPrev

    timeCollocation = timeCollocationPrev ;
    stateCollocation = stateCollocationPrev ;

else

    thrustDataVecPrev = thrustDataVec ;

    nStages = launcher(1);
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
            tf = opt.stage{i}.mProp * (length(thrustDataVec(:,1,i))-1)* opt.stage{i}.engine.ispZero * mission.environment.g0 * 2 ;
            thrustValue = opt.stage{i}.engine.thrustZero;
        else
            m0 = m0-opt.stage{i-1}.mStage;
            x0 = [stateCollocation(1:3,end,i-1); stateCollocation(4:6,end,i-1); m0];
            t0 = timeCollocation(end,i-1);
            tf = opt.stage{i}.mProp * (length(thrustDataVec(:,1,i))-1)* opt.stage{i}.engine.ispVac * mission.environment.g0 * 2 ;
            thrustValue = opt.stage{i}.engine.thrustVacum;
        end


        tf = tf / (opt.stage{i}.nEngines * thrustValue) / (2*sum(thrustDataVec(:,1,i)) - thrustDataVec(1,1,i) - thrustDataVec(end,1,i)) ;
        tSpan = [t0 t0+tf]; %da rivedere nel caso i tempi non vadano bene

        tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));
        fThrust = griddedInterpolant(tVec, thrustDataVec(:,:,i), 'linear', 'none');
        thrustData= @(t) fThrust(t).';

        opt.totalMass = m0;

        [tt,xx] = launcherTrajectory(x0,mission,thrustData,tSpan,nDeval,i,opt,settings.trajectoryOption2D);

        stateCollocation(:,:,i) = xx;
        timeCollocation(:,i) = tt;
    end

    x0Capsule = [stateCollocation(1:3,end,nStages); stateCollocation(4:6,end,nStages)];
    t0Capsule = timeCollocation(end,end-1);
    windDirection = [1 0 0]';

    [ttCaps,xxCaps] = ballisticTrajectory(x0Capsule,mission,windDirection,t0Capsule,nDeval);

    stateCollocation(:,:,end) = [xxCaps;mission.capsule.weight*ones(1,nDeval)];
    timeCollocation(:,end) = ttCaps;

    timeCollocationPrev = timeCollocation ;
    stateCollocationPrev = stateCollocation ;

end



end




