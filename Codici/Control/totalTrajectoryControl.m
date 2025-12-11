function [timeCollocationControlled, stateCollocationControlled] = totalTrajectoryControl(mission,mer,configuration,settings, launcher,thrustDataVec, guidancePoints, guidanceTime, gains)

nStages = launcher(1);
m0Tot = configuration.totalMass;

nDeval = settings.nEvalPointsTraj;
stateCollocationControlled = zeros(14,nDeval, nStages+1);
timeCollocationControlled  = zeros(nDeval, nStages+1);


for i = 1:nStages

        if i == 1
            m0 = m0Tot;
            q0 = [1; 0; 0; 0] ;
            x0 = [mission.launchBase.initialPointECI'; 0; 0; 0; q0; 0; 0; 0; m0];
            t0 = 0;
            tf = configuration.stage{i}.mProp * (length(thrustDataVec(:,1,i))-1)* configuration.stage{i}.engine.ispZero * mission.environment.g0 * 2 ;
            thrustValue = configuration.stage{i}.engine.thrustZero;
            guidancePointsStage = squeeze(guidancePoints(1:end-1, :,i)) ;
            guidanceTimeStage = guidanceTime(:,i) ;
        else
            m0 = m0 - configuration.stage{i-1}.mStage;
            x0 = [stateCollocationControlled(1:end-1,end,i-1); m0];
            t0 = timeCollocationControlled(end,i-1);
            tf = configuration.stage{i}.mProp * (length(thrustDataVec(:,1,i))-1)* configuration.stage{i}.engine.ispVac * mission.environment.g0 * 2 ;
            thrustValue = configuration.stage{i}.engine.thrustVacum;
            guidancePointsStage = squeeze(guidancePoints(1:end-1,:, i)) ;
            guidanceTimeStage = guidanceTime(:,i) ;
        end

    tf = tf / (configuration.stage{i}.nEngines * thrustValue) / (2*sum(thrustDataVec(:,1,i)) - thrustDataVec(1,1,i) - thrustDataVec(end,1,i)) ;
    tSpan = [t0 t0+tf]; 
    
    tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));    
    fThrust = griddedInterpolant(tVec, thrustDataVec(:,:,i), 'linear', 'none'); 
    thrustData = @(t) fThrust(t).';

    configuration.totalMass = m0;
    
    [tt,xx] = launcherTrajectoryControl(x0, mission, mer, configuration, launcher, thrustData,tSpan,nDeval,i, guidancePointsStage, guidanceTimeStage, gains);
    
    stateCollocationControlled(:,:,i) = xx;
    timeCollocationControlled(:,i) = tt;
end

x0Capsule = [stateCollocationControlled(1:3,end,nStages); stateCollocationControlled(4:6,end,nStages)];
t0Capsule = timeCollocationControlled(end,end-1);
windDirection = [1 0 0]';

[ttCaps,xxCaps] = ballisticTrajectory(x0Capsule,mission,windDirection,t0Capsule,nDeval);

stateCollocationControlled(1:6,:,end) = xxCaps ;
stateCollocationControlled(end, :, end) = mission.capsule.weight*ones(1,nDeval) ;
timeCollocationControlled(:,end) = ttCaps;

end