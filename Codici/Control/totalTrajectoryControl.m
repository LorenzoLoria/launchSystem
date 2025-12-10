function [timeCollocationControlled, stateCollocationControlled] = totalTrajectoryControl(mission,mer,configuration,launcher,thrustDataVec, guidancePoints, guidanceTime, gains)

nStages = launcher(1);
m0Tot = configuration.totalMass;

nDeval = setting.nEvalPointsTraj;
stateCollocationControlled = zeros(7,nDeval, nStages);
timeCollocationControlled  = zeros(nDeval, nStages+1);


for i = 1:nStages

        if i == 1
            m0 = m0Tot;
            x0 = [mission.environment.rEarth; 0; 0; 0; 0; 0; 0; 0; m0];
            t0 = 0;
            tf = configuration.stage{i}.mProp * (length(thrustDataVec(:,1,i))-1)* configuration.stage{i}.engine.ispZero * mission.environment.g0 * 2 ;
            thrustValue = configuration.stage{i}.engine.thrustZero;
            guidancePointsStage = guidancePoints(:, 1:nDeval) ;
            guidanceTimeStage = guidanceTime(1:nDeval) ;
        else
            m0 = m0 - configuration.stage{i-1}.mStage;
            x0 = [stateCollocationControlled(1:8,end,i-1); m0];
            t0 = timeCollocationControlled(end,i-1);
            tf = configuration.stage{i}.mProp * (length(thrustDataVec(:,1,i))-1)* configuration.stage{i}.engine.ispVac * mission.environment.g0 * 2 ;
            thrustValue = configuration.stage{i}.engine.thrustVacum;
            guidancePointsStage = guidancePoints(:,(i-1)*nDeval : i*nDeval) ;
            guidanceTimeStage = guidanceTime((i-1)*nDeval : i*nDeval) ;
        end

    tf = tf / (configuration.stage{i}.nEngines * thrustValue) / (2*sum(thrustDataVec(:,1,i)) - thrustDataVec(1,1,i) - thrustDataVec(end,1,i)) ;
    tSpan = [t0 t0+tf]; 
    
    tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));    
    fThrust = griddedInterpolant(tVec, thrustDataVec(:,:,i), 'linear', 'none'); 
    thrustData = @(t) fThrust(t).';

    configuration.totalMass = m0;
    
    [tt,xx] = launcherTrajectoryControl(x0,mission, mer, configuration, launcher, thrustData,tSpan,nDeval,i, guidancePointsStage, guidanceTimeStage, gains);
    
    stateCollocationControlled(:,:,i) = xx;
    timeCollocationControlled(:,i) = tt;
end

x0Capsule = [stateCollocationControlled(1:3,end,nStages); stateCollocationControlled(4:6,end,nStages)];
t0Capsule = timeCollocationControlled(end,end-1);
windDirection = [1 0 0]';

[ttCaps,xxCaps] = ballisticTrajectory(x0Capsule,mission,windDirection,t0Capsule,nDeval);

stateCollocationControlled(:,:,end) = [xxCaps;mission.capsule.weigth*ones(1,nDeval)];
timeCollocationControlled(:,end) = ttCaps;

end