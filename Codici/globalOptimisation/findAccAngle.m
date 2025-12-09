function [angleWrtVelMax,accMax] = findAccAngle(launcher,opt,mission,stateCollocation,timeCollocation,thrustDataVec,settings)


option2D = settings.trajectoryOption2D ; 
nDeval = settings.nEvalPointsTraj;
nEl = settings.nOptPointsTraj;
cosangleWrtVelMaxStage = zeros(launcher(1),1);
accMaxStage = zeros(launcher(1),1);
RIRF = mission.target.Rfinal.';

    for i = 1:launcher(1)
        dtime      = timeCollocation(2,i) - timeCollocation(1,i);
        time2     = timeCollocation(1:nDeval/(nEl):nDeval,i);
        pos = stateCollocation(1:3,:,i);
        vel = stateCollocation(4:6,:,i);
        rNorm = sqrt(pos(1,:).^2+pos(2,:).^2+pos(3,:).^2);
        vNorm = sqrt(vel(1,:).^2+vel(2,:).^2+vel(3,:).^2);
        percThrustVec = interp1(time2,thrustDataVec(:,1,i),timeCollocation(:,i),'linear','extrap');
        gammaGimball = interp1(time2,thrustDataVec(:,2,i),timeCollocation(:,i),'linear','extrap');
        if i == 1
        h   = rNorm-mission.environment.rEarth;
        staticContribution = (101325-mission.environment.pressure(h))*opt.stage{i}.engine.effAreaZero;
        thrustValue = opt.stage{i}.engine.thrustZero;
        else
        staticContribution = 0;
        thrustValue = opt.stage{i}.engine.thrustVacum;
        end
        if option2D == 1
        thetaGimball = zeros(1,settings.nEvalPointsTraj);
        else
        thetaGimball = interp1(time2,thrustDataVec(:,3,i),timeCollocation(:,i),'linear','extrap');
        end
    
        ThrustBRF = percThrustVec' .* opt.stage{i}.nEngines .*(thrustValue+staticContribution).* [cos(thetaGimball).*cos(gammaGimball'); cos(thetaGimball).*sin(gammaGimball'); sin(thetaGimball)];
        ThrustIRF = RIRF*ThrustBRF;
        ThrustIRFNorm = sqrt(ThrustIRF(1,:).^2 + ThrustIRF(2,:).^2 + ThrustIRF(3,:).^2 );
        cosangleWrtVelMaxStage(i) = min(sum(ThrustIRF .* vel, 1)./ThrustIRFNorm./vNorm);

        accx      = diff(vel(1,:))./(dtime);
        accy      = diff(vel(2,:))./(dtime);
        accz      = diff(vel(3,:))./(dtime);
        acc       = sqrt( accx.^2 + accy.^2 + accz.^2);
        accMaxStage(i) = max(acc);
    end
    angleWrtVelMax = acosd(min(cosangleWrtVelMaxStage));
    accMax = max(accMaxStage);

end