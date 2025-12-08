function [angleWrtVelMax,accMax] = findAccAngle(launcher,opt,mission,stateCollocation,timeCollocation,thrustDataVec,settings)

option2D = settings.trajectoryOption2D ; 
    for i = 1:launcher(1)
        time      = linspace(timeCollocation(1,i),timeCollocation(end,i),length(timeCollocation)-1);
        time2     = linspace(timeCollocation(1,i),timeCollocation(end,i),size(thrustDataVec,1));
        pos = [stateCollocation(1:3,:,i)];
        vel = [stateCollocation(4:6,:,i)];
        rNorm = sqrt(pos(1,:).^2+pos(2,:).^2+pos(3,:).^2);
        vNorm = sqrt(vel(1,:).^2+vel(2,:).^2+vel(3,:).^2);
        h   = rNorm-mission.environment.rEarth;  
        percThrustFun = griddedInterpolant(time2,thrustDataVec(:,1,i),'linear','linear');
        angleThrustFun = griddedInterpolant(time2,thrustDataVec(:,2,i),'linear','linear');
        if i == 1
        staticContribution = (101325-mission.environment.pressure(h))*opt.stage{i}.engine.effAreaZero;
        thrustValue = opt.stage{i}.engine.thrustZero;
        else
        staticContribution = 0;
        thrustValue = opt.stage{i}.engine.thrustVacum;
        end
        if option2D == 1
        thetaGimball = zeros(1,settings.nEvalPointsTraj);
        else
        thetaGimballFun = griddedInterpolant(time2,thrustDataVec(:,3,i),'linear','linear');
        thetaGimball = thetaGimballFun(timeCollocation(:,i)) ;
        end
    
        percThrustVec = percThrustFun(timeCollocation(:,i));
        gammaGimball = deg2rad(angleThrustFun(timeCollocation(:,i)));
    
        ThrustBRF = percThrustVec' .* opt.stage{i}.nEngines .*(thrustValue+staticContribution).* [cos(thetaGimball).*cos(gammaGimball'); cos(thetaGimball).*sin(gammaGimball'); sin(thetaGimball)];
        ThrustIRF = mission.target.Rfinal'*ThrustBRF;
        ThrustIRFNorm = sqrt(ThrustIRF(1,:).^2 + ThrustIRF(2,:).^2 + ThrustIRF(3,:).^2 );
        angleWrtVelMaxStage(i) = max(abs(acosd(sum(ThrustIRF .* vel, 1)./ThrustIRFNorm./vNorm)));
    
        accx      = diff(vel(1,:))./(time(2)-time(1));
        accy      = diff(vel(2,:))./(time(2)-time(1));
        accz      = diff(vel(3,:))./(time(2)-time(1));
        acc       = sqrt( accx.^2 + accy.^2 + accz.^2);
        accMaxStage(i) = max(acc);
    
    end
    angleWrtVelMax = max(angleWrtVelMaxStage);
    accMax = max(accMaxStage);

end