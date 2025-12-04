function[cin,ceq] = nlconMultiStagesGA(x,mission,opt,option2D)

thrustDataVec = x;

[timeCollocation,stateCollocation] = totalTrajectory(mission,opt,thrustDataVec,option2D);

latInitial = mission.target.latInitial ;
latFinal = latInitial ;
lonInitial = mission.target.lonInitial ;
if option2D == 1
    omega = 0;

else
    omega = mission.target.omega ;
    
end
lonFinal = lonInitial + omega * timeCollocation(end,end) ;
targetFinalPos = 6371000*[cos(latFinal)*cos(lonFinal); cos(latFinal)*sin(lonFinal); sin(latFinal) ];


for i = 1:opt.nStages
    time      = linspace(timeCollocation(1,i),timeCollocation(end,i),length(timeCollocation)-1);
    time2     = linspace(timeCollocation(1,i),timeCollocation(end,i),size(thrustDataVec,1));
    pos = [stateCollocation(1:3,:,i)];
    vel = [stateCollocation(4:6,:,i)];
    rNorm = sqrt(pos(1,:).^2+pos(2,:).^2+pos(3,:).^2);
    vNorm = sqrt(vel(1,:).^2+vel(2,:).^2+vel(3,:).^2);
    h   = rNorm-mission.environment.rEarth;  
    rho = mission.environment.rhoFun(h);
    percThrustFun = griddedInterpolant(time2,thrustDataVec(:,1,i),'linear','linear');
    angleThrustFun = griddedInterpolant(time2,thrustDataVec(:,2,i),'linear','linear');
    if i == 1
    staticContribution = (101325-mission.environment.pressure(h))*opt.stage{i}.engine.effAreaZero;
    isp = opt.stage{i}.engine.ispZero;
    thrustValue = opt.stage{i}.engine.thrustZero;
    else
    staticContribution = 0;
    isp = opt.stage{i}.engine.ispVac;
    thrustValue = opt.stage{i}.engine.thrustVacum;
    end
    if option2D == 1
    thetaGimball = zeros(1,opt.nDeval);
    else
    thetaGimballFun = griddedInterpolant(time2,thrustDataVec(:,3,i),'linear','linear');
    thetaGimball = thetaGimballFun(timeCollocation(:,i)) ;
    end

    percThrustVec = percThrustFun(timeCollocation(:,i));
    gammaGimball = deg2rad(angleThrustFun(timeCollocation(:,i)));

    ThrustBRF = percThrustVec' .* opt.stage{i}.nEngines .*(thrustValue+staticContribution).* [cos(thetaGimball).*cos(gammaGimball'); cos(thetaGimball).*sin(gammaGimball'); sin(thetaGimball)];
    ThrustIRF = mission.Rfinal'*ThrustBRF;
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
%cin = [(accMax-6*mission.environment.g0)/mission.environment.g0 ; (norm(stateCollocation(1:3,end,end) - targetFinalPos) - 5e3)/1e3];
%ceq = norm(stateCollocation(1:3,end,end) - targetFinalPos);
cin = [(accMax-6*mission.environment.g0)/mission.environment.g0 ];
ceq = [];
end