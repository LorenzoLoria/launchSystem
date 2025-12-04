function [mission] = externalLoads(timeCollocation,stateCollocation,mission,opt,thrustData,alpha)

vel = stateCollocation(4:6,:,1:end-1)-stateCollocation(4:6,1,1);
vel = mission.Rfinal* vel(1:3,:);
absVel = sqrt ( vel(1,:).^2+ vel(2,:).^2 + vel(3,:).^2 );

normRocket = vel./absVel ;

pos = stateCollocation(1:3,:,1:end-1);
pos = pos(1:3,:);
absH = sqrt ( pos(1,:).^2+ pos(2,:).^2 + pos(3,:).^2 )-mission.environment.rEarth;

mass = stateCollocation(7, :, 1:end-1);
mass = mass(:);

rhoVec = mission.environment.rhoFun(absH);

q = 0.5*rhoVec.*(absVel).^2;

timeStage = timeCollocation(:,1:end-1);
timeStage = timeStage(:);


acc = [diff(vel(1,:))'./diff(timeStage),diff(vel(2,:))'./diff(timeStage),diff(vel(3,:))'./diff(timeStage)]';


[maxq,idx] = max(q);

vMaxQ = vel(:,idx);
aMaxQ = acc(:,idx);

v1 = vMaxQ/norm(vMaxQ);

v3 = [0;0;1];

v2 = cross(v3,v1)/norm(cross(v3,v1));

rot = [v1,v2,v3];

rot2 = [cos(alpha),-sin(alpha),0; sin(alpha),cos(alpha),0;0,0,1];

rot3 = rot2*rot;

vMaxQ = rot3'*vMaxQ;     %check
aMaxQ = rot3'*aMaxQ;
massMaxQ = mass(idx);
posMaxQ = pos(:,idx);
hMaxQ = absH(idx);
%dsdt(4:6) = (ThrustIRF + D ) / m + G;  

soundSpeed = mission.aerodynamics.soundspeedFun(hMaxQ);

machNumber = norm(vMaxQ) / soundSpeed;

[~,~,~,~, mainbodyCL, mainbodyCD, finsCL, finsCD] = CLCDcomputation(machNumber,alpha,maxq,1,mission);

dMaxQ = -rot2*maxq * mainbodyCD *mission.capsule.Area* [1;0;0];
lMaxQ = -rot2*maxq * mainbodyCL *mission.capsule.Area * [0;-1;0];
gMaxQ = -rot3'*mission.Rfinal* mission.environment.GM * posMaxQ /norm(posMaxQ)^3;

tMaxQ = (aMaxQ - gMaxQ)*massMaxQ - dMaxQ-lMaxQ;

% ==================== STRUCTURE DA ESTRARRE ==============================

mission.structure.dynamicPressure   = maxq;
mission.structure.dMaxQ             = dMaxQ;
mission.structure.lMaxQ             = lMaxQ;
mission.structure.aMaxQ             = aMaxQ;
mission.structure.tMaxQ             = tMaxQ;
mission.structure.massMaxQ          = massMaxQ;
mission.structure.gMaxQ             = gMaxQ;
mission.structure.hMaxQ             = hMaxQ;
mission.structure.vMaxQ             = vMaxQ;

m1Stage = opt.stage{1}.mStage + massMaxQ - opt.m0Tot;

mission.structure.massMaxQVec = [mission.capsule.weigth];

for ii = opt.nStages:-1:1
    mission.structure.massMaxQVec = [mission.structure.massMaxQVec, mission.structures{ii}.mInterstage, opt.stage{ii}.mStage];
end

mission.structure.massMaxQVec(end) = m1Stage;

% --- DA MODIFICARE

mission.structure.dragFinsMaxQ = -rot2*maxq * finsCD * mission.capsule.Area * [1;0;0];
mission.structure.liftFinsMaxQ = -rot2*maxq * finsCL * mission.capsule.Area * [0;-1;0];

end