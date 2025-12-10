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

mission.structure.machNumber = norm(vMaxQ) / soundSpeed;

[~,~,~,~, mainbodyCL, mainbodyCD, mission.structure.finsCL, mission.structure.finsCD] = CLCDcomputation(mission.structure.machNumber,alpha,maxq,1,mission);

dMaxQ = -rot2*maxq * mainbodyCD *mission.capsule.Area* [1;0;0];
lMaxQ = -rot2*maxq * mainbodyCL *mission.capsule.Area * [0;-1;0];
gMaxQ = -rot3'*mission.Rfinal* mission.environment.GM * posMaxQ /norm(posMaxQ)^3;
alphaFin = 10 * pi / 180;
dFinsMaxQ = -rot2*maxq * mission.structure.finsCD * mission.aerodynamics.finsGeom.Se * [sin(alphaFin);cos(alphaFin);0];
lFinsMaxQ = -rot2*maxq * mission.structure.finsCL * mission.aerodynamics.finsGeom.Se * [cos(alphaFin);-sin(alphaFin);0];

tMaxQ = (aMaxQ - gMaxQ)*massMaxQ - dMaxQ-lMaxQ - dFinsMaxQ - lFinsMaxQ;

% Position computation
mission.structure.hMaxQ = hMaxQ;
mission.structure.vMaxQ = vMaxQ;
mission.structure.massMaxQ = massMaxQ;
xcp = computeXcp(mission, opt);
xcp_a = mission.aerodynamics.rootChord - computeFinXcp(mission); % compute position starting from the launcher bottom
xcg = computeXCG(mission, opt);

% Deflection angle
delta = asin((-(lFinsMaxQ(2) + dFinsMaxQ(2)) * (mission.structure.launcherLength - xcp_a - xcg) + (dMaxQ(2) + lMaxQ(2)) * (xcg - xcp))/(norm(tMaxQ)*(mission.structure.launcherLength - xcg)));
% delta = 0;
% Rotation of thrust
tMaxQ = sqrt(tMaxQ(1)^2 + tMaxQ(2)^2) * [cos(delta); sin(delta); 0];

% Rotation of the acceleration 
aMaxQ = gMaxQ + (dMaxQ+lMaxQ+dFinsMaxQ+lFinsMaxQ+tMaxQ) / massMaxQ;


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
mission.structure.dragFinsMaxQ      = dFinsMaxQ;
mission.structure.liftFinsMaxQ      = lFinsMaxQ;

m1Stage = opt.stage{1}.mStage + massMaxQ - opt.m0Tot;

mission.structure.massMaxQVec = [mission.capsule.weigth];

for ii = opt.nStages:-1:1
    mission.structure.massMaxQVec = [mission.structure.massMaxQVec, mission.structures{ii}.mInterstage, opt.stage{ii}.mStage];
end

mission.structure.massMaxQVec(end) = m1Stage;

end