function [mission] = externalLoads(timeCollocation,stateCollocation,mission,opt,thrustData,launcher,alpha)

vel = stateCollocation(4:6,:,1:end-1)-stateCollocation(4:6,1,1);
vel = mission.target.Rfinal* vel(1:3,:);
absVel = sqrt ( vel(1,:).^2+ vel(2,:).^2 + vel(3,:).^2 );

stageNumber = 1;

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

nStages = length(opt.stage);
geoStages = opt.geometry.stage;
optStage = opt.stage;

stageRadius = geoStages{nStages}.radius;
lCylinder = geoStages{nStages}.tanksLength; 
interstageLength = geoStages{nStages}.interstage.length;

dCylinder = lCylinder *stageRadius*2;

for i=nStages-1:-1:stageNumber
    lCylinder = lCylinder + geoStages{i}.tanksLength + geoStages{i}.interstage.length ;
    dCylinder = dCylinder + geoStages{i}.tanksLength * geoStages{i}.radius*2 +...
                geoStages {i}.interstage.length *...
                (geoStages {i}.radius + geoStages {i+1}.radius);
    
end

lNose = mission.capsule.height + geoStages{nStages}.interstage.length;
Anose = max(mission.capsule.Area , pi*stageRadius^2) ;
stage1Radius = geoStages{1}.radius;
boatTailRadius = geoStages{stageNumber}.radius ;
dimensions = [stageRadius,lCylinder,dCylinder,lNose,Anose,stage1Radius,boatTailRadius,interstageLength];
nEngines = optStage{stageNumber}.nEngines;
aTotZero = optStage{stageNumber}.engine.effAreaZero;
aTotVacum = optStage{stageNumber}.engine.effAreaVac;

engineVec = [nEngines,aTotZero,aTotVacum];

[CL,CD,CN,CA, mainbodyCL, mainbodyCD, finsCL, finsCD] = CLCDcomputation(machNumber,alpha,maxq,1,mission,1,dimensions,engineVec);

dMaxQ = -rot2*maxq * mainbodyCD *mission.capsule.Area* [1;0;0];
lMaxQ = -rot2*maxq * mainbodyCL *mission.capsule.Area * [0;-1;0];
gMaxQ = -rot3'*mission.target.Rfinal* mission.environment.GM * posMaxQ /norm(posMaxQ)^3;

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

m1Stage = opt.stage{1}.mStage + massMaxQ - opt.totalMass;

mission.structure.massMaxQVec = [mission.capsule.weight];

for ii = launcher(1):-1:1
    mission.structure.massMaxQVec = [mission.structure.massMaxQVec, mission.structures{ii}.mInterstage, opt.stage{ii}.mStage];
end

mission.structure.massMaxQVec(end) = m1Stage;

% --- DA MODIFICARE

mission.structure.dragFinsMaxQ = -rot2*maxq * finsCD * mission.aerodynamics.finsGeom.Se * [1;0;0];
mission.structure.liftFinsMaxQ = -rot2*maxq * finsCL * mission.aerodynamics.finsGeom.Se * [0;-1;0];

end