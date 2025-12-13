function [maxQData] = externalLoads(timeCollocation,stateCollocation,mission,configuration,launcher,mer, alpha, staging)

vel = stateCollocation(4:6,:,1:end-1)-stateCollocation(4:6,1,1);
vel = mission.target.Rfinal* vel(1:3,:);
absVel = sqrt ( vel(1,:).^2+ vel(2,:).^2 + vel(3,:).^2 );

stageNumber = 1;


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
% acc(:,end+1) = acc(:,end);
acc = [acc, acc(:,end)] ;
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

% m1Stage = configuration.stage{1}.mStage - mer.stage{1}.interStage + massMaxQ - configuration.totalMass;
mPropUsed = configuration.totalMass - massMaxQ;
mPropRemaining = staging{1}.mProp - mPropUsed;
mFuRemaining = mPropRemaining / (1 + configuration.stage{1}.engine.OF);
mOxRemaining = mPropRemaining * configuration.stage{1}.engine.OF / (1 + configuration.stage{1}.engine.OF);

maxQData.massMaxQVec = [mission.capsule.weight];

for ii = launcher(1):-1:1
    maxQData.massMaxQVec = [maxQData.massMaxQVec, mer.stage{ii}.interStage, ...
        mer.stage{ii}.tankMassFuel+mer.stage{ii}.cryoInsuFuel+staging{ii}.mProp/(1 + configuration.stage{ii}.engine.OF), ...
        mer.stage{ii}.tankMassOx+mer.stage{ii}.cryoInsuOx+staging{ii}.mProp*configuration.stage{ii}.engine.OF/(1 + configuration.stage{ii}.engine.OF), ...
        mer.stage{ii}.thrustFrame + mer.stage{ii}.engineWeight + ...
        mer.stage{ii}.avionics + mer.stage{ii}.wiring + mer.stage{ii}.tvc +  ...
        mer.stage{ii}.battery + mer.stage{ii}.pressurant];
end

maxQData.massMaxQVec(end-2) = mer.stage{1}.tankMassFuel + mFuRemaining + mer.stage{1}.cryoInsuFuel;
maxQData.massMaxQVec(end-1) = mer.stage{1}.tankMassOx + mOxRemaining + mer.stage{1}.cryoInsuOx;

soundSpeed = mission.aerodynamics.soundspeedFun(hMaxQ);

machNumber = norm(vMaxQ) / soundSpeed;

nStages = length(configuration.stage);
geoStages = configuration.geometry.stage;
configurationStage = configuration.stage;

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
Aref  = mission.capsule.Area;
stage1Radius = geoStages{1}.radius;
boatTailRadius = geoStages{stageNumber}.radius ;
dimensions = [stageRadius,lCylinder,dCylinder,lNose,Anose,stage1Radius,boatTailRadius,interstageLength];
nEngines = configurationStage{stageNumber}.nEngines;
aTotZero = configurationStage{stageNumber}.engine.effAreaZero;
aTotVacum = configurationStage{stageNumber}.engine.effAreaVac;
finsVec   = [mission.aerodynamics.finsGeom.rootChord, mission.aerodynamics.finsGeom.tipChord, ...
    mission.aerodynamics.finsGeom.bfin, mission.aerodynamics.finsGeom.Nfins, mission.aerodynamics.finsGeom.Sfin, ...
    mission.aerodynamics.finsGeom.cmac, mission.aerodynamics.finsGeom.delta_le, ...
    mission.aerodynamics.finsGeom.Lambda_le, mission.aerodynamics.finsGeom.tmac];
engineVec = [nEngines,aTotZero,aTotVacum];

[~,~,~,~, mainbodyCL, mainbodyCD, finsCL, finsCD] = CLCDcomputation(machNumber,alpha,maxq,1,mission,1,dimensions,engineVec, finsVec);

dMaxQ = -rot2*maxq * mainbodyCD * Aref * [1;0;0];
lMaxQ = -rot2*maxq * mainbodyCL * Aref * [0;-1;0];
gMaxQ = -rot3'*mission.target.Rfinal* mission.environment.GM * posMaxQ /norm(posMaxQ)^3;
alphaFin = 10 * pi / 180;
dFinsMaxQ = -rot2*maxq * finsCD * Aref * [cos(alphaFin);sin(alphaFin);0];
lFinsMaxQ = -rot2*maxq * finsCL * Aref * [-sin(alphaFin);cos(alphaFin);0];

tMaxQ = (aMaxQ - gMaxQ)*massMaxQ - dMaxQ-lMaxQ - dFinsMaxQ - lFinsMaxQ;

% Position computation
maxQData.hMaxQ = hMaxQ;
maxQData.vMaxQ = vMaxQ;
maxQData.massMaxQ = massMaxQ;

xcp = computeXcp(mission, configuration,launcher);
xcp_a = mission.aerodynamics.finsGeom.rootChord - computeFinXcp(mission, maxQData);

% Length of the element used for structural analysis
h = [mission.capsule.height/2, xcp - mission.capsule.height/2, mission.capsule.height - xcp];

for ii = launcher(1):-1:1
        h = [ h, configuration.geometry.stage{ii}.interstage.length/2, ...
            configuration.geometry.stage{ii}.interstage.length/2, ...
            configuration.stage{ii}.fuelTankH/2, configuration.stage{ii}.fuelTankH/2, ...
            configuration.stage{ii}.oxTankH/2, configuration.stage{ii}.oxTankH/2, ...
            configuration.geometry.stage{ii}.thrustFrame/2,configuration.geometry.stage{ii}.thrustFrame/2];
end

%CG evaluation
h4CG = cumsum(h);
h4CG = h4CG([1,3:end]);
h4cgFinal = h4CG(1:2:end);
xcg = centerOfGravity(maxQData.massMaxQVec, h4cgFinal );


% Deflection angle
delta = asin((-(lFinsMaxQ(2) + dFinsMaxQ(2)) * (configuration.geometry.totalLength - xcp_a - xcg) + (dMaxQ(2) + lMaxQ(2)) * (xcg - xcp))/(norm(tMaxQ)*(configuration.geometry.totalLength - xcg)));

% Rotation of thrust
tMaxQ = sqrt(tMaxQ'* tMaxQ) * [cos(delta); sin(delta); 0];

% Rotation of the acceleration 
aMaxQ = gMaxQ + (dMaxQ+lMaxQ+dFinsMaxQ+lFinsMaxQ+tMaxQ) / massMaxQ;

% ==================== STRUCTURE DA ESTRARRE ==============================

maxQData.dynamicPressure   = maxq;
maxQData.dMaxQ             = dMaxQ;
maxQData.lMaxQ             = lMaxQ;
maxQData.aMaxQ             = aMaxQ;
maxQData.tMaxQ             = tMaxQ;
maxQData.massMaxQ          = massMaxQ;
maxQData.gMaxQ             = gMaxQ;
maxQData.hMaxQ             = hMaxQ;
maxQData.vMaxQ             = vMaxQ;
maxQData.liftFinsMaxQ      = lFinsMaxQ;
maxQData.dragFinsMaxQ      = dFinsMaxQ;
end