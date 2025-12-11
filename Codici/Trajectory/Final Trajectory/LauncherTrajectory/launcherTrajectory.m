function [tt,xx] = launcherTrajectory(x0,mission,thrustData,tSpan,nDeval,stageNumber,opt,option2D)

%persistent copyThrustData copytt copyxx

%if isempty(copyThrustData)
%copyThrustData = zeros(mission.optimisation.GA.variables,2);
%end

%if thrustDataVec == copyThrustData 
%tt = copytt;
%xx = copyxx;
%else 

% ThrustDataVec deve essere una matrice n*3

% Manca la propEvent che deve riferirsi a quando finisce il carburante
% dello stadio o qualcosa di simile.

options = odeset('RelTol',1e-6,'AbsTol',1e-2 ,'Events',@(t,x) propEvent(t,x, mission,opt.stage{stageNumber}.mProp,opt,stageNumber));


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

finsVec   = [mission.aerodynamics.finsGeom.rootChord, mission.aerodynamics.finsGeom.tipChord, ...
    mission.aerodynamics.finsGeom.bfin, mission.aerodynamics.finsGeom.Nfins, mission.aerodynamics.finsGeom.Sfin, ...
    mission.aerodynamics.finsGeom.cmac, mission.aerodynamics.finsGeom.delta_le, ...
    mission.aerodynamics.finsGeom.Lambda_le, mission.aerodynamics.finsGeom.tmac];

solution = ode113(@(t,x) launcherDynamicsECI(t, x,thrustData, mission,stageNumber,opt,option2D,dimensions,engineVec, finsVec), tSpan, x0,options);

tt = linspace(solution.x(1),solution.x(end),nDeval);

xx = deval(solution,tt);

%copytt = tt;
%copyxx = xx;
%copyThrustData = thrustDataVec;
%end
end