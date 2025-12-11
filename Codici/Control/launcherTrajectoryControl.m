function [tt,xx] = launcherTrajectoryControl(x0,mission, mer, configuration, launcher, ~,tSpan,nDeval, stageNumber, guidancePoints, guidanceTime, gains)
 
nStages = launcher(1);

stageRadius = configuration.geometry.stage{nStages}.radius;
lCylinder = configuration.geometry.stage{nStages}.tanksLength; 
interstageLength = configuration.geometry.stage{nStages}.interstage.length;

dCylinder = lCylinder *stageRadius*2;

for i=nStages-1:-1:stageNumber
    lCylinder = lCylinder + configuration.geometry.stage{i}.tanksLength + configuration.geometry.stage{i}.interstage.length ;
    dCylinder = dCylinder + configuration.geometry.stage{i}.tanksLength * configuration.geometry.stage{i}.radius*2 +...
                configuration.geometry.stage {i}.interstage.length *...
                (configuration.geometry.stage {i}.radius + configuration.geometry.stage {i+1}.radius);
    
end

lNose = mission.capsule.height + configuration.geometry.stage{nStages}.interstage.length;
Anose = max(mission.capsule.Area , pi*stageRadius^2) ;
stage1Radius = configuration.geometry.stage{1}.radius;
boatTailRadius = configuration.geometry.stage{stageNumber}.radius ;
dimensions = [stageRadius,lCylinder,dCylinder,lNose,Anose,stage1Radius,boatTailRadius,interstageLength];
nEngines  = configuration.stage{stageNumber}.nEngines;
aTotZero  = configuration.stage{stageNumber}.engine.effAreaZero;
aTotVacum = configuration.stage{stageNumber}.engine.effAreaVac;

engineVec = [nEngines,aTotZero,aTotVacum];

finsVec   = [mission.aerodynamics.finsGeom.rootChord, mission.aerodynamics.finsGeom.tipChord, ...
    mission.aerodynamics.finsGeom.bfin, mission.aerodynamics.finsGeom.Nfins, mission.aerodynamics.finsGeom.Sfin, ...
    mission.aerodynamics.finsGeom.cmac, mission.aerodynamics.finsGeom.delta_le, ...
    mission.aerodynamics.finsGeom.Lambda_le, mission.aerodynamics.finsGeom.tmac];


options = odeset('RelTol',1e-6,'AbsTol',1e-2 ,'Events',@(t,x) propEvent(t,x, mission,configuration.stage{stageNumber}.mProp,configuration,stageNumber));


solution = ode113(@(t,x) launcherDynamicsAndControlECI(t, x, mission, mer, configuration, launcher, stageNumber, guidancePoints, guidanceTime, gains, dimensions, engineVec, finsVec), tSpan, x0, options);

tt = linspace(solution.x(1),solution.x(end), nDeval);

xx = deval(solution,tt);

end