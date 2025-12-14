function dxdt = launcherDynamicsAndControlECI2(t, x, mission, mer, configuration, launcher, stageNumber, guidancePoints, guidanceTime, gains, dimensions, engineVec, finsVec,thrustData)

% STATE:
% x(1:3)   = r_ECI   [m]
% x(4:6)   = v_ECI   [m/s]
% x(7:10)  = q       [-]   quaternion (scalar-first, body->inertial)
% x(11:13) = omega   [rad/s] body angular rates in BRF
% x(14)    = m       [kg]
global angleErr
r     = x(1:3);
v     = x(4:6);
angle = x(7);
omega = x(8);
m     = x(9);

rMag = norm(r);
vMag = norm(v);

A   = mission.capsule.Area;
g0  = mission.environment.g0;
GM  = mission.environment.GM;

totalStageNumber = launcher(1);

% Approximate launcher length for thrust arm
launcherLength = mission.capsule.height;
for ii = totalStageNumber:-1:stageNumber
    launcherLength = launcherLength + ...
        configuration.geometry.stage{ii}.tanksLength + ...
        configuration.geometry.stage{ii}.interstage.length;
end

BRFtoIRF = mission.target.Rfinal' * [1 0 0 ; 0 0 -1 ; 0 1 0] * [cos(angle) 0 -sin(angle) ; 0 1 0 ; sin(angle) 0 cos(angle)];
IRFtoBRF = BRFtoIRF';

[I,xCG] = InertiaEvaluation(mission, configuration, mer, stageNumber, m, launcher);

% Desire Position evaluation
rDes = interp1(guidanceTime, guidancePoints(1:3,:)', t+5, 'linear', 'extrap');
vDes = interp1(guidanceTime, guidancePoints(4:6,:)', t+5, 'linear', 'extrap');

h   = rMag - mission.environment.rEarth;
rho = mission.environment.rhoFun(h);
dynamicPressure = 0.5 * rho * vMag^2;
soundspeed      = mission.aerodynamics.soundspeedFun(h);
Mach            = vMag / soundspeed ;

% Angle of attack: angle between body x-axis and velocity in BRF
vBRF  = IRFtoBRF * v;
alpha = atan2(vBRF(1),vBRF(3));
if norm(vBRF) == 0
    alpha = 0;
end
% Aerodynamic coefficients
if Mach == 0
    CdBody = 0.01;
    ClBody = 0.0;
else
    [~,~,~,~, ClBody, CdBody, ~, ~] = CLCDcomputation(Mach,alpha,dynamicPressure,1,mission,stageNumber,dimensions,engineVec, finsVec) ;

end

% Engine and thrust data
if stageNumber == 1
    staticContribution = (101325 - mission.environment.pressure(h)) * ...
                         configuration.stage{stageNumber}.engine.effAreaZero;
    nominalThrust = configuration.stage{stageNumber}.engine.thrustZero;
    iSp           = configuration.stage{stageNumber}.engine.ispZero;
else
    staticContribution = 0;
    nominalThrust = configuration.stage{stageNumber}.engine.thrustVacum;
    iSp           = configuration.stage{stageNumber}.engine.ispVac;
end

nEngines = configuration.stage{stageNumber}.nEngines;

% Center of Pressure
xCP = computeXcpGlobal(mission, configuration, launcher, alpha, stageNumber);        % [m]
xCpCg = xCP - xCG;
thrustArm = launcherLength - xCG;

% Drag Evaluation
if vMag > 1e-6
    dragIRF = -0.5 * rho * vMag * A * CdBody * v;
else
    dragIRF = [0; 0; 0];
end

% Simple lift model: perpendicular to v and some reference axis (here z)
if vMag > 1e-6 
    liftDir = cross([0;1;0],vBRF/norm(vBRF));
    if norm(liftDir) > 1e-6 
        liftDir = liftDir / norm(liftDir);
    else
        liftDir = [0;0;0];
    end
    liftMag = 0.5 * rho * vMag^2 * A * ClBody;
    liftBRF = liftMag * liftDir;
else
    liftBRF = [0;0;0];
end

% Gravity
gravityIRF = -GM * r / rMag^3;

if t<1
    thrust = thrustData(t);
    thrustBRF = thrust(1)* (nominalThrust + staticContribution) * configuration.stage{stageNumber}.nEngines * [0;0;1];
    thrustIRF = BRFtoIRF * thrustBRF;
    torqueCommandBRF = [0;0;0];
else

% --- ATTITUDE: align body z-axis with vDes direction in inertial frame ---

if norm(vDes)<1e-20
    refDirz = IRFtoBRF * (mission.launchBase.initialPointECI/norm(mission.launchBase.initialPointECI))';
else
    refDirz = vBRF/norm(vBRF);
end

refDiry = IRFtoBRF * mission.target.Rfinal(:,3);
refDirx = cross(refDiry,refDirz);


% Build desired rotation matrix R_des (columns = body axes in inertial frame)
% Choose a reference vector not almost parallel to ezDes

% exDes = cross(yRef, ezDes); 
% exDes = exDes / norm(exDes);  % body x-axis
% eyDes = cross(ezDes, exDes);  % body y-axis
% 
% % Ensure exDes, eyDes, and ezDes are column vectors
% if isrow(exDes)
%     exDes = exDes';
% end
% if isrow(eyDes)
%     eyDes = eyDes';
% end
% if isrow(ezDes)
%     ezDes = ezDes';
% end

Rdes = [refDirx, refDiry,refDirz ];  % in teroia la seconda colonna deve essere [0;1;0]
vRefBRF = IRFtoBRF * vDes';
if acos(dot([1;0;0],vRefBRF/norm(vRefBRF))) <= 0
    theta = acos(dot([0;0;1],vRefBRF/norm(vRefBRF)));
else
    theta = -acos(dot([0;0;1],vRefBRF/norm(vRefBRF)));
end


% Longitudinal thrust to realise commanded acceleration along body x
% Guidance Law
errPos = rDes' - r;
errVel = vDes' - v;

% PD torque command (body frame)
gainsStage = gains([1:8]+(stageNumber-1)*8);
% if acosd(dot([1;0;0],IRFtoBRF*errPos/norm(errPos))) < 0
%     angleErr =  acosd(dot([0;0;1],IRFtoBRF*errPos/norm(errPos)));
% else
%     angleErr = -acosd(dot([0;0;1],IRFtoBRF*errPos/norm(errPos)));
% end
targetPoint = guidancePoints(1:3,end)-x(1:3);
u = IRFtoBRF * (targetPoint / norm(targetPoint));   % direzione target in BRF

% errore di pitch nel piano x-z (firmato)
% se body-z = [0;0;1], l'angolo per puntare verso u è:
angleErr = atan2(u(1), u(3));    % radianti, firmato [-pi, pi]
torqueCommandBRF = (-100*gainsStage(7) * angleErr - gainsStage(8) * omega) .* [0;1;0] ;%+ gains(11:13) .* (0 - omega);

thrustGuidanceIRF = m .* (gainsStage(1:3)'.*errPos + gainsStage(4:6)'.*errVel);   


% Command thrust from Guidance and Attitude
attitudeForce = cross(torqueCommandBRF , (launcherLength-xCG).*[0;0;1]) ;

% Saturate thrust magnitude

TMax = nEngines * (nominalThrust + staticContribution);
T_norm = norm(thrustGuidanceIRF);

if T_norm > TMax && T_norm > 0
    thrustGuidanceIRF = thrustGuidanceIRF / T_norm * TMax;
end

if attitudeForce > 5e5
    attitudeForce = 5e5;
end
%torqueCommandBRF = cross(attitudeForce,(launcherLength-xCG).*[0;0;1]);

thrustIRF = thrustGuidanceIRF;
end

% scatter3(x(1),x(2),x(3),'b*')
% hold on
% scatter3(rDes(1),rDes(2),rDes(3),'r*')
% zIRF = BRFtoIRF * [0;0;1];
% xIRF = BRFtoIRF * [1;0;0];
% yIRF = BRFtoIRF * [0;1;0];
% quiver3(x(1),x(2),x(3), 1000000*zIRF(1), 1000000*zIRF(2), 1000000*zIRF(3), 0);
% plot3(guidancePoints(1,end),guidancePoints(2,end),guidancePoints(3,end),'go')
% axis equal
% hold off
% drawnow



thrustBRF = IRFtoBRF*thrustIRF;
if acosd(dot(thrustBRF/norm(thrustBRF),[1;0;0])) <= 0
    thrustAngle = acosd(dot(thrustBRF/norm(thrustBRF),[0;0;1]));
else
    thrustAngle = -acosd(dot(thrustBRF/norm(thrustBRF),[0;0;1]));
end
% 

% if thrustAngle > 30 && thrustAngle>=0
%     thrustAngle = 30;
%     thrustBRF = norm(thrustBRF) * [sind(thrustAngle) ; 0 ; cosd(thrustAngle)];
%     thrustIRF = BRFtoIRF * thrustBRF;
% elseif thrustAngle < -30 && thrustAngle<0
%     thrustAngle = -30;
%     thrustBRF = norm(thrustBRF) * [sind(thrustAngle) ; 0 ; cosd(thrustAngle)];
%     thrustIRF = BRFtoIRF * thrustBRF;
% end
% Rotational dynamics 
dragBRF = IRFtoBRF * dragIRF;
liftIRF = BRFtoIRF * liftBRF;

aeroTorqueBRF   = cross( dragBRF + liftBRF, [0;0;xCpCg]);
thrustTorqueBRF = cross([0;0;-thrustArm], [1;0;0]);

%torqueCommandBRF = torqueCommandBRF - aeroTorqueBRF;
if norm(torqueCommandBRF-aeroTorqueBRF) > 1e7
    torqueCommandBRF = ( 1e7 + norm(aeroTorqueBRF) ) .* torqueCommandBRF/norm(torqueCommandBRF);
end
totalTorqueBRF = torqueCommandBRF; 

omegaDot = totalTorqueBRF(2)/I(2,2);

% Assemble dxdt 
dxdt          = zeros(9,1);

dxdt(1:3)     = v;
dxdt(4:6)     = (thrustIRF + dragIRF + liftIRF) / m + gravityIRF;
dxdt(7)    = omega;
dxdt(8)   = omegaDot;
dxdt(9)      = - norm(thrustIRF) / (g0 * iSp);



% hold on
% scatter3(x(1)-mission.launchBase.initialPointECI(1),x(2)-mission.launchBase.initialPointECI(2),x(3)-mission.launchBase.initialPointECI(3),'b*')

end