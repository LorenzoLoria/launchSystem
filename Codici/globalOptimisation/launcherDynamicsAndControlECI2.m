function dxdt = launcherDynamicsAndControlECI2(t, x, mission, mer, configuration, launcher, stageNumber, guidancePoints, guidanceTime, gains, dimensions, engineVec, finsVec,thrustData)

% STATE:
% x(1:3)   = r_ECI   [m]
% x(4:6)   = v_ECI   [m/s]
% x(7:10)  = q       [-]   quaternion (scalar-first, body->inertial)
% x(11:13) = omega   [rad/s] body angular rates in BRF
% x(14)    = m       [kg]

r     = x(1:3);
v     = x(4:6);
q     = x(7:10);
omega = x(11:13);
m     = x(14);

% Always keep quaternion normalised (cheap, helps numerics)
q     = q / norm(q);

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

BRFtoIRF = mission.target.Rfinal * [ 1-2*(q(3)^2+q(4)^2),       2*(q(2)*q(3) - q(4)*q(1)),   2*(q(2)*q(4) + q(3)*q(1));...
             2*(q(2)*q(3) + q(4)*q(1)), 1-2*(q(2)^2+q(4)^2),         2*(q(3)*q(4) - q(2)*q(1));...
             2*(q(2)*q(4) - q(3)*q(1)), 2*(q(3)*q(4) + q(2)*q(1)),   1-2*(q(2)^2+q(3)^2) ];...

IRFtoBRF = BRFtoIRF';   % inertial -> body

[I,xCG] = InertiaEvaluation(mission, configuration, mer, stageNumber, m, launcher);

% Desire Position evaluation
rDes = interp1(guidanceTime, guidancePoints(1:3,:)', t, 'linear', 'extrap');
vDes = interp1(guidanceTime, guidancePoints(4:6,:)', t, 'linear', 'extrap');

h   = rMag - mission.environment.rEarth;
rho = mission.environment.rhoFun(h);
dynamicPressure = 0.5 * rho * vMag^2;
soundspeed      = mission.aerodynamics.soundspeedFun(h);
Mach            = vMag / soundspeed ;

% Angle of attack: angle between body x-axis and velocity in BRF
vBRF  = IRFtoBRF * v;
alpha = -atan2(norm(vBRF(1:2)), vBRF(3));

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
    liftDir = cross([0;1;0],vBRF);
    if norm(liftDir) > 1e-6 
        liftDir = liftDir / norm(liftDir);
    else
        liftDir = [0;0;0];
    end
    liftMag = 0.5 * rho * vMag^2 * A * ClBody;
    liftIRF = liftMag * liftDir;
else
    liftIRF = [0;0;0];
end

% Gravity
gravityIRF = -GM * r / rMag^3;

if t<1
    thrust = thrustData(t);
    thrustBRF = thrust(1)* (nominalThrust + staticContribution) * configuration.stage{stageNumber}.nEngines * [0;0;1];
    thrustIRF = BRFtoIRF * thrustBRF;
    torqueCommandBRF = 0;
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
theta = acos(dot([0;0;1],refDirz));
% % Quaternion from desired DCM 
% qDes = dcm2quat_shepperd(Rdes)';  
% 
% % Quaternion error: q_err = q_des ⊗ conj(q_curr)
% qCurrConj = [q(1); -q(2); -q(3); -q(4)];
% 
% Q_des = [ qDes(1) -qDes(2) -qDes(3) -qDes(4);
%           qDes(2)  qDes(1) -qDes(4)  qDes(3);
%           qDes(3)  qDes(4)  qDes(1) -qDes(2);
%           qDes(4) -qDes(3)  qDes(2)  qDes(1) ];
% 
% qErrBRF = Q_des * qCurrConj;

% Short-rotation convention
% if qErrIRF(1) < 0
%     qErrIRF = -qErrIRF;
% end

% PD torque command (body frame)
torqueCommandBRF = gains(7) * theta * [0;1;0] ;%+ gains(11:13) .* (0 - omega);

% Longitudinal thrust to realise commanded acceleration along body x
% Guidance Law
errPos = rDes' - r;
errVel = vDes' - v;

thrustGuidanceIRF = m .* (gains(1:3)'.*errPos + gains(4:6)'.*errVel);   


% Command thrust from Guidance and Attitude
attitudeForce = cross(torqueCommandBRF , (launcherLength-xCG).*[0;0;1]) ;

% Saturate thrust magnitude

TMax = nEngines * (nominalThrust + staticContribution);
T_norm = norm(thrustGuidanceIRF);
if T_norm > TMax && T_norm > 0
    thrustGuidanceIRF = thrustGuidanceIRF / T_norm * TMax;
end

if attitudeForce > 1e5
    attitudeForce = 1e5;
end
torqueCommandBRF = cross((launcherLength-xCG).*[0;0;1],attitudeForce);

thrustIRF = thrustGuidanceIRF;
end
% Rotational dynamics 
dragBRF = IRFtoBRF * dragIRF;
liftBRF = IRFtoBRF * liftIRF;

aeroTorqueBRF   = cross([0;0;xCpCg], dragBRF + liftBRF);
thrustTorqueBRF = cross([0;0;-thrustArm], [1;0;0]);

totalTorqueBRF = thrustTorqueBRF*0 + aeroTorqueBRF + torqueCommandBRF; % (if you want torqueCommandBRF only via thrust, leave *0)

omegaDot = I \ ( totalTorqueBRF - cross(omega, I * omega) );

% Quaternion kinematics 
OmegaMat = [ 0        -omega(1) -omega(2) -omega(3);
             omega(1)  0         omega(3) -omega(2);
             omega(2) -omega(3)  0         omega(1);
             omega(3)  omega(2) -omega(1)  0        ];

% Assemble dxdt 
dxdt          = zeros(14,1);

dxdt(1:3)     = v;
dxdt(4:6)     = (thrustIRF + dragIRF + liftIRF) / m + gravityIRF;
dxdt(7:10)    = 0.5 * OmegaMat * q;
dxdt(11:13)   = omegaDot;
dxdt(14)      = - norm(thrustIRF) / (g0 * iSp);



% hold on
% scatter3(x(1)-mission.launchBase.initialPointECI(1),x(2)-mission.launchBase.initialPointECI(2),x(3)-mission.launchBase.initialPointECI(3),'b*')

end