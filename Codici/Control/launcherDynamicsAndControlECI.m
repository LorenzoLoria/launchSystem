function dxdt = launcherDynamicsAndControlECI(t, x, mission, mer, configuration, launcher, stageNumber, guidancePoints, guidanceTime, gains)

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
        configuration.geometry.stage{ii}.length + ...
        configuration.geometry.stage{ii}.interstage.length;
end

% Limits
deltaGimbalMax = 20;   % [deg]
maxTiltAngle   = 20;   % [deg]

BRFtoIRF = [ 1-2*(q(3)^2+q(4)^2),       2*(q(2)*q(3) - q(4)*q(1)),   2*(q(2)*q(4) + q(3)*q(1));...
             2*(q(2)*q(3) + q(4)*q(1)), 1-2*(q(2)^2+q(4)^2),         2*(q(3)*q(4) - q(2)*q(1));...
             2*(q(2)*q(4) - q(3)*q(1)), 2*(q(3)*q(4) + q(2)*q(1)),   1-2*(q(2)^2+q(3)^2) ];...

IRFtoBRF = BRFtoIRF';   % inertial -> body

I = InertiaEvaluation(mission, configuration, mer, launcher, totalStageNumber);

% Desire Position evaluation
xDes = interp1(guidanceTime, guidancePoints(1,:), t, 'pchip', 'extrap');
yDes = interp1(guidanceTime, guidancePoints(2,:), t, 'pchip', 'extrap');
zDes = interp1(guidanceTime, guidancePoints(3,:), t, 'pchip', 'extrap');
rDes = [xDes; yDes; zDes];

vxDes = interp1(guidanceTime, guidancePoints(4,:), t, 'pchip', 'extrap');
vyDes = interp1(guidanceTime, guidancePoints(5,:), t, 'pchip', 'extrap');
vzDes = interp1(guidanceTime, guidancePoints(6,:), t, 'pchip', 'extrap');
vDes = [vxDes; vyDes; vzDes];

h   = rMag - mission.environment.rEarth;
rho = mission.environment.rhoFun(h);
dynamicPressure = 0.5 * rho * vMag^2;
soundspeed      = mission.aerodynamics.soundspeedFun(h);
Mach            = vMag / soundspeed ;

% Angle of attack: angle between body x-axis and velocity in BRF
vBRF  = IRFtoBRF * v;
alpha = atan2( sqrt(vBRF(2)^2 + vBRF(3)^2), vBRF(1) );

% Aerodynamic coefficients
if Mach == 0
    Cd = 0.01;
    Cl = 0.0;
else
    [Cl, Cd, ~, ~] = CLCDcomputation(Mach, alpha, dynamicPressure, 1, mission, stageNumber, configuration);
end

% Engine and thrust data
if stageNumber == 1
    staticContribution = (101325 - mission.environment.pressure(h)) * ...
                         configuration.stage{stageNumber}.engine.effAreaZero;
    nominalThrust = configuration.stage{stageNumber}.engine.thrustZero;
    iSp           = configuration.stage{stageNumber}.engine.ispZero;
else
    staticContribution = 0;
    nominalThrust = configuration.stage{stageNumber}.engine.thrustVacuum;
    iSp           = configuration.stage{stageNumber}.engine.ispVacuum;
end

nEngines = configuration.stage{stageNumber}.nEngines;

throttling = 1.0;

TMax = throttling * nEngines * (nominalThrust + staticContribution);

% Center of Gravity
xCG = computeXcgGlobal(mission, configuration, launcher, mer, m);       % [m]

% Center of Pressure
xCP = computeXcpGlobal(mission, configuration, launcher, alpha);                  % [m]
xCpCg = xCP - xCG;
thrustArm = launcherLength - xCG;

% Drag Evaluation
if vMag > 1e-6
    dragIRF = -0.5 * rho * vMag * A * Cd * v;
else
    dragIRF = [0; 0; 0];
end

% Simple lift model: perpendicular to v and some reference axis (here z)
if vMag > 1e-6 
    liftDir = cross(v, [0;0;1]);
    if norm(liftDir) > 1e-6 
        liftDir = liftDir / norm(liftDir);
    else
        liftDir = [0;0;0];
    end
    liftMag = 0.5 * rho * vMag^2 * A * Cl;
    liftIRF = liftMag * liftDir;
else
    liftIRF = [0;0;0];
end

% Gravity
gravityIRF = -GM * r / rMag^3;

% Guidance Law
aCommandIRF = gains(1:3) .* (rDes - r) + ...
              gains(4:6) .* (vDes - v) - gravityIRF;

%  ATTITUDE: align body x-axis with vDes direction
if norm(vDes) < 1e-6
    vRef = aCommandIRF;
else
    vRef = vDes;
end

if norm(vRef) < 1e-8
    exDes = BRFtoIRF(:,1);      % fallback: keep current attitude
else
    exDes = vRef / norm(vRef);  % desired body x-axis in inertial frame
end

% Limit tilt angle wrt inertial Z
tiltAngle = acosd( exDes(3) / norm(exDes) );
if tiltAngle > maxTiltAngle
    exDes(3) = norm(exDes(1:2)) / tand(maxTiltAngle);
    exDes    = exDes / norm(exDes);
end

% Build desired rotation matrix 
yRef = [0; 1; 0];
if abs(dot(exDes, yRef)) > 0.99
    yRef = [0; 0; 1];
end
ezDes = cross(exDes, yRef); 
ezDes = ezDes / norm(ezDes);
eyDes = cross(ezDes, exDes);

Rdes = [exDes, eyDes, ezDes];

% Quaternion from desired DCM 
qDes = dcm2quat_shepperd(Rdes)';  

% Quaternion error: q_err = q_des ⊗ conj(q_curr)
qCurrConj = [q(1); -q(2); -q(3); -q(4)];

Q_des = [ qDes(1) -qDes(2) -qDes(3) -qDes(4);
          qDes(2)  qDes(1) -qDes(4)  qDes(3);
          qDes(3)  qDes(4)  qDes(1) -qDes(2);
          qDes(4) -qDes(3)  qDes(2)  qDes(1) ];

qErr = Q_des * qCurrConj;

% Short-rotation convention
if qErr(1) < 0
    qErr = -qErr;
end

% Error vector in body frame 
qErrVec = qErr(2:4);

% PD torque command (body frame)
torqueCommandBRF = gains(7:9) .* qErrVec + gains(10:12) .* (0 - omega);

% Longitudinal thrust to realise commanded acceleration along body x
aCommandBRF = IRFtoBRF * aCommandIRF;
thrustGuidanceBRF = m * aCommandBRF;   

% Command thrust from Guidance and Attitude
thrustBRF = [torqueCommandBRF(2) / thrustArm; -torqueCommandBRF(1) / thrustArm; thrustGuidanceBRF(3)];

% Gimbal limit: angle between thrust and body z (or x, depending on convention)
gimbalAngle = acosd( thrustBRF(3) / norm(thrustBRF));

if gimbalAngle > deltaGimbalMax
    % Rescale lateral components to max gimbal
    T_lat = norm(thrustBRF(1:2));
    Tz    = T_lat / tand(deltaGimbalMax);
    scale = sqrt(Tz^2 + T_lat^2) / max(norm(thrustBRF), 1e-6);
    thrustBRF = thrustBRF / scale;
end

% Saturate thrust magnitude
T_norm = norm(thrustBRF);
if T_norm > TMax && T_norm > 0
    thrustBRF = thrustBRF / T_norm * TMax;
end

% Thrust in inertial frame
thrustIRF = BRFtoIRF * thrustBRF;

% Rotational dynamics 
dragBRF = IRFtoBRF * dragIRF;
liftBRF = IRFtoBRF * liftIRF;

aeroTorqueBRF   = cross([0;0;xCpCg], dragBRF + liftBRF);
thrustTorqueBRF = cross([0;0;-thrustArm], thrustBRF);

totalTorqueBRF = thrustTorqueBRF + aeroTorqueBRF + torqueCommandBRF*0; % (if you want torqueCommandBRF only via thrust, leave *0)

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
dxdt(14)      = -norm(thrustBRF) / (g0 * iSp);

end