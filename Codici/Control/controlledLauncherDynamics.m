% function dxdt = controlledLauncherDynamics(t, x, mission, windDirection,targetPoint)
%     % -------------------------------------------------
%     %  State convention:
%     %  x(1:3)   = r_I  (position in the inertial frame)
%     %  x(4:6)   = v_I  (velocity in the inertial frame)
%     %  x(7:10)  = q    (quaternion from body to inertial frame [q0;q1;q2;q3])
%     %  x(11:13) = w_B  (angular velocity in the body frame)
%     % -------------------------------------------------
% 
%     rI = x(1:3); % position of launcher
%     vI = x(4:6); % velocity of launcher
%     q  = x(7:10); % quaternion from body to inertial frame
%     wB = x(11:13); % angular velocity in body frame
% 
%     errPos = ;
%     errVel = ;
%     errAngl = ;
%     errAnglVel = ;
% 
%     forceActuator = 0;
% 
% 
%     % quaternionNormalisation
%     q = q / norm(q);
%     GM     = mission.environment.GM;
%     rEarth = mission.environment.rEarth;
% 
%     % altitude&density
% 
%     h   = norm(rI) - rEarth;              % [m]
%     rho = mission.environment.rhoFun(h);
% 
%     % ------------ Wind------------------
%     windIntensity = 0;
%     windVelocityI = vI - windIntensity * windDirection;
% 
%     % cd and area informations (to note that this might change wrt attitude)
%     Cd   = mission.capsule.supersonicCD;
%     Aref = mission.capsule.Area;
% 
%     vRelNorm = norm(windVelocityI);
%     if vRelNorm > 0
%         aeroForceI = -0.5 * rho * vRelNorm * windVelocityI * Cd * Aref;
%     else
%         aeroForceI = [0;0;0];
%     end
%     %aeroForceI = [0; 0; 0];
%     % ------------ First cardinal equation ------------
%     m = mission.capsule.weigth;  
% 
%     % gravity acceleration
%     rNorm = norm(rI);
%     aGravI = -GM * rI / (rNorm^3);
% 
%     % Accelerazione totale
%     aI = aGravI + aeroForceI / m + Fcmd/m;
% 
%     % ------------ Quaternion Rotation Matrix ------------
%     q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
% 
%     % ------------ quaternion evolution via kinematic ------------
%     wx = wB(1); wy = wB(2); wz = wB(3);
%     Omega = [  0   wz  -wy  wx;
%                -wz    0   wx  wy;
%                wy  -wx    0   wz;
%                -wx   -wy  -wz    0 ];
% 
%     qdot = 0.5 * Omega * q;
% 
%     % Rbi: body -> inertial
%     Rbi = [q0^2-q1^2-q2^2+q3^2,  2*(q0*q1+q3*q2),  2*(q0*q2-q1*q3);...
%            2*(q0*q1-q3*q2),  -q0^2+q1^2-q2^2+q3^2,  2*(q1*q2+q0*q3);...
%            2*(q0*q2+q1*q3),  2 * (q1*q2-q0*q3), -q0^2-q1^2+q2^2+q3^2]';
% 
%     targetVec =windVelocityI/norm(windVelocityI);
% 
%     % Rib: inerziale -> body
%     Rib = Rbi';
% 
%     ySystem = Rib(2,:);
%     zSystem = Rib(3,:);
% 
%     alphaY = zSystem*targetVec;
%     alphaZ = -ySystem*targetVec;
%     % AeroForce In body reference frame
%     aeroForceB = Rib * aeroForceI;
% 
%     % ------------ aerodynamic torque wrt center of mass (2ª cardinale) ------------
%     % cp e cg defined in body frame (euler equations)
%     vBody = Rib * windVelocityI;
%     alpha = acos(dot([1;0;0], vBody)/norm(vBody) );
%     xCgTop = 20; 
%     cpB = [( 0.5* 25 *sin(alpha)^2 ).*(alpha<=pi/2) + (25-25/2*sin(alpha)^2).*(alpha>pi/2) - xCgTop;0; 0];   
%     %scatter(t,alpha)
%     %hold on
%     cgB = [0; 0; 0];
%     % scatter(cgB,0)
%     % hold on
%     % scatter(cpB,0)
% 
%     rCpCgB = cpB-cgB;   % torque arm
%     torque = cross(rCpCgB, aeroForceB) -50*[wx;wy;wz] - 10*[0; alphaY; alphaZ] + Mcmd;   % M = r x F  (nel body)
% 
%     % Matrice di inerzia nel body (3x3)
%     Ix = 1000;
%     Iy = 5000;
%     Iz = 5000;
%     J = diag([Ix,Iy Iz]);
% 
%     % 2ª cardinal: J * wdot + w x (J*w) = M
%     Jw   = J * wB;
%     wdot = J \ (torque - cross(wB, Jw));
% 
%     % ------------ build dxdt ------------
%     dxdt         = zeros(13,1);
%     dxdt(1:3)    = vI;    % dr/dt = v
%     dxdt(4:6)    = aI;    % dv/dt
%     dxdt(7:10)   = qdot;  % dq/dt
%     dxdt(11:13)  = wdot;  % dw/dt
% end


function dxdt = controlledLauncherDynamics(t, x, mission, windDirection, targetPoint)

    %% ---------------- STATE UNPACKING ----------------
    rI = x(1:3); 
    vI = x(4:6);
    q  = x(7:10);   % quaternion body→inertial
    wB = x(11:13);

    q = q / norm(q);

    GM     = mission.environment.GM;
    rEarth = mission.environment.rEarth;
    m      = mission.capsule.weigth;

    %% ---------------- POSITION / VELOCITY ERRORS ----------------
    % You want vertical control only → target is a point
    err_pos = ;            % 3×1
    err_vel = ;                         % desire zero velocity at touchdown

    %% ---------------- ATTITUDE TARGET ----------------
    % Desired body z-axis is aligned with the relative wind
  

    R_d = [ex_d ey_d ez_d];              % desired orientation matrix

    %% quaternion from body→inertial actual: q
    Rbi = quatToRot(q);

    %% quaternion error
    q_d = rotToQuat(R_d);
    q_err = quatMultiply(q_d, quatConjugate(q));

    % small-angle attitude error vector
    err_theta = q_err(2:4) * sign(q_err(1));

    % angular velocity error (want zero)
    err_omega = -wB;

    %% ---------------- CALL CONTROLLER ----------------
    gains = mission.control.gains;
    [Fcmd, Mcmd] = controller(err_pos, err_vel, err_theta, err_omega, ...
                              gains.KPos, gains.KVel, gains.KpAtt, gains.KdAtt, m, [0;0;-9.81]);

    %% ---------------- GRAVITY ----------------
    rNorm = norm(rI);
    aGravI = -GM * rI / (rNorm^3);

    %% ---------------- AERODYNAMICS ----------------
    h   = rNorm - rEarth;
    rho = mission.environment.rhoFun(h);

    Cd   = mission.capsule.supersonicCD;
    Aref = mission.capsule.Area;

    vRel = windVelocityI;
    vRelNorm = norm(vRel);

    if vRelNorm > 1e-6
        aeroForceI = -0.5 * rho * Cd * Aref * vRelNorm * vRel;
    else
        aeroForceI = [0;0;0];
    end

    %% ---------------- ROTATION MATRICES ----------------
    Rib = Rbi';

    %% ---------------- FORCE RESULTANT ----------------
    aI = aGravI + aeroForceI/m + Fcmd/m;

    %% ---------------- QUATERNION KINEMATICS ----------------
    qdot = 0.5 * [ 0      -wB';
                   wB  -skew(wB) ] * q;

    %% ---------------- AERODYNAMIC TORQUE + CONTROL TORQUE ----------------
    aeroForceB = Rib * aeroForceI;
    torqueB = Mcmd + cross([0;0;2], aeroForceB);   % if cp is 2m ahead

    %% ---------------- ROTATIONAL DYNAMICS ----------------
    J = diag([1000, 5000, 5000]);
    wdot = J \ (torqueB - cross(wB, J*wB));

    %% ---------------- BUILD OUTPUT ----------------
    dxdt = zeros(13,1);
    dxdt(1:3)   = vI;
    dxdt(4:6)   = aI;
    dxdt(7:10)  = qdot;
    dxdt(11:13) = wdot;

end
