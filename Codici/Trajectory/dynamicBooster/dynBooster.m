function dxdt = dynBooster(~, x, mission, windDirection)

    % -------------------------------------------------
    %  State convention:
    %  x(1:3)   = r_I  (position in the inertial frame)
    %  x(4:6)   = v_I  (velocity in the inertial frame)
    %  x(7:10)  = q    (quaternion from body to inertial frame [q0;q1;q2;q3])
    %  x(11:13) = w_B  (angular velocity in the body frame)
    % -------------------------------------------------

    rI = x(1:3);
    vI = x(4:6);
    q  = x(7:10);
    wB = x(11:13);

    % quaternionNormalisation
    q = q / norm(q);

    GM     = mission.environment.GM;
    rEarth = mission.environment.rEarth;

    % altitude&density

    h   = norm(rI) - rEarth;              % [m]
    rho = mission.environment.gridInterp(h);

    % ------------ Wind------------------
    windIntensity = 0;
    windVelocityI = vI - windIntensity * windDirection;
    
    % cd and area informations (to note that this might change wrt attitude)
    Cd   = mission.capsule.supersonicCD;
    Aref = mission.capsule.Area;

    vRelNorm = norm(windVelocityI);
    if vRelNorm > 0
        aeroForceI = -0.5 * rho * vRelNorm * windVelocityI * Cd * Aref;
    else
        aeroForceI = [0;0;0];
    end

    % ------------ First cardinal equation ------------
    m = mission.capsule.weigth;  

    % gravity acceleration
    rNorm = norm(rI);
    aGravI = -GM * rI / (rNorm^3);

    % Accelerazione totale
    aI = aGravI + aeroForceI / m;

    % ------------ Quaternion Rotation Matrix ------------
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);

    % Rbi: body -> inertial

    Rbi = [ 1-2*(q2^2+q3^2),   2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2);
            2*(q1*q2 + q0*q3), 1-2*(q1^2+q3^2),   2*(q2*q3 - q0*q1);
            2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1-2*(q1^2+q2^2) ];

    % Rib: inerziale -> body
    Rib = Rbi.';

    % AeroForce In body reference frame
    aeroForceB = Rib * aeroForceI;

    % ------------ aerodynamic torque wrt center of mass (2ª cardinale) ------------
    % cp e cg defined in body frame (euler equations)
    
    cpB = [0; 0; -2];   
    cgB = [0; 0; 0];

    rCpCgB = cpB - cgB;   % torque arm
    torque = cross(rCpCgB, aeroForceB);   % M = r x F  (nel body)

    % Matrice di inerzia nel body (3x3)
    Ix = 10000;
    Iy = 10000;
    Iz = 1000;
    J = diag([Ix,Iy Iz]);

    % 2ª cardinal: J * wdot + w x (J*w) = M
    Jw   = J * wB;
    wdot = J \ (torque - cross(wB, Jw));

    % ------------ quaternion evolution via kinematic ------------
    wx = wB(1); wy = wB(2); wz = wB(3);
    Omega = [  0   -wx  -wy  -wz;
               wx    0   wz  -wy;
               wy  -wz    0   wx;
               wz   wy  -wx    0 ];

    qdot = 0.5 * Omega * q;

    % ------------ build dxdt ------------
    dxdt         = zeros(13,1);
    dxdt(1:3)    = vI;    % dr/dt = v
    dxdt(4:6)    = aI;    % dv/dt
    dxdt(7:10)   = qdot;  % dq/dt
    dxdt(11:13)  = wdot;  % dw/dt

end
