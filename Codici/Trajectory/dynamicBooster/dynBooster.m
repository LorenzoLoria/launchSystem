function dxdt = dynBooster(~, x, mission, windDirection)

    % -------------------------------------------------
    %  Convenzione stato:
    %  x(1:3)   = r_I  (posizione in inerziale)
    %  x(4:6)   = v_I  (velocità in inerziale)
    %  x(7:10)  = q    (quaternione body->inerziale [q0;q1;q2;q3])
    %  x(11:13) = w_B  (velocità angolare nel frame corpo)
    % -------------------------------------------------

    rI = x(1:3);
    vI = x(4:6);
    q  = x(7:10);
    wB = x(11:13);

    % Normalizza quaternione per sicurezza numerica
    q = q / norm(q);

    GM     = mission.environment.GM;
    rEarth = mission.environment.rEarth;

    % Altitudine e densità
    h   = norm(rI) - rEarth;              % [m]
    rho = mission.environment.gridInterp(h);

    % ------------ Wind------------------
    windIntensity = 0;
    windVelocityI = vI - windIntensity * windDirection;

    Cd   = mission.capsule.supersonicCD;
    Aref = mission.capsule.Area;

    vrel_norm = norm(windVelocityI);
    if vrel_norm > 0
        aeroForceI = -0.5 * rho * vrel_norm * windVelocityI * Cd * Aref;
    else
        aeroForceI = [0;0;0];
    end

    % ------------ DINAMICA TRASLAZIONALE (1ª cardinale) ------------
    m = mission.capsule.weigth;   % ATTENZIONE: usa la MASSA, non il peso!

    % Gravità newtoniana
    r_norm = norm(rI);
    aGravI = -GM * rI / (r_norm^3);

    % Accelerazione totale
    aI = aGravI + aeroForceI / m;

    % ------------ MATRICE DI ROTAZIONE DA QUATERNIONE ------------
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);

    % Rbi: body -> inerziale
    Rbi = [ 1-2*(q2^2+q3^2),   2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2);
            2*(q1*q2 + q0*q3), 1-2*(q1^2+q3^2),   2*(q2*q3 - q0*q1);
            2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1-2*(q1^2+q2^2) ];

    % Rib: inerziale -> body
    Rib = Rbi.';

    % Forza aero in frame corpo
    aeroForceB = Rib * aeroForceI;

    % ------------ MOMENTO AERODINAMICO (2ª cardinale) ------------
    % cp e cg definiti nel frame corpo
    % (metti questi nei parametri del missione, qui li definisco come nel tuo esempio)
    cpB = [0; 0; -2];   
    cgB = [0; 0; 0];

    r_cp_cg_B = cpB - cgB;   % vettore dal CG al CP nel body
    MB = cross(r_cp_cg_B, aeroForceB);   % M = r x F  (nel body)

    % Matrice di inerzia nel body (3x3)
    Ix = 10000;
    Iy = 10000;
    Iz = 1000;
    J = diag([Ix,Iy Iz]);

    % 2ª cardinale: J * wdot + w x (J*w) = M
    Jw   = J * wB;
    wdot = J \ (MB - cross(wB, Jw));

    % ------------ CINEMATICA DEL QUATERNIONE ------------
    wx = wB(1); wy = wB(2); wz = wB(3);
    Omega = [  0   -wx  -wy  -wz;
               wx    0   wz  -wy;
               wy  -wz    0   wx;
               wz   wy  -wx    0 ];

    qdot = 0.5 * Omega * q;

    % ------------ COSTRUZIONE dxdt ------------
    dxdt         = zeros(13,1);
    dxdt(1:3)    = vI;    % dr/dt = v
    dxdt(4:6)    = aI;    % dv/dt
    dxdt(7:10)   = qdot;  % dq/dt
    dxdt(11:13)  = wdot;  % dw/dt

end
