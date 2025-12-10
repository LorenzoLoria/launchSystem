function controlPlots(t, x, thrustData, mission, configuration, launcher, stageNumber, guidancePoints)
% =========================================================================
%   Plots position, velocity, attitude errors + thrust magnitude & gimbal
%   for FREE trajectory vs CONTROLLED trajectory.
% =========================================================================

tspan = t(:);

N = length(tspan);

% Allocate logs
err_pos_free  = zeros(N,2);
err_pos_ctrl  = zeros(N,2);

err_vel_free  = zeros(N,2);
err_vel_ctrl  = zeros(N,2);

attErr_free   = zeros(N,1);
attErr_ctrl   = zeros(N,1);

Flog          = zeros(N,1);
GimbalLog     = zeros(N,1);


for k = 1:N

    % -------- Extract free-flight states --------
    r   = x(k,1:2).';
    v   = x(k,3:4).';
    th  = x(k,5);
    om  = x(k,6);
    m   = x(k,7);


    % -------- Desired values from guidance --------
    r_des = guidancePoints(1:2);
    v_des = guidancePoints(4:5);    
    th_des = atan2(a_des(2), a_des(1));

    % ==========================================================
    %   ERRORS
    % ==========================================================

    err_pos(k,:) = (r_des - r).';
    err_vel(k,:) = (v_des - v).';
    attErr(k)  = angleDiff(th_des, th_);


    % ==========================================================
    %   THRUST + GIMBAL (ONLY CONTROLLED CASE)
    % ==========================================================
    opt = thrustData(tspan(k));

    throttling  = opt(1);
    gimbal      = deg2rad(opt(2));
    nEng        = configuration.stage{stageNumber}.nEngines;

    % Determine thrust level
    if stageNumber == 1
        staticContr  = (101325 - mission.environment.pressure(r_c)) ...
                           * configuration.stage{stageNumber}.engine.effAreaZero;
        nominalTh    = configuration.stage{stageNumber}.engine.thrustZero;
    else
        staticContr  = 0;
        nominalTh    = configuration.stage{stageNumber}.engine.thrustVacuum;
    end

    thrustMag = throttling * nEng * (nominalTh + staticContr);

    Flog(k)      = thrustMag;
    GimbalLog(k) = gimbal;

end

%% ===========================
%    POSITION ERROR
% =============================

figure;
plot(tspan, vecnorm(err_pos,2,2), 'r', 'LineWidth',1.5); hold on;
xlabel('Time [s]'); ylabel('|| position error || [m]');
title('Position Error vs Time');
grid on;


%% ===========================
%    VELOCITY ERROR
% =============================

figure;
plot(tspan, vecnorm(err_vel,2,2), 'r', 'LineWidth',1.5); hold on;
xlabel('Time [s]'); ylabel('|| velocity error || [m/s]');
title('Velocity Error vs Time');
grid on;


%% ===========================
%    ATTITUDE ERROR
% =============================

figure;
plot(tspan, abs(attErr), 'r', 'LineWidth',1.5); hold on;
xlabel('Time [s]'); ylabel('| attitude error | [rad]');
title('Attitude Error vs Time'); 
grid on;


%% ===========================
%    THRUST MAGNITUDE + GIMBAL
% =============================

figure;
subplot(2,1,1);
plot(tspan, Flog, 'k','LineWidth',1.5);
xlabel('Time [s]'); ylabel('Thrust [N]');
title('Thrust Magnitude'); grid on;

subplot(2,1,2);
plot(tspan, GimbalLog * 180/pi, 'm','LineWidth',1.5);
xlabel('Time [s]'); ylabel('Gimbal [deg]');
title('Gimbal Angle'); grid on;



end
