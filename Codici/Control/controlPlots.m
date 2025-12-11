function controlPlots(thrustData, configuration, mer, stageNumber)
% =========================================================================
%   Plots position, velocity, attitude errors + thrust magnitude & gimbal
%   for FREE trajectory vs CONTROLLED trajectory.
% =========================================================================


launcher = [2,1,4,4,0.4056,0.4016,0.7];
[mission,settings] = dataStructGlobal;


tspan = [0 500];



% control on

[timeCollocationControlled, stateCollocationControlled] = totalTrajectoryControl(mission,mer,configuration,settings, launcher,thrustDataVec, guidancePoints, guidanceTime, gains);
x_control = stateCollocationControlled;
% control off

gains = [0 0 0 0 0 0 0 0 0 0 0 0]';

[timeCollocationUnControlled, stateCollocationUnControlled] = totalTrajectoryControl(mission,mer,configuration,settings, launcher,thrustDataVec, guidancePoints, guidanceTime, gains);
x_nocontrol = stateCollocationUnControlled;

N = length(t);

% Allocate logs
err_pos  = zeros(N,2);
err_vel  = zeros(N,2);
attErr   = zeros(N,1);
Flog     = zeros(N,1);
GimbalLog = zeros(N,1);

for k = 1:N
    % -------- Extract free-flight states --------
    r     = x(1:3);
    v     = x(4:6);
    q     = x(7:10);
    omega = x(11:13);
    m     = x(14); % Attitude angle

    r_control     = x_control(1:3);
    v_control     = x_control(4:6);
    q_control     = x_control(7:10);
    omega_control = x_control(11:13);
    m_control     = x_control(14);

    r_nocontrol     = x_nocontrol(1:3);
    v_nocontrol     = x_nocontrol(4:6);
    q_nocontrol     = x_nocontrol(7:10);
    omega_nocontrol = x_nocontrol(11:13);
    m_nocontrol     = x_nocontrol(14);

    % -------- Desired values from guidance --------
    r_des = guidancePoints(1:2);
    v_des = guidancePoints(3:4);
    th_des = atan2(v_des(2), v_des(1));

    % ==========================================================
    %   ERRORS
    % ==========================================================
    err_pos(k,:) = (r_des - r).';
    err_vel(k,:) = (v_des - v).';
    attErr(k) = angleDiff(th_des, theta);

    % ==========================================================
    %   THRUST + GIMBAL (ONLY CONTROLLED CASE)
    % ==========================================================
    throttling = thrustData(k,1);
    gimbal = deg2rad(thrustData(k,2));
    nEng = configuration.stage{stageNumber}.nEngines;

    % Determine thrust level
    if stageNumber == 1
        % Use current position to calculate altitude
        rMag = norm(x(k,1:2));
        h = rMag - mission.environment.rEarth;
        staticContr = (101325 - mission.environment.pressure(h)) * ...
                      configuration.stage{stageNumber}.engine.effAreaZero;
        nominalTh = configuration.stage{stageNumber}.engine.thrustZero;
    else
        staticContr = 0;
        nominalTh = configuration.stage{stageNumber}.engine.thrustVacuum;
    end

    thrustMag = throttling * nEng * (nominalTh + staticContr);
    Flog(k) = thrustMag;
    GimbalLog(k) = gimbal;
end

% ===========================
%    POSITION ERROR
% ===========================
figure;
plot(t, vecnorm(err_pos,2,2), 'r', 'LineWidth',1.5);
xlabel('Time [s]'); ylabel('|| position error || [m]');
title('Position Error vs Time'); grid on;

% ===========================
%    VELOCITY ERROR
% ===========================
figure;
plot(t, vecnorm(err_vel,2,2), 'r', 'LineWidth',1.5);
xlabel('Time [s]'); ylabel('|| velocity error || [m/s]');
title('Velocity Error vs Time'); grid on;

% ===========================
%    ATTITUDE ERROR
% ===========================
figure;
plot(t, abs(attErr), 'r', 'LineWidth',1.5);
xlabel('Time [s]'); ylabel('| attitude error | [rad]');
title('Attitude Error vs Time'); grid on;

% ===========================
%    THRUST MAGNITUDE + GIMBAL
% ===========================
figure;
subplot(2,1,1);
plot(t, Flog, 'k','LineWidth',1.5);
xlabel('Time [s]'); ylabel('Thrust [N]');
title('Thrust Magnitude'); grid on;

subplot(2,1,2);
plot(t, GimbalLog*180/pi, 'm','LineWidth',1.5);
xlabel('Time [s]'); ylabel('Gimbal [deg]');
title('Gimbal Angle'); grid on;
end
