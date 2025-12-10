function controlPlots(thrustData, configuration, stageNumber, timeCollocation, stateCollocation, y0)

% =========================================================================
%   Plots position, velocity, attitude errors + thrust magnitude & gimbal
%   for FREE trajectory vs CONTROLLED trajectory.
% =========================================================================

tspan = timeCollocation(:,1);   
guidancePoints = stateCollocation(1:6,1,stageNumber); 
launcher = [2,1,4,4,0.4056,0.4016,0.7];

[mission,settings] = dataStructGlobal;

% Gains for PD guidance
gains.Kp_pos = 0.01;
gains.Kd_vel = 0.01;
gains.Kp_theta = 0.1;
gains.Kp_omega = 0.05;


% Integration of Dynamics


options = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [t, x] = ode113(@(t,x) launcherDynamicsAndControlECI(t, x, thrustData, mission, configuration, launcher, stageNumber, guidancePoints, gains), ...
                    tspan, y0, options);



N = length(t);

% Allocate logs
err_pos  = zeros(N,2);
err_vel  = zeros(N,2);
attErr  = zeros(N,1);
Flog          = zeros(N,1);
GimbalLog     = zeros(N,1);


for k = 1:N

    % -------- Extract free-flight states --------
    r = x(k,1:3).';
    v = x(k,4:6).';
    m = x(k,7);
    th = zeros(N,1);


    % -------- Desired values from guidance --------
    r_des = guidancePoints(1:2);
    v_des = guidancePoints(4:5);    
    th_des = atan2(v_des(2), v_des(1));


    % ==========================================================
    %   ERRORS
    % ==========================================================

    err_pos(k,:) = (r_des - r).';
    err_vel(k,:) = (v_des - v).';
    attErr(k)  = angleDiff(th_des, th);


    % ==========================================================
    %   THRUST + GIMBAL (ONLY CONTROLLED CASE)
    % ==========================================================
    opt = thrustData(k,:);
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