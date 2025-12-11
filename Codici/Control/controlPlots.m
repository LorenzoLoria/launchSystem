% function controlPlots(thrustData, configuration, mer, guidancePoints, guidanceTime)
% % =========================================================================
% %   Plots position, velocity, attitude errors + thrust magnitude & gimbal
% %   for FREE (no control) trajectory vs CONTROLLED trajectory.
% % =========================================================================
% launcher = [2,1,4,4,0.4056,0.4016,0.7];
% 
% stageNumber = launcher(1);
% 
% [mission,settings] = dataStructGlobal;
% gains = [10 10 10 10 10 10 10 10 10 10 10 10]' ;
% gains0 = [50 0 10 1 0 60 100 10 4 2 8 8]' ;
% 
% % -------------------------------------------------------------------------
% %     Compute CONTROLLED trajectory
% % -------------------------------------------------------------------------
% 
% 
% [time_control, x_control] = totalTrajectoryControl( ...
%         mission, mer, configuration, settings, launcher, ...
%         thrustData, guidancePoints, ...
%         guidanceTime, gains );
% 
% % -------------------------------------------------------------------------
% %     Compute NO-CONTROL trajectory
% % -------------------------------------------------------------------------
% 
% 
% [time_nocontrol, x_nocontrol] = totalTrajectoryControl( ...
%         mission, mer, configuration, settings, launcher, ...
%         thrustData, guidancePoints, ...
%         guidanceTime, gains0 );
% 
% 
% % -------------------------------------------------------------------------
% %     Check consistency
% % -------------------------------------------------------------------------
% if ~isequal(time_control, time_nocontrol)
%     error('Time vectors for controlled and uncontrolled sims do not match.');
% end
% 
% t = time_control;
% N = length(t);
% 
% % Allocate logs
% err_pos_control  = zeros(N,2);
% err_vel_control  = zeros(N,2);
% attErr_control   = zeros(N,1);
% 
% err_pos_nocontrol  = zeros(N,2);
% err_vel_nocontrol  = zeros(N,2);
% attErr_nocontrol   = zeros(N,1);
% 
% Flog     = zeros(N,1);
% GimbalLog = zeros(N,1);
% 
% 
% % -------------------------------------------------------------------------
% %     LOOP OVER TIME
% % -------------------------------------------------------------------------
% for k = 1:N
% 
%     % ===== Controlled states =================================================
%     r_control     = x_control(k,1:2).';
%     v_control     = x_control(k,3:4).';
%     q_control     = x_control(k,7:10).';  
% 
% 
%     % ===== Uncontrolled states ===============================================
%     r_nc = x_nocontrol(k,1:2).';
%     v_nc = x_nocontrol(k,3:4).';
%     q_nc = x_nocontrol(k,7:10).';
% 
%     % ===== Desired values from guidance =====================================
%     r_des = guidancePoints(1:2,k);
%     v_des = guidancePoints(3:4,k);
% 
%     % ===== Desired attitude = pointing along v_des ===========================
% 
%     th_des = atan2(v_des(2), v_des(1));
% 
%     q_des = [cos(th_des/2); 0; 0; sin(th_des/2)];
% 
%     % -------- Convert current attitude to equivalent yaw angle ---------------
%     %  q = [w x y z]
%     yaw_control = atan2( 2*(q_control(1)*q_control(4) + q_control(2)*q_control(3)), ...
%                          1 - 2*(q_control(3)^2 + q_control(4)^2) );
% 
%     yaw_nc      = atan2( 2*(q_nc(1)*q_nc(4) + q_nc(2)*q_nc(3)), ...
%                          1 - 2*(q_nc(3)^2 + q_nc(4)^2) );
% 
%     % Desired yaw:
%     yaw_des = th_des;
% 
%     % Attitude errors (scalar angle)
%     attErr_control(k)   = wrapToPi(yaw_des - yaw_control);
%     attErr_nocontrol(k) = wrapToPi(yaw_des - yaw_nc);
% 
% 
%     % ==========================================================
%     %   POSITION and VELOCITY ERRORS
%     % ==========================================================
%     err_pos_control(k,:) = (r_des - r_control).';
%     err_vel_control(k,:) = (v_des - v_control).';
% 
%     err_pos_nocontrol(k,:) = (r_des - r_nc).';
%     err_vel_nocontrol(k,:) = (v_des - v_nc).';
% 
% 
%     % ==========================================================
%     %   THRUST + GIMBAL
%     % ==========================================================
%     throttling = thrustData(k,1);
%     gimbal     = deg2rad(thrustData(k,2));
%     nEng       = configuration.stage{stageNumber}.nEngines;
% 
%     if stageNumber == 1
%         rMag = norm(r_control);
%         h = rMag - mission.environment.rEarth;
%         staticContr = (101325 - mission.environment.pressure(h)) * ...
%                       configuration.stage{stageNumber}.engine.effAreaZero;
%         nominalTh = configuration.stage{stageNumber}.engine.thrustZero;
%     else
%         staticContr = 0;
%         nominalTh = configuration.stage{stageNumber}.engine.thrustVacuum;
%     end
% 
%     thrustMag = throttling * nEng * (nominalTh + staticContr);
% 
%     Flog(k) = thrustMag;
%     GimbalLog(k) = gimbal;
% 
% end
% 
% 
% % =========================================================================
% %   PLOTS
% % =========================================================================
% 
% % ===========================
% %    POSITION ERROR
% % ===========================
% figure;
% plot(t, vecnorm(err_pos_control,2,2), 'b','LineWidth',1.5); hold on;
% plot(t, vecnorm(err_pos_nocontrol,2,2), 'r','LineWidth',1.5);
% xlabel('Time [s]'); ylabel('|| position error || [m]');
% title('Position Error vs Time'); grid on;
% legend('Controlled','No Control');
% 
% % ===========================
% %    VELOCITY ERROR
% % ===========================
% figure;
% plot(t, vecnorm(err_vel_control,2,2), 'b','LineWidth',1.5); hold on;
% plot(t, vecnorm(err_vel_nocontrol,2,2), 'r','LineWidth',1.5);
% xlabel('Time [s]'); ylabel('|| velocity error || [m/s]');
% title('Velocity Error vs Time'); grid on;
% legend('Controlled','No Control');
% 
% % ===========================
% %    ATTITUDE ERROR
% % ===========================
% figure;
% plot(t, abs(attErr_control), 'b','LineWidth',1.5); hold on;
% plot(t, abs(attErr_nocontrol),'r','LineWidth',1.5);
% xlabel('Time [s]'); ylabel('| attitude error | [rad]');
% title('Attitude Error vs Time'); grid on;
% legend('Controlled','No Control');
% 
% % ===========================
% %    THRUST + GIMBAL
% % ===========================
% figure;
% subplot(2,1,1);
% plot(t, Flog, 'k','LineWidth',1.5);
% xlabel('Time [s]'); ylabel('Thrust [N]');
% title('Thrust Magnitude'); grid on;
% 
% subplot(2,1,2);
% plot(t, GimbalLog*180/pi, 'm','LineWidth',1.5);
% xlabel('Time [s]'); ylabel('Gimbal [deg]');
% title('Gimbal Angle'); grid on;
% 
% end
% 
% 

function controlPlots(thrustData, configuration, mer, guidancePoints, guidanceTime)
% =========================================================================
%   Plots position, velocity, attitude errors + thrust magnitude & gimbal
%   for FREE (no control) trajectory vs CONTROLLED trajectory.
% =========================================================================
launcher = [2,1,4,4,0.4056,0.4016,0.7];

stageNumber = launcher(1);

[mission,settings] = dataStructGlobal;
gains  = [10 10 10 10 10 10 10 10 10 10 10 10]' ;
gains0 = [50 0 10 1 0 60 100 10 4 2 8 8]' ;

% -------------------------------------------------------------------------
%     Compute CONTROLLED trajectory
% -------------------------------------------------------------------------
[timeCollocation_control, x_control] = totalTrajectoryControl( ...
        mission, mer, configuration, settings, launcher, ...
        thrustData, guidancePoints, ...
        guidanceTime, gains );

% -------------------------------------------------------------------------
%     Compute NO-CONTROL trajectory
% -------------------------------------------------------------------------
[timeCollocation_nocontrol, x_nocontrol] = totalTrajectoryControl( ...
        mission, mer, configuration, settings, launcher, ...
        thrustData, guidancePoints, ...
        guidanceTime, gains0 );

time_control = timeCollocation_control(:);
time_nocontrol = timeCollocation_nocontrol(:);
% -------------------------------------------------------------------------
%     NEW: Build a common time vector & interpolate
% -------------------------------------------------------------------------
t_end = min(time_control(end), time_nocontrol(end));
Npts  = min(length(time_control), length(time_nocontrol));
t     = linspace(0, t_end, Npts)';

% Interpolate states
x_control    = interp1(time_control,    x_control,    t);
x_nocontrol  = interp1(time_nocontrol,  x_nocontrol,  t);

% Interpolate guidance points (2 pos + 2 vel)
gt = linspace(guidanceTime(1), guidanceTime(end), size(guidancePoints,2))';
guidance_interp = interp1(gt, guidancePoints.', t, 'linear', 'extrap')';

% Extract desired along interpolated grid
r_des_all = guidance_interp(1:2,:);
v_des_all = guidance_interp(3:4,:);

% -------------------------------------------------------------------------
%     Continue unchanged from here
% -------------------------------------------------------------------------
N = length(t);

% Allocate logs
err_pos_control   = zeros(N,2);
err_vel_control   = zeros(N,2);
attErr_control    = zeros(N,1);

err_pos_nocontrol = zeros(N,2);
err_vel_nocontrol = zeros(N,2);
attErr_nocontrol  = zeros(N,1);

Flog      = zeros(N,1);
GimbalLog = zeros(N,1);

% -------------------------------------------------------------------------
%     LOOP OVER TIME
% -------------------------------------------------------------------------
for k = 1:N

    % ===== Controlled states ===============================================
    r_control = x_control(k,1:2).';
    v_control = x_control(k,3:4).';
    q_control = x_control(k,7:10).';

    % ===== Uncontrolled states =============================================
    r_nc = x_nocontrol(k,1:2).';
    v_nc = x_nocontrol(k,3:4).';
    q_nc = x_nocontrol(k,7:10).';

    % ===== Desired values (interpolated) ===================================
    r_des = r_des_all(:,k);
    v_des = v_des_all(:,k);

    % ===== Desired yaw from velocity =======================================
    th_des = atan2(v_des(2), v_des(1));
    yaw_des = th_des;

    % Current yaw from quaternion
    yaw_control = atan2( 2*(q_control(1)*q_control(4) + q_control(2)*q_control(3)), ...
                         1 - 2*(q_control(3)^2 + q_control(4)^2) );

    yaw_nc      = atan2( 2*(q_nc(1)*q_nc(4) + q_nc(2)*q_nc(3)), ...
                         1 - 2*(q_nc(3)^2 + q_nc(4)^2) );

    % Attitude errors
    attErr_control(k)   = wrapToPi(yaw_des - yaw_control);
    attErr_nocontrol(k) = wrapToPi(yaw_des - yaw_nc);

    % ==========================================================
    %   POSITION and VELOCITY ERRORS
    % ==========================================================
    err_pos_control(k,:) = (r_des - r_control).';
    err_vel_control(k,:) = (v_des - v_control).';

    err_pos_nocontrol(k,:) = (r_des - r_nc).';
    err_vel_nocontrol(k,:) = (v_des - v_nc).';

    % ==========================================================
    %   THRUST + GIMBAL
    % ==========================================================
    throttling = thrustData(k,1);
    gimbal     = deg2rad(thrustData(k,2));
    nEng       = configuration.stage{stageNumber}.nEngines;

    if stageNumber == 1
        rMag = norm(r_control);
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

% =========================================================================
%   PLOTS (unchanged)
% =========================================================================

figure;
plot(t, vecnorm(err_pos_control,2,2), 'b','LineWidth',1.5); hold on;
plot(t, vecnorm(err_pos_nocontrol,2,2), 'r','LineWidth',1.5);
xlabel('Time [s]'); ylabel('|| position error || [m]');
title('Position Error vs Time'); grid on; legend('Controlled','No Control');

figure;
plot(t, vecnorm(err_vel_control,2,2), 'b','LineWidth',1.5); hold on;
plot(t, vecnorm(err_vel_nocontrol,2,2), 'r','LineWidth',1.5);
xlabel('Time [s]'); ylabel('|| velocity error || [m/s]');
title('Velocity Error vs Time'); grid on; legend('Controlled','No Control');

figure;
plot(t, abs(attErr_control), 'b','LineWidth',1.5); hold on;
plot(t, abs(attErr_nocontrol),'r','LineWidth',1.5);
xlabel('Time [s]'); ylabel('| attitude error | [rad]');
title('Attitude Error vs Time'); grid on; legend('Controlled','No Control');

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
