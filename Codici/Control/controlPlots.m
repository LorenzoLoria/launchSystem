%% SIMULATION WITH AND WITHOUT CONTROL


tspan = linspace(0,200,2001);

%% PLOTS


% --- position error ---s
err_pos_free = ;
err_pos_ctrl = ;

% --- velocity error ---
err_vel_free = ;
err_vel_ctrl = ;

% --- attitude error ---
attErr_free = ;
attErr_ctrl = ;




%%  PLOTS: POSITION ERROR

figure;
plot(tspan, vecnorm(err_pos_free,2,2), 'r', 'LineWidth',1.5); hold on;
plot(tspan, vecnorm(err_pos_ctrl,2,2), 'b', 'LineWidth',1.5);
xlabel('Time [s]'); ylabel('||position error|| [m]');
title('Position Error vs Time');
legend('Without Control','With Control');
grid on;


%%  PLOTS: VELOCITY ERROR

figure;
plot(tspan, vecnorm(err_vel_free,2,2), 'r', 'LineWidth',1.5); hold on;
plot(tspan, vecnorm(err_vel_ctrl,2,2), 'b', 'LineWidth',1.5);
xlabel('Time [s]'); ylabel('||velocity error|| [m/s]');
title('Velocity Error vs Time');
legend('Without Control','With Control');
grid on;


%%  PLOTS: ATTITUDE ERROR

figure;
plot(tspan, vecnorm(attErr_free,2,2), 'r', 'LineWidth',1.5); hold on;
plot(tspan, vecnorm(attErr_ctrl,2,2), 'b', 'LineWidth',1.5);
xlabel('Time [s]'); ylabel('||attitude error vector||');
title('Attitude Error vs Time');
legend('Without Control','With Control');
grid on;


%%  PLOTS: THRUST COMMAND & GIMBAL ANGLE

Fmag = vecnorm(Flog,2,2);
gimbal_deg = GimbalLog * 180/pi;

figure;
subplot(2,1,1);
plot(tspan, Fmag, 'k','LineWidth',1.5);
xlabel('Time [s]'); ylabel('Thrust Magnitude [N]');
title('Thrust Command');
grid on;

subplot(2,1,2);
plot(tspan, gimbal_deg, 'm','LineWidth',1.5);
xlabel('Time [s]'); ylabel('Gimbal Angle [deg]');
title('Gimbal Angle');
grid on;