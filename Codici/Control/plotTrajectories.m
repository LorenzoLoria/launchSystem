function plotTrajectories(thrustData, configuration, mer, guidancePoints, guidanceTime)
% =========================================================================
%   Plots POSITION and VELOCITY for Controlled, No-Control, and Ideal
% =========================================================================
launcher = [2,1,4,4,0.4056,0.4016,0.7];
stageNumber = launcher(1);

[mission,settings] = dataStructGlobal;
gains  = [10 10 10 10 10 10 10 10 10 10 10 10]';
gains0 = [50 0 10 1 0 60 100 10 4 2 8 8]';

% -------------------------------------------------------------------------
%   Compute trajectories
% -------------------------------------------------------------------------
[time_control, x_control] = totalTrajectoryControl( ...
        mission, mer, configuration, settings, launcher, ...
        thrustData, guidancePoints, guidanceTime, gains );

[time_nocontrol, x_nocontrol] = totalTrajectoryControl( ...
        mission, mer, configuration, settings, launcher, ...
        thrustData, guidancePoints, guidanceTime, gains0 );

% Use first column of collocation times for plotting
time_control = time_control(:);
time_nocontrol = time_nocontrol(:);

% Guidance trajectory
r_guidance = guidancePoints(1:2,:);  % position
v_guidance = guidancePoints(3:4,:);  % velocity

% -------------------------------------------------------------------------
%   PLOTS
% -------------------------------------------------------------------------

% ---------- POSITION ----------
figure;
plot(time_control, x_control(:,1), 'b', 'LineWidth',1.5); hold on;
plot(time_control, x_control(:,2), 'b--', 'LineWidth',1.5);

plot(time_nocontrol, x_nocontrol(:,1), 'r', 'LineWidth',1.5);
plot(time_nocontrol, x_nocontrol(:,2), 'r--', 'LineWidth',1.5);

plot(guidanceTime, r_guidance(1,:), 'k', 'LineWidth',1.5);
plot(guidanceTime, r_guidance(2,:), 'k--', 'LineWidth',1.5);

xlabel('Time [s]');
ylabel('Position [m]');
title('Position vs Time');
grid on;
legend('Controlled X','Controlled Y', 'No-Control X','No-Control Y', 'Ideal X','Ideal Y');

% ---------- VELOCITY ----------
figure;
plot(time_control, x_control(:,3), 'b', 'LineWidth',1.5); hold on;
plot(time_control, x_control(:,4), 'b--', 'LineWidth',1.5);

plot(time_nocontrol, x_nocontrol(:,3), 'r', 'LineWidth',1.5);
plot(time_nocontrol, x_nocontrol(:,4), 'r--', 'LineWidth',1.5);

plot(guidanceTime, v_guidance(1,:), 'k', 'LineWidth',1.5);
plot(guidanceTime, v_guidance(2,:), 'k--', 'LineWidth',1.5);

xlabel('Time [s]');
ylabel('Velocity [m/s]');
title('Velocity vs Time');
grid on;
legend('Controlled X','Controlled Y', 'No-Control X','No-Control Y', 'Ideal X','Ideal Y');

end
