%% === INPUT DATA ===

time = (0:100:2000)';     % Time steps 

alpha = 5 * ones(size(time));   % Angle of attack [deg]
Lb = 1.2;                       % Vehicle length [m]
h  = 0.15;                      % Cone length [m]

%% === COMPUTE Xcp for each time step ===

Xcp = zeros(size(time));

for i = 1:length(time)
    Xcp(i) = computeXcp(alpha(i), Lb, h);
end

%% === PLOT ===

figure;
plot(time, Xcp, 'LineWidth', 2);
grid on;

xlabel('Time [s]');
ylabel('X_{cp} [m]');
title('Center of Pressure X_{cp} vs Time');
