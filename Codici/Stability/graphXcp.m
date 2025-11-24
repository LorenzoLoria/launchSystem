%% === INPUT DATA ===

time = (0:100:2000)';     % Time steps 

% Launcher data

% Structure
d = 2*ones(size(time);   % Reference diameter of launcher wrt time [m]
h  = 0.15;   % Cone length [m]
Lb = 1.2*ones(size(time));    % Vehicle length wrt time [m]
Xcg = 1.2*ones(size(time));   % Center of gravity wrt time [m]

% Flight
alpha = 5 * ones(size(time));   % Angle of attack wrt time [deg]
istage = 9;   % Time of first staging

% Data with flare 

hf = 2;   % Flare height [m]
db = 3;   % Flare diameter [m]

%% === COMPUTE Xcp with no flare and no fins for each time step ===

Xcp = zeros(size(time));

for i = 1:length(time)
    Xcp(i) = computeXcp(alpha(i), Lb(i), h);
end


%% === COMPUTE Xcp with flare for each time step ===

Xcp_with_flare = zeros(size(time));

for i = 1:length(time)
    if i<istage
        Xcp_with_flare(i) = computeFlareXcp(Lb(i), d(i), h, hf, db(i));
    else
        Xcp_with_flare(i) = 2/3*h;
end


%% === COMPUTE Xcp with fins for each time step ===


%% === PLOT ===

figure;
plot(time, Xcp, 'LineWidth', 2); hold on;
plot(time, Xcp_with_flare, 'LineWidth', 2);
grid on;

xlabel('Time [s]');
ylabel('X_{cp} [m]');
title('Center of Pressure X_{cp} vs Time');
legend('No flare', 'With flare');