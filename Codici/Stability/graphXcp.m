%% === INPUT DATA ===



% Structure
d = 2*ones(size(time));   % Reference diameter of launcher wrt time [m]
h  = 2;   % Cone length [m]
Lb = 1.2*ones(size(time));    % Vehicle length wrt time [m]


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



%% %% PLOT XCG

% Parameters 
lc1 = 10;      % CG of stage 1 (m)
lc2 = 5;     % CG of stage 2 (m)
lco = 1;      % CG of nose cone (m)
mc1 = 1000;   % Mass of stage 1 (kg)
mc2 = 500;    % Mass of stage 2 (kg)
mco = 100;    % Mass of nose cone (kg)
m_dot1 = 50;  % Propellant burn rate stage 1 (kg/s)
m_dot2 = 20;  % Propellant burn rate stage 2 (kg/s)
t1 = 15;      % Stage 1 separation time (s)
t2 = 35;      % Stage 2 separation time (s)



% Preallocate CG array
Xcg = zeros(size(t));

% Compute CG for each time
for i = 1:length(t)
    Xcg(i) = computeXCG(lc1, lc2, lco, mc1, mc2, mco, m_dot1, m_dot2, t(i), t1, t2);
end

% Plot
figure;
plot(t, Xcg, 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Center of Gravity (m)');
title('Evolution of Center of Gravity over Time');
ylim([0 max(Xcg)*1.1]);


%% PLOT XCG AND XCP 

% Time vector

time = (0:100:2000)';     % Time steps 
N = 

% Parameters XCG 
lc1 = 10;      % CG of stage 1 (m)
lc2 = 5;     % CG of stage 2 (m)
lco = 1;      % CG of nose cone (m)
m = 


% Parameters XCP
d = 2*ones(size(time));   % Reference diameter of launcher wrt time [m]
h  = 2;   % Cone length [m]
Lb = 1.2*ones(size(time));    % Vehicle length wrt time [m]


% Flight
alpha = 5 * ones(size(time));   % Angle of attack wrt time [deg]

