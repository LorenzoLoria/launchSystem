%% Clear environment
close all
clear all

%% PARAMETERS
lc1 = 20;  % Stage 1 length [m]
lc2 = 10;  % Stage 2 length [m]
lco = 3;   % Base length [m]

m0  = 210.4e3; % Initial mass [kg]
mc1 = 162.5e3; % Stage 1 mass [kg]
mc2 = 40.6e3;  % Stage 2 mass [kg]
mco = m0 - mc1 - mc2; % Payload / structure mass

d  = 3;           % Diameter [m]
hf = 0.806 * d;   % Flare height [m]
db = 1.094 * d;   % Base diameter [m]

flare_ratio = (db - d)/(2*hf); % Flare ratio

Sref = pi*(d/2)^2; % Reference area [m^2] (cylinder cross-section)

t = 0:0.5:50; % Time vector [s]

%% Preallocate curves
Xcg_curve             = zeros(size(t));
Xcp_cone_curve        = zeros(size(t));
Xcp_flare_curve       = zeros(size(t));
Xcp_fin_flare_curve   = zeros(size(t));
Xcp_fin_curve         = zeros(size(t));

%% DYNAMIC PARAMETERS
% Stage number: 2 = stage 1, 1 = stage 2, 0 = no stage
N = ones(size(t));
N(t <= 20) = 2;
N(t > 20 & t <= 40) = 1;
N(t > 40) = 0;

% Mass profile
m = m0 - (m0-mc2-mco)/20 .* min(t,20) ...
      - (mc2/20) .* max(0,t-20);
m = max(m, mco);

% Angle of attack
alpha = 0.05 * sin(0.1*t);  % [rad]

% Fin geometry & velocity (for computeTotalFinXcp)
vlauncher = linspace(0, 1500, length(t));  % Example speed profile [m/s]
vsound    = 340;                            % Speed of sound [m/s]

cmac = 2;   % Fin mean chord [m]
be   = 1.5; % 2*fin axial base [m]
Se   = 3;   % 2*fin surface [m^2]

%% COMPUTE CG AND CP
for i = 1:length(t)

    % --- CG
    Xcg_curve(i) = computeXCG(N(i), lc1, lc2, lco, m(i), mc2, mco);

    % --- Cone CP
    Xcp_cone_curve(i) = computeXcp(N(i), alpha(i), lc1, lc2, lco);

    % --- Flare CP
    if N(i) == 2
        Xcp_flare_curve(i) = computeFlareXcp(lc1, lc2, lco, d, hf, db);
    else
        Xcp_flare_curve(i) = computeXcp(N(i), alpha(i), lc1, lc2, lco);
    end

    % --- Fin + flare CP
    Xcp_fin_flare_curve(i) = computeFinFlareXcp(N(i), lc2, lc1, lco, alpha(i), d, hf, db);

    % --- Total fin CP (continuous over all stages)
    Xcp_fin_curve(i) = computetotalFinXcp(N(i), alpha(i), lc1, lc2, lco, ...
                                           vlauncher(i), vsound, cmac, be, Se, Sref);
end

%% PLOT: CG + CP
figure;
plot(t, Xcg_curve, 'LineWidth', 2, 'DisplayName', 'X_{cg}');
hold on;
plot(t, Xcp_cone_curve, '--', 'LineWidth', 2, 'DisplayName', 'X_{cp} cone');
plot(t, Xcp_flare_curve, ':', 'LineWidth', 2, 'DisplayName', 'X_{cp} flare');
plot(t, Xcp_fin_flare_curve, '-.', 'LineWidth', 2, 'DisplayName', 'X_{cp} fin + flare');
plot(t, Xcp_fin_curve, 'k-', 'LineWidth', 2, 'DisplayName', 'X_{cp} total fins');
grid on;
xlabel('Time [s]');
ylabel('Position [m]');
title('Center of Gravity and Center of Pressure Evolution');
legend('Location', 'best');
ylim([0 max([Xcg_curve, Xcp_cone_curve, Xcp_flare_curve, Xcp_fin_flare_curve, Xcp_fin_curve])*1.1]);

%% DISTANCE FROM BASE
L_total = lco + lc1 + lc2;

Xcg_from_base             = L_total - Xcg_curve;
Xcp_cone_from_base        = L_total - Xcp_cone_curve;
Xcp_flare_from_base       = L_total - Xcp_flare_curve;
Xcp_fin_flare_from_base   = L_total - Xcp_fin_flare_curve;
Xcp_fin_from_base         = L_total - Xcp_fin_curve;

figure;
plot(t, Xcg_from_base, 'LineWidth', 2, 'DisplayName', 'L - X_{cg}');
hold on;
plot(t, Xcp_cone_from_base, '--', 'LineWidth', 2, 'DisplayName', 'L - X_{cp} cone');
plot(t, Xcp_flare_from_base, ':', 'LineWidth', 2, 'DisplayName', 'L - X_{cp} flare');
plot(t, Xcp_fin_flare_from_base, '-.', 'LineWidth', 2, 'DisplayName', 'L - X_{cp} fin + flare');
plot(t, Xcp_fin_from_base, 'k-', 'LineWidth', 2, 'DisplayName', 'L - X_{cp} total fins');
grid on;
xlabel('Time [s]');
ylabel('Distance from base [m]');
title('Distance from Base vs Time');
legend('Location', 'best');

%% STABILITY MARGIN
margin_cone      = (Xcp_cone_curve  - Xcg_curve) / d;
margin_flare     = (Xcp_flare_curve - Xcg_curve) / d;
margin_fin_flare = (Xcp_fin_flare_curve - Xcg_curve) / d;
margin_fin        = (Xcp_fin_curve - Xcg_curve) / d;

figure;
plot(t, margin_cone, 'LineWidth', 2, 'DisplayName', 'Margin (cone CP)');
hold on;
plot(t, margin_flare, '--', 'LineWidth', 2, 'DisplayName', 'Margin (flare CP)');
plot(t, margin_fin_flare, '-.', 'LineWidth', 2, 'DisplayName', 'Margin (fin + flare)');
plot(t, margin_fin, 'k-', 'LineWidth', 2, 'DisplayName', 'Margin (total fins)');
yline(0, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Neutral Stability');
grid on;
xlabel('Time [s]');
ylabel('(X_{cp} - X_{cg}) / d');
title('Static Stability Margin Over Time');
legend('Location', 'best');
