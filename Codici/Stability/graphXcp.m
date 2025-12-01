close all
clear all

%% PARAMETERS
lc1 = 20;
lc2 = 10;
lco = 3;

m0  = 210.4e3;
mc1 = 162.5e3;
mc2 = 40.6e3;
mco = m0 - mc1 - mc2;

d  = 3;
hf = 0.806 * d;
db = 1.094 * d;

flare_ratio = (db - d)/(2*hf);

t = 0:0.5:50;

Xcg_curve       = zeros(size(t));
Xcp_cone_curve  = zeros(size(t));
Xcp_flare_curve = zeros(size(t));
Xcp_fin_curve   = zeros(size(t));   % <-- NEW

%% DYNAMIC PARAMETERS
N = ones(size(t));
N(t <= 20) = 2;
N(t > 20 & t <= 40) = 1;
N(t > 40) = 0;

m = m0 - (m0-mc2-mco)/20 .* min(t,20) ...
      - (mc2/20) .* max(0,t-20);
m = max(m, mco);

alpha = 0.05 * sin(0.1*t); 

%% COMPUTE CG AND CP
for i = 1:length(t)

    % CG
    Xcg_curve(i) = computeXCG(N(i), lc1, lc2, lco, m(i), mc2, mco);

    % Cone CP (always available)
    Xcp_cone_curve(i) = computeXcp(N(i), alpha(i), lc1, lc2, lco);

    % Flare CP: only when stage 1 attached
    if N(i) == 2
        Xcp_flare_curve(i) = computeFlareXcp(lc1, lc2, lco, d, hf, db);
    else
        Xcp_flare_curve(i) = computeXcp(N(i), alpha(i), lc1, lc2, lco);
    end

    % NEW — Fin CP model (always computed)
    Xcp_fin_curve(i) = computeFinXcp(N(i), lc2, lc1, lco, alpha(i), d, hf, db);

end

%% PLOT CG + CP
figure;
plot(t, Xcg_curve, 'LineWidth', 2, 'DisplayName', 'X_{cg}');
hold on;
plot(t, Xcp_cone_curve, '--', 'LineWidth', 2, 'DisplayName', 'X_{cp} cone');
plot(t, Xcp_flare_curve, ':', 'LineWidth', 2, 'DisplayName', 'X_{cp} flare');
plot(t, Xcp_fin_curve, '-.', 'LineWidth', 2, 'DisplayName', 'X_{cp} fin');   % <-- NEW

grid on;
xlabel('Time [s]');
ylabel('Position [m]');
title('Center of Gravity and Center of Pressure Evolution');
legend('Location', 'best');
ylim([0 max([Xcg_curve Xcp_cone_curve Xcp_flare_curve Xcp_fin_curve])*1.1]);

%% FROM BASE
L_total = lco + lc1 + lc2;

Xcg_from_base       = L_total - Xcg_curve;
Xcp_cone_from_base  = L_total - Xcp_cone_curve;
Xcp_flare_from_base = L_total - Xcp_flare_curve;
Xcp_fin_from_base   = L_total - Xcp_fin_curve;    % <-- NEW

figure;
plot(t, Xcg_from_base, 'LineWidth', 2, 'DisplayName', 'L - X_{cg}');
hold on;
plot(t, Xcp_cone_from_base, '--', 'LineWidth', 2, 'DisplayName', 'L - X_{cp} cone');
plot(t, Xcp_flare_from_base, ':', 'LineWidth', 2, 'DisplayName', 'L - X_{cp} flare');
plot(t, Xcp_fin_from_base, '-.', 'LineWidth', 2, 'DisplayName', 'L - X_{cp} fin');  % <-- NEW

grid on;
xlabel('Time [s]');
ylabel('Distance from base [m]');
title('Distance from Base vs Time');
legend('Location', 'best');

%% STABILITY MARGIN
margin_cone  = (Xcp_cone_curve  - Xcg_curve) / d;
margin_flare = (Xcp_flare_curve - Xcg_curve) / d;
margin_fin   = (Xcp_fin_curve   - Xcg_curve) / d;  % <-- NEW

figure;
plot(t, margin_cone, 'LineWidth', 2, 'DisplayName', 'Margin (cone CP)');
hold on;
plot(t, margin_flare, '--', 'LineWidth', 2, 'DisplayName', 'Margin (flare CP)');
plot(t, margin_fin, '-.', 'LineWidth', 2, 'DisplayName', 'Margin (fin CP)'); % <-- NEW

yline(0, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Neutral Stability');

grid on;
xlabel('Time [s]');
ylabel('(X_{cp} - X_{cg}) / d');
title('Static Stability Margin Over Time');
legend('Location', 'best');
