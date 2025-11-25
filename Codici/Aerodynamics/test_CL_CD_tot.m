% test_CL_CD_tot.m
% Calcolo CL_tot(alpha) e CD_tot(M) per razzo con body + fins

clear; clc; close all;

%% 1) Parametri geometrici BODY (dati arbitrari)
geom.d     = 0.3;                     % diametro [m]
geom.l     = 3.0;                     % lunghezza ell [m]
geom.L     = geom.l;                  % per CA_body
geom.Lnose = 0.6;                     % lunghezza nose [m]

geom.Aref  = pi*(geom.d/2)^2;         % area di riferimento A [m^2]
geom.Ab    = geom.Aref;               % area caratteristica A_b (sezione max)
geom.Ap    = geom.l * geom.d;         % area proiettata in crossflow A_p ~ L*d

geom.Anose = geom.Aref;               % nose stesso diametro
geom.Abase = geom.Aref;               % base circolare piena
geom.Aexit = [];                      % nessun ugello
geom.phi   = deg2rad(0);              % giunzione smooth

Sref       = geom.Aref;

%% 2) Parametri geometrici FINS (singola aletta triangolare, dati arbitrari)
cr   = 0.35;                      % root chord [m]
ct   = 0.0;                       % tip chord [m] (triangolo puro)
s    = 0.20;                      % semispan [m]

Se   = 0.5 * cr * s;              % area della fin [m^2]
be   = s;                         % span equivalente ~ semispan reale

% MAC triangolo
cmac = (2/3) * cr;                % mean aerodynamic chord [m]
tmac = 0.08 * cmac;               % spessore max ≈ 8% MAC [m]

b    = 2 * cr;                    % 2 * base della fin [m]

delta_le  = 45;                   % leading edge sweep [deg]
lambda_le = 0;                    % fin base angle [deg]

Nfins = 1;                        % NUMERO DI ALETTE

%% 3) Parametri aerodinamici / modello
aero.Cdn = 1.2;                   % crossflow drag coeff cilindro equivalente

% Parametri subsonici per wave drag body (CA_body)
a_sub = 0.0;
b_sub = 1.0;

%% 4) Atmosfera e suono
vsound = 340;                     % m/s
rho    = 1.225;                   % kg/m^3

%% 5) CL_tot(alpha) per due Mach (0.8 e 5)
M_vec = [0.8 5];                  % Mach di riferimento per CL(alpha)
alpha_deg = [0 : 0.01 : 15];
alpha_rad = deg2rad(alpha_deg);

CL_tot_alpha = zeros(numel(alpha_rad), numel(M_vec));

for j = 1:numel(M_vec)

    M = M_vec(j);
    vlauncher = M * vsound;
    q = 0.5 * rho * vlauncher^2;

    flow.q = q; % struct per CA_body

    for i = 1:numel(alpha_rad)
        alpha = alpha_rad(i);

        % --- Body: CN e CA
        CN_body_val = CN_body(M, alpha, geom, aero);
        CA_body_val = CA_body(M, alpha, geom, flow, false, a_sub, b_sub);

        % --- Fins: CN e CA per UNA ALETTTA
        % Mach normale alla fin (approx con sweep bordo d'attacco)
        M_ale = M * cosd(delta_le);

        [CN_fin_single, CD0_fric_single, CD0_wave_single] = AerodynCoefFins_new(alpha, vlauncher, vsound, be, Se, q, Sref, cmac, delta_le, lambda_le, b, tmac);

        CN_fins_tot = Nfins * CN_fin_single;
        CD0_fins_tot = Nfins * (CD0_fric_single + CD0_wave_single);
        CA_fins_tot = CD0_fins_tot * cos(alpha)^2;

        % Somma totali
        CN_tot = CN_body_val + CN_fins_tot;
        CA_tot = CA_body_val + CA_fins_tot;

        % CL / CD totali
        [CL_val, ~] = calculate_CL_CD(CN_tot, CA_tot, rad2deg(alpha));

        CL_tot_alpha(i, j) = CL_val;
    end
end

%% 6) CD_tot(M) per alpha fisso (5 deg)
alpha_CD_deg = 5;                  % angolo di attacco per CD(M)
alpha_CD_rad = deg2rad(alpha_CD_deg);

M_vec = linspace(0.01, 5.0, 1000); % evito M=0 per q=0
CD_tot_M = zeros(size(M_vec));

for k = 1:numel(M_vec)

    M = M_vec(k);
    vlauncher = M * vsound;
    q = 0.5 * rho * vlauncher^2;
    flow.q = q;

    % --- Body
    CN_body_val = CN_body(M, alpha_CD_rad, geom, aero);
    CA_body_val = CA_body(M, alpha_CD_rad, geom, flow, false, a_sub, b_sub);

    % --- Fins
    M_ale = M * cosd(delta_le);

    [CN_fin_single, CD0_fric_single, CD0_wave_single] = AerodynCoefFins_new(alpha_CD_rad, vlauncher, vsound, be, Se, q, Sref, cmac, delta_le, lambda_le, b, tmac);

    CN_fins_tot = Nfins * CN_fin_single;
    CD0_fins_tot = Nfins * (CD0_fric_single + CD0_wave_single);
    CA_fins_tot = CD0_fins_tot * cos(alpha_CD_rad)^2;

    CN_tot = CN_body_val + CN_fins_tot;
    CA_tot = CA_body_val + CA_fins_tot;

    [~, CD_val] = calculate_CL_CD(CN_tot, CA_tot, alpha_CD_deg);

    CD_tot_M(k) = CD_val;
end

%% 7) PLOT: CL_tot(alpha) per i due Mach
figure;
plot(alpha_deg, CL_tot_alpha(:,1), 'LineWidth', 1.5); hold on;
plot(alpha_deg, CL_tot_alpha(:,2), 'LineWidth', 1.5);
grid on;
xlabel('\alpha [deg]');
ylabel('C_L^{tot}');
title('Total lift coefficient vs \alpha (body + fins)');
legend(sprintf('M = %.1f', M_vec(1)), ...
       sprintf('M = %.1f', M_vec(2)), 'Location', 'northwest');

%% 8) PLOT: CD_tot(M) per alpha fisso
figure;
plot(M_vec, CD_tot_M, 'LineWidth', 1.5);
grid on;
xlabel('Mach number M [-]');
ylabel('C_D^{tot}');
title(sprintf('Total drag coefficient vs Mach (body + fins, \\alpha = %.1f^\\circ)', alpha_CD_deg));

