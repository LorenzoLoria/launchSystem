% test_CN_body_alpha.m
% Test del coefficiente di forza normale del body C_N al variare di alpha

clear; clc; close all;

%% 1) Definizione geometria del corpo (esempio)
geom.d     = 0.3;                     % diametro [m]
geom.l     = 3.0;                     % lunghezza ell [m]
geom.Aref  = pi*(geom.d/2)^2;         % area di riferimento A [m^2]
geom.Ab    = geom.Aref;               % area caratteristica A_b (sezione max)
geom.Ap    = geom.l * geom.d;         % area proiettata in crossflow A_p ~ L*d

%% 2) Parametri aerodinamici (esempio)
aero.Cdn = 1.2;   % crossflow drag coefficient del cilindro equivalente

%% 3) Mach fissato (puoi cambiarlo)
M1 = 0.8;          % subsonico
M2 = 2.0;          % supersonico

%% 4) Vettore di angoli di attacco alpha
alpha_deg = linspace(0, 15, 200);         % [gradi]
alpha_rad = deg2rad(alpha_deg);           % [rad]

%% 5) Calcolo di C_N(alpha) usando CN_body (loop su alpha)
CN1 = zeros(size(alpha_rad));
CN2 = zeros(size(alpha_rad));

for i = 1:numel(alpha_rad)
    CN1(i) = CN_body(M1, alpha_rad(i), geom, aero);
    CN2(i) = CN_body(M2, alpha_rad(i), geom, aero);
end

%% 6) Plot C_N vs alpha
figure;
plot(alpha_deg, CN1, 'LineWidth', 1.5);
grid on;
hold on;
xlabel('\alpha [deg]');
plot(alpha_deg, CN2, 'LineWidth', 1.5);
ylabel('C_{N}^{body}');
title('Normal force coefficient of the body vs alpha');
legend('M = 0.8', 'M = 2.0');
