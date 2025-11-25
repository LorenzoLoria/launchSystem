% test_CA_body.m
% Script di esempio per calcolare e plottare C_A,body al variare del Mach

clear; clc; close all;

%% 1) Definizione geometria del corpo
geom.L     = 3.0;          % lunghezza totale [m]
geom.d     = 0.3;          % diametro di riferimento [m]
geom.Lnose = 0.6;          % lunghezza nose [m]

geom.Aref  = pi*(geom.d/2)^2;   % area di riferimento [m^2]
geom.Anose = geom.Aref;         % se nose ha stesso diametro
geom.Abase = geom.Aref;         % base circolare piena

geom.Aexit = [];               % nessun ugello (metti valori se vuoi motore acceso)
geom.phi   = deg2rad(0);       % giunzione smooth -> phi ~ 0 rad

%% 2) Parametri di flusso (esempio)
flow.q = 20000;   % pressione dinamica [Pa] (es. ~ 20 kPa)

%% 3) Parametri modello subsonico per wave drag
% (da scegliere per raccordare bene con il modello transonico/supersonico)
a_sub = 0.0;
b_sub = 1.0;

%% 4) Vettore di Mach e angolo di attacco
M_vec  = linspace(0.1, 5.0, 1000);   % vettore di Mach
alpha  = 0.0;                       % [rad] -> asse del corpo allineato al flusso

%% 5) Fase propulsa o coasting
isPowered = false;   % false = coasting, true = motore acceso

%% 6) Calcolo del C_A,body(M)
CA_vec = CA_body(M_vec, alpha, geom, flow, isPowered, a_sub, b_sub);

%% 8) Plot C_A vs Mach
figure;
plot(M_vec, CA_vec, 'LineWidth', 1.5);
grid on;
xlabel('Mach number M');
ylabel('C_{A}^{body}');
title('Axial force coefficient of the body vs Mach');


