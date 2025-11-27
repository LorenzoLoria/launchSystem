clear all;
clc;

% ================= DATI MISSIONE (Costanti) =================
mission.environment.g0 = 9.81; 
mission.launcher.nx = 5;               % Max accelerazione longitudinale (g)
mission.launcher.structures.SF = 1.25; % Safety Factor standard
mission.launcher.structures.ultimate = 480e6; % Alluminio 7075-T6 (Pa)
mission.launcher.structures.E = 71e9;         % Modulo Young (Pa)
mission.launcher.engines.oxDens = 1141;       % Densità LOX (kg/m^3)
mission.launcher.structures.rho = 2800;       % Densità Alluminio (kg/m^3)

% ================= DATI VETTORIALI (5 Componenti) =================
nComponents = 5;

% 1. GEOMETRIA
% [Adapter, S2_Tank, Interstage, S1_Tank_Ox, S1_Tank_Fuel]
% Nota: Definisco il diametro come vettore perché varia tra gli stadi
mission.launcher.diameter = [3.0, 3.0, 3.7, 3.7, 3.7]; % [m]
mission.launcher.length =   [1.5, 5.0, 2.0, 20, 10]; % [m]

% 2. PRESSIONI INTERNE [Pa]
% 0 = Strutture a secco, ~3-4 bar = Serbatoi
mission.launcher.tankPressure = [0, 3.5e5, 0, 4.0e5, 4.0e5]; 

% 3. CARICHI (Azioni Interne Cumulative)
% N: Compressione (aumenta scendendo verso il basso per il peso accumulato)
mission.launcher.structures.N = [50e3, 250e3, 400e3, 2.5e6, 4.0e6]; % [N]

% T: Taglio (picchi nelle zone interstadio o max-q)
mission.launcher.structures.T = [5e3, 20e3, 80e3, 150e3, 100e3]; % [N]

% Mb: Momento Flettente (Max al centro del razzo durante le raffiche)
mission.launcher.structures.Mb = [2e4, 5e5, 2.5e6, 4.0e6, 1.5e6]; % [Nm]


% ================= ESECUZIONE =================
fprintf('--- Inizio Sizing Strutturale (%d Componenti) ---\n', nComponents);

tic
% Chiamata alla funzione (Assicurati che thicknessFunction.m sia nella stessa cartella)
[massStruct, t, stressMatrix] = thicknessFunction(mission, nComponents);
timeElapsed = toc;

% ================= OUTPUT RISULTATI =================
fprintf('Calcolo completato in %.4f secondi.\n\n', timeElapsed);

% Creazione Tabella per visualizzare i risultati in modo leggibile
Componenti = {'Payload Adapter'; 'S2 Tank (Pres)'; 'Interstage (Unp)'; 'S1 LOX (Pres)'; 'S1 Fuel (Pres)'};
Spessore_mm = t * 1000;
Massa_kg = massStruct;
Sigma_Buckling_MPa = stressMatrix(:,6) / 1e6; % Colonna 6 è il Buckling Stress
Sigma_Axial_MPa = stressMatrix(:,1) / 1e6;    % Colonna 1 è lo stress Longitudinale

ResultsTable = table(Componenti, Spessore_mm, Massa_kg, Sigma_Axial_MPa, Sigma_Buckling_MPa);
disp(ResultsTable);

% fprintf('Spessore Massimo del Lanciatore: %.2f mm\n', tMax * 1000);
fprintf('Massa Totale Strutturale: %.2f kg\n', sum(massStruct));