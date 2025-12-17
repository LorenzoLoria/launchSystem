clear all;
clc;
close all

% ==========================  DATI ========================================

addpath(genpath('..\..\'))

[mission, settings] = dataStructGlobal;

mission.structure.alphaQmax = 0;

% launcher = [nStages, nMotore1, nMotore2, nMotore3, %massa1, %massa2,
% %massa3];
launcher = [2,2,3,4,0.459952176990556, 0.753370531158904, 0.634795741885559];

for i = 1:launcher(1)
    configuration.stage{i}.engine = mission.engines{launcher(1+i)};
end

[mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);
thrustDataGA = load('thrustdataVecTraj.mat','xGATraj');

thrustDataVecFMC(:,:,1) = [0.902082365568723	1.480898931628005
                            0.999984156345040	23.253294859564580
                            0.900002678098914	52.979241033086943
                            0.900000000000007	59.571815331701984
                            0.903941814015555	55.058714159781090];


thrustDataVecFMC(:,:,2 ) = [0.400917809388214	65.122710138507202
                            0.964494359624014	79.658359202140389
                            0.975968800776448	91.801043507018605
                            0.992714640706230	89.085172454227390
                            0.993244065056187	99.345740209598944];



[timeCollocation, stateCollocation] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustDataVecFMC);
%% Validation 1 --> Jorgensen-Allen model (body)
% - Figure 1) Jorgensen-Allen: fineness ratio effect (cone-cylinder, M_inf = 2.86)
% - Figure 2) Jorgensen-Allen: M_inf effect (Body 9, ogive-cylinder, L/d = 11)
validate_JorgensenAllen(settings)


%% Validation 2 --> Jerger/Fleeman model for friction drag (body)
% Validazione del modello di Jerger/Fleeman per la friction drag del corpo:
%
%   (C_D0)_f = 0.053 * (L/d) * ( M / (q_psf * L_ft) )^0.2   (unità anglosassoni)
%
%   (C_D0)_f ≈ 0.091 * (L/d) * ( M / (q_SI * L_SI) )^0.2    (unità SI)
%
% Casi considerati:
%  1) "Rocket baseline" di Fleeman (Tactical Missile Design)
%  2) Mortaio guidato FOI PGMM 120 mm (FOI-R--2618--SE)

% validate_JergerFleeman


%% Validation 3 --> Fleeman model for x_CP (body)
% Funzione per la validazione della formula di Fleeman per la posizione del centro di pressione
% con valori di L_body/L_nose = 1, 2, 5, 10 e angolo di attacco alpha da 0° a 90°.
validate_Fleeman_xCP(settings)