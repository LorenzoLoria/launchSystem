clear;
clc;

%% Validation 1 --> Jorgensen-Allen model (body)
% - Figure 1) Jorgensen-Allen: fineness ratio effect (cone-cylinder, M_inf = 2.86)
% - Figure 2) Jorgensen-Allen: M_inf effect (Body 9, ogive-cylinder, L/d = 11)
validate_JorgensenAllen


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
validate_JergerFleeman


%% Validation 3 --> Fleeman model for x_CP (body)
% Funzione per la validazione della formula di Fleeman per la posizione del centro di pressione
% con valori di L_body/L_nose = 1, 2, 5, 10 e angolo di attacco alpha da 0° a 90°.
validate_Fleeman_xCP