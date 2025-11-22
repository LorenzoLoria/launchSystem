function validate_JergerFleeman
% VALIDATE_JERGERFLEEMAN
%
% Validazione del modello di Jerger/Fleeman per la friction drag del corpo:
%
%   (C_D0)_f = 0.053 * (L/d) * ( M / (q_psf * L_ft) )^0.2   (unità anglosassoni)
%
%   (C_D0)_f ≈ 0.091 * (L/d) * ( M / (q_SI * L_SI) )^0.2    (unità SI)
%
% Casi considerati:
%  1) "Rocket baseline" di Fleeman (Tactical Missile Design)
%  2) Mortaio guidato FOI PGMM 120 mm (FOI-R--2618--SE)
%
% Riferimenti:
% - Fleeman, "Tactical Missile Design", AIAA, 2006.
%
% - FOI-R--2618--SE, "Flight Mechanical Modeling of a Precision
%   Guided Mortar Munition", FOI, 2008.

clc;

%% Costanti di conversione
PA_PER_PSF  = 47.880258;   % 1 psf = 47.880258 Pa
FT_PER_M    = 3.280839895; % 1 m   = 3.28084 ft

fprintf('=============================================\n');
fprintf(' Validazione modello Jerger/Fleeman\n');
fprintf('=============================================\n\n');

%% ==========================================================
%  CASO 1: "Rocket baseline" di Fleeman (Mach 2, high L/d)
% ===========================================================

% Dati dal caso di esempio di Fleeman:
L_over_d1 = 18.0;          % fineness ratio L/d
L_ft1     = 12.0;          % lunghezza in ft
M1        = 2.0;           % numero di Mach
q_psf1    = 2725.0;        % pressione dinamica in psf

% Valore di riferimento riportato da Fleeman (body friction drag)
CD0_ref1  = 0.14;

% Calcolo in unità anglosassoni
CD0_1_imp = CD0_JergerFleeman_imperial(M1, q_psf1, L_over_d1, L_ft1);

% Calcolo equivalente in unità SI (stesse condizioni convertite)
L_SI1   = L_ft1 / FT_PER_M;             % [m]
q_SI1   = q_psf1 * PA_PER_PSF;          % [Pa]
CD0_1_SI = CD0_JergerFleeman_SI(M1, q_SI1, L_over_d1, L_SI1);

% Errori rispetto al riferimento Fleeman
err1_abs_imp = CD0_1_imp - CD0_ref1;
err1_rel_imp = 100 * err1_abs_imp / CD0_ref1;

fprintf('Caso 1: Fleeman "rocket baseline"\n');
fprintf('  Input: L/d = %.1f, L = %.1f ft, M = %.2f, q = %.0f psf\n', ...
        L_over_d1, L_ft1, M1, q_psf1);
fprintf('  (C_D0)_f (imperial) = %.4f\n', CD0_1_imp);
fprintf('  (C_D0)_f (SI)       = %.4f\n', CD0_1_SI);
fprintf('  Riferimento         = %.4f\n', CD0_ref1);
fprintf('  Errore assoluto     = %+ .4f\n', err1_abs_imp);
fprintf('  Errore percentuale  = %+ .2f %%\n\n', err1_rel_imp);


%% ==========================================================
%  CASO 2: FOI PGMM 120 mm (subsonico, Mach ~ 0.37)
% ===========================================================

% Dati geometrici del mortaio guidato (FOI-R--2618--SE)
d2      = 0.12;    % diametro [m] (120 mm)
L_SI2   = 0.718;   % lunghezza aerodinamica corpo [m]
L_over_d2 = L_SI2 / d2;

% Condizioni di volo di esempio (FOI): V ≈ 125 m/s, h ≈ 1500 m
V2      = 125.0;   % velocità [m/s]
h2      = 1500.0;  % quota [m]

% Atmosfera ISA semplificata per ricavare M2 e q2
[rho2, a2] = ISA_atmosphere(h2);
M2   = V2 / a2;
q_SI2 = 0.5 * rho2 * V2^2;            % [Pa]
q_psf2 = q_SI2 / PA_PER_PSF;         % [psf]
L_ft2  = L_SI2 * FT_PER_M;           % [ft]

% Valore di riferimento (ordine di grandezza) da FOI (curva C_D0,bodyfric)
CD0_ref2 = 0.08;   % ≈ valore letto da FOI-R--2618--SE a V=125 m/s

% Calcolo con formula Jerger/Fleeman in unità SI
CD0_2_SI  = CD0_JergerFleeman_SI(M2, q_SI2, L_over_d2, L_SI2);

% Calcolo equivalente in unità anglosassoni (consistenza)
CD0_2_imp = CD0_JergerFleeman_imperial(M2, q_psf2, L_over_d2, L_ft2);

% Errori rispetto al riferimento FOI
err2_abs = CD0_2_SI - CD0_ref2;
err2_rel = 100 * err2_abs / CD0_ref2;

fprintf('Caso 2: FOI PGMM 120 mm (mortar shell)\n');
fprintf('  Input: L/d ≈ %.2f, L = %.3f m, d = %.3f m\n', L_over_d2, L_SI2, d2);
fprintf('         V = %.1f m/s, h = %.0f m\n', V2, h2);
fprintf('         M ≈ %.3f, q ≈ %.1f Pa (%.1f psf)\n', M2, q_SI2, q_psf2);
fprintf('  (C_D0)_f (SI)       = %.4f\n', CD0_2_SI);
fprintf('  (C_D0)_f (imperial) = %.4f\n', CD0_2_imp);
fprintf('  Riferimento (FOI)   ≈ %.4f\n', CD0_ref2);
fprintf('  Errore assoluto     = %+ .4f\n', err2_abs);
fprintf('  Errore percentuale  = %+ .2f %%\n\n', err2_rel);

fprintf('=============================================\n');
fprintf(' Fine validazione Jerger/Fleeman\n');
fprintf('=============================================\n');

end


% ==========================================================
%        SOTTOFUNZIONI: Jerger/Fleeman + atmosfera ISA
% ==========================================================

function CD0 = CD0_JergerFleeman_imperial(M, q_psf, L_over_d, L_ft)
% CD0_JERGERFLEEMAN_IMPERIAL
%   Formula originale Jerger/Fleeman in unità anglosassoni:
%
%   (C_D0)_f = 0.053 * (L/d) * ( M / (q_psf * L_ft) )^0.2
%
%   dove:
%     - M      : Mach
%     - q_psf  : pressione dinamica [psf]
%     - L_over_d = L/d
%     - L_ft   : lunghezza del corpo [ft]

CD0 = 0.053 .* L_over_d .* ( M ./ (q_psf .* L_ft) ).^0.2;

end


function CD0 = CD0_JergerFleeman_SI(M, q_SI, L_over_d, L_SI)
% CD0_JERGERFLEEMAN_SI
%   Formula Jerger/Fleeman riadattata per uso pratico in unità SI:
%
%   (C_D0)_f ≈ 0.091 * (L/d) * ( M / (q_SI * L_SI) )^0.2
%
%   dove:
%     - M      : Mach
%     - q_SI   : pressione dinamica [Pa]
%     - L_over_d = L/d
%     - L_SI   : lunghezza del corpo [m]

CD0 = 0.091 .* L_over_d .* ( M ./ (q_SI .* L_SI) ).^0.2;

end


function [rho, a] = ISA_atmosphere(h)
% ISA_ATMOSPHERE  Atmosfera standard semplificata (strato troposfera, h < 11 km)
% Restituisce densità rho [kg/m^3] e velocità del suono a [m/s]
% in funzione della quota geometrica h [m].

% Costanti ISA
T0 = 288.15;       % K
p0 = 101325;       % Pa
L  = 0.0065;       % K/m (gradiente troposfera)
g  = 9.80665;      % m/s^2
R  = 287.05;       % J/(kg K)
gamma = 1.4;

if h < 11000
    T = T0 - L*h;
    p = p0 * (T/T0)^(g/(R*L));
else
    % Per questo esercizio non ci serve oltre, ma mettiamo un fallback semplice
    T = 216.65;
    p = 22632 * exp(-g*(h-11000)/(R*T));
end

rho = p/(R*T);
a   = sqrt(gamma * R * T);

end
