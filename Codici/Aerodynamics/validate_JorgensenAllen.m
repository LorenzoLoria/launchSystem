function [CN_ref_cone_Fleeman] = validate_JorgensenAllen
%
% VALIDATE_JORGENSENALLEN
% Validazione del modello Jørgensen–Allen per C_N usando i dati
% di Jorgensen, NASA TN D-6996 (1973):
%    - bodies 3–5 (cone–cylinder)
%    - body 9 (ogive-cylinder)
%
% Riferimento:
% L. H. Jorgensen, "Prediction of Static Aerodynamic Characteristics
% for Space-Shuttle-Like and Other Bodies at Angles of Attack from
% 0 to 180 deg", NASA TN D-6996, 1973.
%
% Questo script:
%  - calcola C_N(alpha) per i cone–cylinder bodies 3–5 a M=2.86, Re_d=1.25e6
%  - studia l'effetto di M per il body 9 (L/d=11)

clc;

%% Parametri di riferimento (da TN D-6996)
M_ref = 2.86;        % Mach di prova per i 9 corpi (Moo = 2.86)
Re_d  = 1.25e6;      % Reynolds basato sul diametro (Re = 1.25e6)
alpha_deg = linspace(0,180,181);  % range di attacco [deg] per la validazione

% Area di riferimento A = A_cyl = pi d^2 / 4, quindi A/d^2:
A_over_d2 = pi/4;

%% Geometria dei cone–cylinder bodies (fig. 9, riga "Cone")
% Tabella di fig. 9 (estratto per bodies 3–5):
% Body  L/d   ...   Ap/d^2   V/d^3    Xc/d    As/d^2   Nose
%   3     7   ...   5.500    3.925    4.183   17.34    Cone
%   4     9   ...   7.500    5.495    5.200   23.62    Cone
%   5    11   ...   9.500    7.065    6.211   29.91    Cone
%   9    11   ...   9.340    6.811    6.255   29.38    Ogive

bodies = struct([]);

% Body 3
bodies(1).id         = 3;
bodies(1).L_over_d   = 7.0;
bodies(1).Ap_over_d2 = 5.500;   % Ap/d^2
bodies(1).nose       = 'cone';

% Body 4
bodies(2).id         = 4;
bodies(2).L_over_d   = 9.0;
bodies(2).Ap_over_d2 = 7.500;   % Ap/d^2
bodies(2).nose       = 'cone';

% Body 5
bodies(3).id         = 5;
bodies(3).L_over_d   = 11.0;
bodies(3).Ap_over_d2 = 9.500;   % Ap/d^2
bodies(3).nose       = 'cone';

% Body 9
bodies(4).id         = 9;
bodies(4).L_over_d   = 11.0;
bodies(4).Ap_over_d2 = 9.340;   % Ap/d^2
bodies(4).nose       = 'ogive';

% Completa i coefficienti di area normalizzati con A = pi d^2/4
for k = 1:numel(bodies)
    bodies(k).Ab_over_A = 1.0;  % cilindro: area di base = area di riferimento
    bodies(k).Ap_over_A = bodies(k).Ap_over_d2 / A_over_d2; % (Ap/d^2)/(A/d^2)
end

%% VALIDAZIONE 1: M = 2.86, variazione fineness ratio (L/d = 7, 9, 11)

figure; hold on; grid on; box on;
colors = lines(numel(bodies));

for k = 1:(numel(bodies)-1)
    CN = CN_JA(alpha_deg, M_ref, Re_d, bodies(k));
    plot(alpha_deg, CN, 'LineWidth', 3, ...
         'Color', colors(k,:), ...
         'DisplayName', sprintf('Body %d, L/d = %.0f', ...
                      bodies(k).id, bodies(k).L_over_d));

    if k == 1
        CN_ref_cone_Fleeman.b7 = CN;
    elseif k == 2
        CN_ref_cone_Fleeman.b9 = CN;
    else
        CN_ref_cone_Fleeman.b11 = CN;
    end

end

xlabel('\alpha [deg]');
ylabel('C_N');
ylim([0 20]);          % Limiti dell’asse y
yticks(0:5:20);        % Tick ogni 5
xlim([0 180])
%title('Jorgensen-Allen: fineness ratio effect (cone-cylinder, M_\infty = 2.86)');
lgd = legend('Location','best');  
lgd.FontSize = 15;               




%% VALIDAZIONE 2: L/d = 11 (Body 9), effetto del Mach

bodyLd11 = bodies(4);   % Body 9, L/d = 11, ogive–cylinder

% Valori di Mach per lo sweep (puoi modificarli a piacere)
Mach_list = [0.3 1.5 2.9 7.0];

figure; hold on; grid on; box on;
colors = lines(numel(Mach_list));

for i = 1:numel(Mach_list)
    M = Mach_list(i);
    CN_M = CN_JA(alpha_deg, M, Re_d, bodyLd11);
    plot(alpha_deg, CN_M, 'LineWidth', 3, ...
         'Color', colors(i,:), ...
         'DisplayName', sprintf('M_\\infty = %.2g', M));
end

xlabel('\alpha [deg]');
ylabel('C_N');
ylim([0 20]);          % Limiti dell’asse y  
yticks(0:4:20);        % Tick ogni 4
xlim([0 180])
%title(sprintf('Jorgensen-Allen: M_\\infty effect (Body %d, ogive-cylinder, L/d = %.0f)', ...
%      bodyLd11.id, bodyLd11.L_over_d));
lgd = legend('Location','best');  
lgd.FontSize = 15;   



end 





% =====================================================================
%                 SOTTOFUNZIONI
% =====================================================================

function CN = CN_JA(alpha_deg, M_inf, Re_d, body)
% CN_JA  Calcola C_N(alpha) con il modello Jørgensen–Allen
%
%  CN(alpha) = (Ab/A) sin(2a) cos(a/2) + eta_eff * Cdn * (Ap/A) sin^2(a')
%  con a' = a   per 0 <= a <= 90 deg
%      a' = 180-a per 90 < a <= 180 deg
%
% alpha_deg : vettore di angoli d'attacco [deg], tipicamente 0–180
% M_inf     : Mach di freestream
% Re_d      : Reynolds basato sul diametro
% body      : struct con campi:
%             L_over_d, Ab_over_A, Ap_over_A

alpha_deg = alpha_deg(:).';      % riga
a_rad     = deg2rad(alpha_deg);  % [rad]

% Angolo ridotto a' in [0, 90] gradi per il termine di crossflow
alpha_prime_deg = alpha_deg;
mask = alpha_prime_deg > 90;
alpha_prime_deg(mask) = 180 - alpha_prime_deg(mask);
a_p_rad = deg2rad(alpha_prime_deg);

% Componenti di Mach e Reynolds in crossflow
M_n  = M_inf .* sin(a_p_rad);
Re_n = Re_d  .* sin(a_p_rad);

% Crossflow drag coefficient Cdn(Mn, Ren)
Cdn = Cdn_circ_cyl(M_n, Re_n);

% Fattore di proporzionalità base eta(L/d, M_inf) (come prima)
eta0 = eta_crossflow(body.L_over_d, M_inf);

% NUOVO: fattore di scala in funzione del Mach
gM   = Mach_factor(M_inf);

% Fattore efficace per il termine di crossflow
eta_eff = eta0 * gM;

% Termine slender-body (potenziale)
slender = body.Ab_over_A .* sin(2*a_rad) .* cos(a_rad/2);

% Termine crossflow (viscoso) scalato con il Mach
cross = eta_eff .* Cdn .* body.Ap_over_A .* (sin(a_p_rad)).^2;

% Somma: modello Jørgensen–Allen
CN = slender + cross;

end





function Cdn = Cdn_circ_cyl(Mn, ~)
Cdn_high = 1.3;
Mn = abs(Mn);

Mc = 0.05;  % scala "piccola" di Mach normale

% Cdn sale come ~ (Mn/Mc)^2 per Mn << Mc, e tende a Cdn_high per Mn >> Mc
Cdn = Cdn_high * (1 - exp(-(Mn./Mc).^2));
end




function eta = eta_crossflow(L_over_d, M_inf)
% ETA_CROSSFLOW  Fattore di proporzionalità eta per il termine di crossflow.
%
% Da Jorgensen (fig. 4) e letteratura associata:
%  - a Mach supersonici/ipersonici si può assumere eta ~ 1
%  - a Mach subsonici, eta cresce con L/d (drag di cilindro finito vs infinito)
%
% Qui uso una forma semplice:
%  eta = min( 0.05 (L/d) + 0.52 , 1.0 ) per M_inf < 1.2
%  eta = 1.0 per M_inf >= 1.2
%
% Se vuoi ancora più regolarità rispetto a M, puoi rendere anche questa
% transizione liscia (es. con una sigmoide in funzione di M_inf).

if M_inf >= 1.2
    eta = 1.0;
else
    eta = 0.05 * L_over_d + 0.52;
    eta = min(eta, 1.0);
end

end




function gM = Mach_factor(M_inf)
% MACH_FACTOR  Fattore di scala in funzione del Mach per il termine di crossflow.

if M_inf <= 0.3
    % per M <= 0.3 fisso il valore desiderato
    gM = 0.64;
elseif M_inf <= 1.5
    % interpolo linearmente da 0.64 (M=0.3) a 1.15 (M=1.5)
    M1 = 0.3;  g1 = 0.64;
    M2 = 1.5;  g2 = 1.15;
    t  = (M_inf - M1)/(M2 - M1);
    gM = g1 + (g2 - g1)*t;
elseif M_inf <= 2.86
    % interpolo da 1.15 (M=1.5) a 1.00 (M=2.86)
    M1 = 1.5;   g1 = 1.15;
    M2 = 2.9;   g2 = 1.00;
    t  = (M_inf - M1)/(M2 - M1);
    gM = g1 + (g2 - g1)*t;
else
    % per M >= 2.86 scendo da 1.00 (M=2.86) a 0.96 (M=7.0)
    M1 = 2.9;   g1 = 1.00;
    M2 = 7.0;   g2 = 0.96;
    t  = (M_inf - M1)/(M2 - M1);
    t  = max(min(t,1.0),0.0);
    gM = g1 + (g2 - g1)*t;
end

end
