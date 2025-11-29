function CA = CA_body(M, alpha, geom, flow, isPowered, a_sub, b_sub)
% CA_body  Calcolo del coefficiente assiale del body C_A,body
%          usando i modelli di Jørgensen–Allen + Jerger/Fleeman + base drag.
%
%   CA = CA_body(M, alpha, geom, flow, isPowered, a_sub, b_sub)
%
% INPUT:
%   M         : numero di Mach (scalare o vettore)
%   alpha     : angolo d'attacco [rad] (stessa size di M o scalare)
%
%   geom      : struct con i parametri geometrici:
%               geom.l      = lunghezza totale corpo (ell_SI) [m]
%               geom.d      = diametro di riferimento [m]
%               geom.Lnose  = lunghezza del nose l_N [m]
%               geom.Aref   = area di riferimento A_ref [m^2]
%               geom.Anose  = area del nose A_nose [m^2]
%               geom.Abase  = area di base A_base [m^2]
%               geom.Aexit  = vettore aree ugelli A_exit,nozzle [m^2]
%               geom.phi    = angolo di giunzione nose-body φ [rad]
%
%   flow      : struct con i parametri di flusso:
%               flow.q      = pressione dinamica q_SI [Pa]
%
%   isPowered : boolean
%               true  -> motore acceso (usa A_base,eff)
%               false -> coasting (usa A_base)
%
%   a_sub, b_sub : parametri per wave drag subsonico:
%                  (C_A)_W(M <= 0.8) = a_sub * M^b_sub + (C_A)_{W,M≈0}
%                  Tipicamente scelti per raccordare con il transonico.
%
% OUTPUT:
%   CA : coefficiente assiale del body (stessa size di M)
%
% NOTE:
%   - Valido indicativamente per 0.3 <= M <= 3 per il termine di friction drag.
%   - alpha in radianti.
%   - Tutte le lunghezze e aree in unità SI.

% -----------------------------
% Controllo argomenti opzionali
% -----------------------------
if nargin < 6
    % se non specificati, niente correzione M^b: solo termine M≈0
    a_sub = 0.0;
    b_sub = 1.0;
end

% Per sicurezza, forziamo M e alpha a colonne compatibili
[M, alpha] = compatDims(M, alpha);

% Estrazione parametri geometrici
l     = geom.l;       % lunghezza totale [m]
d     = geom.d;       % diametro [m]
Lnose = geom.Lnose;   % lunghezza nose [m]
Aref  = geom.Aref;    % area di riferimento [m^2]
Anose = geom.Anose;   % area nose [m^2]
Abase = geom.Abase;   % area base [m^2]
phi   = geom.phi;     % angolo di giunzione [rad]

if isfield(geom, "Aexit") && ~isempty(geom.Aexit)
    Aexit_tot = sum(geom.Aexit);
else
    Aexit_tot = 0.0;
end

q = flow.q;           % pressione dinamica [Pa]

% Fineness ratio
lambda = l / d;

% -----------------------------
% 1) Wave drag (C_A)_W
% -----------------------------

% 1.1) Contributi supersonici (M >= 1.3)
CAW_sharp = (1.586 + 1.834 ./ (M.^2)) .* ...
            (atan(0.5 ./ (Lnose ./ d))).^1.69;

CAW_hemi  = 0.665 * (1.586 + 1.834 ./ (M.^2));

CAW_sup   = CAW_sharp .* ((Aref - Anose) ./ Aref) + ...
            CAW_hemi  .* (Anose ./ Aref);

% 1.2) Limite fortemente subsonico (M ≈ 0)
CAW_M0 = 0.8 * sin(phi).^2;   % scalare

% 1.3) Modello subsonico fino a M <= 0.8
CAW_sub = a_sub .* (M.^b_sub) + CAW_M0;

% 1.4) Combinazione piecewise con transonico e supersonico
CAW = zeros(size(M));

% regioni:
%   M <= 0       : poniamo CAW = CAW_M0
%   0 < M <= 0.8 : CAW_sub
%   0.8 < M < 1.3: interpolazione lineare tra CAW_sub(0.8) e CAW_sup(1.3)
%   M >= 1.3     : CAW_sup

% M <= 0
idx = (M <= 0);
CAW(idx) = CAW_M0;

% 0 < M <= 0.8
idx = (M > 0 & M <= 0.8);
CAW(idx) = CAW_sub(idx);

% 0.8 < M < 1.3  -> interpolazione lineare
idx = (M > 0.8 & M < 1.3);
if any(idx)
    % valore a M=0.8 (subsonico)
    M1   = 0.8;
    CA1  = a_sub * (M1^b_sub) + CAW_M0;
    % valore a M=1.3 (supersonico)
    M2   = 1.3;
    CAW_sharp_13 = (1.586 + 1.834 / M2^2) * ...
                   (atan(0.5 / (Lnose / d)))^1.69;
    CAW_hemi_13  = 0.665 * (1.586 + 1.834 / M2^2);
    CA2  = CAW_sharp_13 * ((Aref - Anose) / Aref) + ...
           CAW_hemi_13  * (Anose / Aref);
    % interp lineare in funzione di M
    CAW(idx) = CA1 + (CA2 - CA1) .* (M(idx) - M1) / (M2 - M1);
end

% M >= 1.3
idx = (M >= 1.3);
CAW(idx) = CAW_sup(idx);

% -----------------------------
% 2) Base drag (C_A)_B
% -----------------------------

% Coefficiente base drag riferito all'area di base
CD0B = zeros(size(M));
idx  = (M < 1);
CD0B(idx) = 0.12 + 0.13 .* (M(idx).^2);
idx  = (M >= 1);
CD0B(idx) = 0.25 ./ M(idx);

% Area di base efficace (motore acceso)
Abase_eff = Abase - Aexit_tot;
if Abase_eff < 0
    warning("CA_body: Abase_eff < 0. Controllare A_base e A_exit.");
    Abase_eff = 0;
end

if isPowered
    CA_B = CD0B .* (Abase_eff ./ Aref);
else
    CA_B = CD0B .* (Abase ./ Aref);
end

% -----------------------------
% 3) Friction drag (C_A)_f  (Jerger/Fleeman in SI)
% -----------------------------
% (C_A)_f ≈ 0.091 * (L/d) * ( M / (q * L) )^0.2
% q = pressione dinamica [Pa], L in [m]

CA_f = 0.091 * lambda .* ( M ./ (q * l) ).^0.2;

% -----------------------------
% 4) C_A,alpha=0 e dipendenza in alpha
% -----------------------------
CA0 = CAW + CA_B + CA_f;

% C_A,body(M, alpha) = C_A,alpha=0(M) * cos^2(alpha)
CA = CA0 .* cos(alpha).^2;

end

% ---------------------------------------------------------
function [M_out, alpha_out] = compatDims(M, alpha)
% Rende M e alpha di dimensioni compatibili per operazioni
% element-wise. Se uno dei due è scalare, lo espande alla size dell'altro.

if isscalar(M) && ~isscalar(alpha)
    M_out = M * ones(size(alpha));
    alpha_out = alpha;
elseif ~isscalar(M) && isscalar(alpha)
    M_out = M;
    alpha_out = alpha * ones(size(M));
elseif isscalar(M) && isscalar(alpha)
    M_out = M;
    alpha_out = alpha;
else
    % controllo stessa size
    if ~isequal(size(M), size(alpha))
        error('CA_body: M e alpha devono avere la stessa size oppure essere scalari.');
    end
    M_out = M;
    alpha_out = alpha;
end

end
