function CN = CN_body(M, alpha, geom, aero)
% CN_body  Calcolo del coefficiente di forza normale del body C_N,body
%          usando il modello "slender-body + crossflow" (Allen/Jørgensen).
%
%   CN = CN_body(M, alpha, geom, aero)
%
% INPUT:
%   M       : numero di Mach (scalare o vettore)
%   alpha   : angolo d'attacco [rad]
%             - scalare -> applicato a tutti i Mach
%             - vettore -> deve avere stesso numero di elementi di M
%
%   geom    : struct con parametri geometrici:
%             geom.Aref = area di riferimento A [m^2]
%             geom.Ab   = body area caratteristica A_b [m^2]
%             geom.Ap   = area proiettata in crossflow A_p [m^2]
%             geom.l    = lunghezza caratteristica ell [m]
%             geom.d    = diametro caratteristico d [m]
%
%   aero    : struct con parametri aerodinamici:
%             aero.Cdn  = crossflow drag coefficient C_dn
%                         - scalare -> stesso per tutti i Mach
%                         - vettore -> stessa size di M
%
% OUTPUT:
%   CN : coefficiente di forza normale del body (stessa size di M)
%
% MODELLO:
%   C_N(alpha) = (A_b/A) sin(2 alpha) cos(alpha/2)
%                + eta * C_dn * (A_p/A) sin^2(alpha)
%
%   con
%     eta =
%       0.05 (ell/d) + 0.52        M <= 1  (subsonico)
%       1                          M > 1   (supersonico)
%
% NOTE:
%   - alpha in radianti.
%   - l/d (slenderness) tipicamente >= 5.

% Salvo la forma di M per ripristinarla in uscita
sizeM = size(M);
Mvec  = M(:);       % colonna

% -----------------------------
% Gestione di alpha (scalare o vettore)
% -----------------------------
if isscalar(alpha)
    alphavec = alpha * ones(size(Mvec));
else
    if numel(alpha) ~= numel(M)
        error('CN_body: alpha deve essere scalare oppure avere lo stesso numero di elementi di M.');
    end
    alphavec = alpha(:);
end

% -----------------------------
% Geometria
% -----------------------------
Aref = geom.Aref;
Ab   = geom.Ab;
Ap   = geom.Ap;
ell  = geom.l;
d    = geom.d;

Ab_over_A = Ab / Aref;
Ap_over_A = Ap / Aref;

lambda = ell / d;                   % slenderness ratio ell/d

% -----------------------------
% Crossflow drag coefficient C_dn
% -----------------------------
if ~isfield(aero, "Cdn")
    error('CN_body: aero.Cdn mancante (crossflow drag coefficient).');
end

Cdn = aero.Cdn;

if isscalar(Cdn)
    Cdnvec = Cdn * ones(size(Mvec));
else
    if numel(Cdn) ~= numel(M)
        error('CN_body: aero.Cdn deve essere scalare oppure avere la stessa size di M.');
    end
    Cdnvec = Cdn(:);
end

% -----------------------------
% Fattore di scala eta (subsonico/supersonico)
% -----------------------------
eta_sub = 0.05 * lambda + 0.52;     % valore subsonico

eta = ones(size(Mvec));
idx_sub = (Mvec <= 1);
eta(idx_sub) = eta_sub;

% -----------------------------
% Calcolo dei termini slender-body e crossflow
% -----------------------------
CN_slender   = Ab_over_A .* sin(2 .* alphavec) .* cos(alphavec ./ 2);
CN_crossflow = eta .* Cdnvec .* Ap_over_A .* (sin(alphavec).^2);

% Somma totale
CN_vec = CN_slender + CN_crossflow;

% Ripristino la shape di M in uscita
CN = reshape(CN_vec, sizeM);

end
