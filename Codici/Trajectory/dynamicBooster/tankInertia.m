function [Iz, Ix] = tankInertia(M, R, H)
% Calcola i momenti d'inerzia (Iz e Ix=Iy) di una cisterna sottile omogenea.
%
% INPUT:
% M: Massa totale del guscio sottile.
% R: Raggio del cilindro e dei tappi sferici.
% H: Altezza (lunghezza) della sezione cilindrica.
%
% OUTPUT:
% Iz: Momento d'inerzia rispetto all'asse di simmetria (z).
% Ix: Momento d'inerzia rispetto all'asse trasversale (x o y).
%
% NOTE: Assume che il CM sia al centro geometrico e che la massa sia distribuita
%       in base all'area superficiale.

% --- 1. Definizione delle masse componenti in base all'area ---
K = H + 2*R; % Denominatore comune (proporzionale all'Area Totale / (2*pi*R))

% Massa del Cilindro (Mc) e di un Tappo Sferico (Ms)
Mc = M * (H / K);
Ms = M * (R / K);

% --- 2. Calcolo di Iz (Asse di Simmetria) ---
% Iz = ( Mc + 4/3 * Ms ) * R^2
Iz = (Mc + (4/3) * Ms) * R^2;


% --- 3. Calcolo di Ix = Iy (Assi Trasversali) ---
% La formula usa il Teorema di Huygens-Steiner:
% Ix = I_cilindro,x + 2 * I_tappo,x

% A. Contributo del Cilindro (rispetto al suo CM)
I_cilindro_x = Mc * ( (R^2/2) + (H^2/12) );

% B. Contributo dei due Tappi Sferici
% Distanza d tra il CM del tappo e il CM della cisterna
d = (H + R) / 2;

% I_tappo_x = I_cm,tappo + Ms * d^2
I_cm_tappo = (5/12) * Ms * R^2;
I_tappo_x = I_cm_tappo + Ms * (d^2);

% Somma finale per Ix
Ix = I_cilindro_x + 2 * I_tappo_x;

end