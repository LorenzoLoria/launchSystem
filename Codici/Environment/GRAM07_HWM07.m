function v_wind = GRAM07_HWM07(h, season)
%
% INPUT
% - h      : altitude [0-100 km]
% - season : options: autumn (1), winter (2), spring (3), summer (4)
%
% OUTPUT
% - v_wind : wind speed [m/s]
%
%
% NOTE:
%   - Modello climatologico stagionale 0–100 km ispirato a GRAM07/HWM07
%     Ref: https://ntrs.nasa.gov/api/citations/20090004463/downloads/20090004463.pdf
%     Ref: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008JA013668?utm_source=chatgpt.com
%
%   - 0–2 km: power-law nel boundary layer (alpha ~ 0.15, moderately rough terrain) [ref: Spera, D. A. (2009). Wind Turbine Technology. 2nd ed. ASME Press]
%   - 2–16 km: jet troposferico, modellato con parabola che interpola (2 km, V_BL(2)), (z_j, V_j), (16 km, V_16)
%   - 16–50 km: stratosfera, decadimento esponenziale tra V_16 e V_50
%   - 50–80 km: mesosfera bassa, interpolazione lineare tra V_50 e V_80
%   - 80–100 km: transizione/plateau verso V_100
%


% -------------------------- Parametri stagionali -------------------------

alpha      = 0.15;  % gradiente verticale nel boundary layer (power-law)
zj         = 11.0;  % quota del jet [km]
H_strat    = 10.0;  % scala di decadimento stratosfera [km]
H_upper    = 10.0;  % scala transizione 80–100 km [km]

switch season
    case 1 % autumn (MAM) 
        v10   = 10;   % [m/s] vento a 10 m
        Vj    = 50;   % [m/s] massimo del jet
        V16   = 30;   % [m/s] a 16 km
        V50   = 18;   % [m/s] a 50 km
        V80   = 30;   % [m/s] a 80 km
        V100  = 35;   % [m/s] a 100 km

    case 2 % winter (JJA) 
        v10   = 8;    % [m/s]
        Vj    = 60;   % [m/s] jet più intenso in inverno
        V16   = 40;   % [m/s]
        V50   = 25;   % [m/s]
        V80   = 40;   % [m/s]
        V100  = 45;   % [m/s]

    case 3 % spring (SON) 
        v10   = 9;    % [m/s]
        Vj    = 55;   % [m/s]
        V16   = 35;   % [m/s]
        V50   = 20;   % [m/s]
        V80   = 35;   % [m/s]
        V100  = 40;   % [m/s]

    case 4 % summer (DJF)
        v10   = 7;    % [m/s]
        Vj    = 38;   % [m/s] jet più debole in estate
        V16   = 25;   % [m/s]
        V50   = 15;   % [m/s]
        V80   = 25;   % [m/s]
        V100  = 30;   % [m/s]

    otherwise
        error('season must be 1 (autumn), 2 (winter), 3 (spring), or 4 (summer).');
end

% Valore BL a 2 km (per continuità nel tratto 2–16 km)
v_bl_2 = v10 * (2/10)^alpha;

% Coefficienti della parabola nel tratto 2–16 km:
% V(z) = A*z^2 + B*z + C che interpola: (2 km, v_bl_2), (zj, Vj), (16 km, V16)
z1 = 2;   y1 = v_bl_2;
z2 = zj;  y2 = Vj;
z3 = 16;  y3 = V16;
M  = [z1^2 z1 1;
      z2^2 z2 1;
      z3^2 z3 1];
Y  = [y1; y2; y3];
abc = M \ Y;
A   = abc(1);
B   = abc(2);
C   = abc(3);


% ------------------- Calcolo v_wind in funzione di h ---------------------

if h == 0
    % vento al suolo (idealizzato nullo in questo modello)
    v_wind = 0;

elseif h > 0 && h <= 2
    % Boundary layer (0–2 km) – power-law
    v_wind = v10 * (h / 10)^alpha;

elseif h > 2 && h <= 16
    % Troposfera media–alta / jet – parabola
    v_wind = A*h^2 + B*h + C;

elseif h > 16 && h <= 50
    % Stratosfera (16–50 km) – decadimento esponenziale da V16 a V50
    v_wind = V16 * exp(-(h - 16)/H_strat) + V50 * (1 - exp(-(h - 16)/H_strat));

elseif h > 50 && h <= 80
    % Mesosfera bassa (50–80 km) – interpolazione lineare
    v_wind = V50 + (V80 - V50) * (h - 50) / (80 - 50);

elseif h > 80 && h <= 100
    % 80–100 km – transizione/plateau verso V100
    v_wind = V80 + (V100 - V80) * (1 - exp(-(h - 80)/H_upper));

else
    % scatta solo se h>100
    v_wind = NaN;
end

end
