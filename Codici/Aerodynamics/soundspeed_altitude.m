function a = soundspeed_altitude(p)
% 
% Calcola e plotta la velocità del suono al variare della quota tra 0 e 50
% km, con passo "p" dato in input --> h_range = 0 : p : 50000
%
% Temperatura secondo atmosfera standard (ISA)
% Strati usati:
%   --> 0 –11 km:    gradiente -6.5 K/km
%   --> 11–20 km:    isoterma 216.65 K
%   --> 20–32 km:    gradiente +1   K/km
%   --> 32–47 km:    gradiente +2.8 K/km
%   --> 47–50 km:    isoterma 270.65 K (parte della fascia 47–51 km)


% Costanti fisiche
gamma = 1.4;          % [-] rapporto dei calori (aria)
R     = 287.05;       % [J/(kg·K)] costante specifica dell'aria

% Quota
h_range = (0:p:50000)';   % [m] da 0 a 50 km con passo "p"

% Temperatura secondo atmosfera standard (ISA)
T = zeros(size(h_range));

idx1 = h_range <= 11000;                           % 0 –11 km
idx2 = h_range > 11000 & h_range <= 20000;         % 11–20 km
idx3 = h_range > 20000 & h_range <= 32000;         % 20–32 km
idx4 = h_range > 32000 & h_range <= 47000;         % 32–47 km
idx5 = h_range > 47000 & h_range <= 50000;         % 47–50 km

% Strato 1: 0–11 km
T(idx1) = 288.15 - 0.0065 * h_range(idx1);

% Strato 2: 11–20 km (isotermo)
T(idx2) = 216.65;

% Strato 3: 20–32 km, +1 K/km
T(idx3) = 216.65 + 0.001 * (h_range(idx3) - 20000);

% Strato 4: 32–47 km, +2.8 K/km
T(idx4) = 228.65 + 0.0028 * (h_range(idx4) - 32000);

% Strato 5: 47–50 km (isotermo)
T(idx5) = 270.65;

% Velocità del suono: a = sqrt(gamma * R * T)
a = sqrt(gamma * R .* T);    % [m/s]

end