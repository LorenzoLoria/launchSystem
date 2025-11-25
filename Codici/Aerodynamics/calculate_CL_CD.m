function [C_L, C_D] = calculate_CL_CD(C_N, C_A, alpha_deg)
% Calcola i coefficienti di lift e drag a partire dai coefficienti
% di forza normale (C_N) e assiale (C_A), dati un angolo di attacco
% alpha in gradi.

% Converti l'angolo di attacco in radianti
alpha_rad = deg2rad(alpha_deg);

% Calcola il coefficiente di lift (C_L)
C_L = C_N * cos(alpha_rad) - C_A * sin(alpha_rad);

% Calcola il coefficiente di drag (C_D)
C_D = C_A * cos(alpha_rad) + C_N * sin(alpha_rad);

end
