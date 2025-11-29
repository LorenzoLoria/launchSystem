function [CL_tot_alpha, CD_tot_mach] = computePlot_CL_CD_tot_2(vector, bodyGeom, finsGeom, bodyInfo, finsInfo, gasProp)
%
% INPUT:
%
% - vector      : struct con range di M e alpha
%                 vector.M                = vettore di Mach [-]
%                 vector.M_cases          = vettore di valori di Mach per CL plot[-]
%                 vector.alpha_deg        = vettore di alpha [deg]
%                 vector.alpha_deg_cases  = vettore di valori di alpha per CD plot [deg]
%
% - bodyGeom    : struct con parametri geometrici del body
%                 bodyGeom.l          = lunghezza totale corpo (ell_SI) [m]
%                 bodyGeom.d          = diametro di riferimento [m]
%                 bodyGeom.Lnose      = lunghezza del nose l_N [m]
%                 bodyGeom.Aref       = area di riferimento A_ref [m^2]
%                 bodyGeom.Anose      = area del nose A_nose [m^2]
%                 bodyGeom.Abase      = area di base A_base [m^2]
%                 bodyGeom.Aexit      = vettore aree ugelli A_exit,nozzle [m^2]
%                 bodyGeom.phi        = angolo di giunzione nose-body φ [rad]
%                 bodyGeom.Ab         = body area caratteristica A_b [m^2]
%                 bodyGeom.Ap         = area proiettata in crossflow A_p [m^2]
%
% - finsGeom    : struct con parametri geometrici delle fins
%                 finsGeom.Nfins      = number of fins [-]
%                 finsGeom.be         = equivalent span of the fin [m]
%                 finsGeom.Se         = surface of the fin [m^2]
%                 finsGeom.cmac       = mean aerodynamic chord of fin [m]
%                 finsGeom.delta_le   = leading edge sweep [deg]
%                 finsGeom.lambda_le  = fin base angle [deg]
%                 finsGeom.b          = 2*length of base of fin [m]
%                 finsGeom.tmac       = max thickness of MAC [m]
%
% - bodyInfo    : struct con info per body CA e CN
%                 bodyInfo.isPowered  = true  --> motore acceso (usa A_base,eff)
%                                       false --> coasting (usa A_base)
%                 bodyInfo.a_sub      = parametro per wave drag subsonico
%                 bodyInfo.b_sub      = parametro per wave drag subsonico
%                 bodyInfo.Cdn        = crossflow drag coefficient C_dn
%
% - finsInfo    : struct con info per fins CA e CN
%                 finsInfo.rho        = density [kg/m^3]
%                 finsInfo.vsound     = speed of sound at time t [m/s]
%
% - gasProp     : struct con info relative al fluido
%                 gasProp.gamma       = gamma
%
% OUTPUT:
%
% - CL_tot_alpha : body+fins CL matrix (rows: alpha_vec -- columns: M_cases)
% - CD_tot_mach  : body+fins CD matrix (rows: M_vec     -- columns: alpha_cases)
%

% DATA
% Vectors in input
M_vec           = vector.M;
M_cases         = vector.M_cases;
alpha_deg_vec   = vector.alpha_deg;
alpha_deg_cases = vector.alpha_deg_cases;   % <-- corretto

% Geometry of the body
l     = bodyGeom.l;
d     = bodyGeom.d;
Lnose = bodyGeom.Lnose;
Aref  = bodyGeom.Aref;
Anose = bodyGeom.Anose;
Abase = bodyGeom.Abase;
Aexit = bodyGeom.Aexit;
phi   = bodyGeom.phi;
Ab    = bodyGeom.Ab;
Ap    = bodyGeom.Ap;

% Geometry of the fins
Nfins     = finsGeom.Nfins;
be        = finsGeom.be;
Se        = finsGeom.Se;
cmac      = finsGeom.cmac;
delta_le  = finsGeom.delta_le;
lambda_le = finsGeom.lambda_le;
b        = finsGeom.b;
tmac     = finsGeom.tmac;

% Info body
isPowered = bodyInfo.isPowered;
a_sub     = bodyInfo.a_sub;
b_sub     = bodyInfo.b_sub;
Cdn       = bodyInfo.Cdn;

% Info fins
rho    = finsInfo.rho;
vsound = finsInfo.vsound;

% Gas properties (per modello Newtoniano ipersonico)
gamma = gasProp.gamma; 

% COMPUTATION
alpha_rad_vec   = deg2rad(alpha_deg_vec);
alpha_rad_cases = deg2rad(alpha_deg_cases);

CL_tot_alpha = zeros(length(alpha_rad_vec),  length(M_cases));
CD_tot_mach  = zeros(length(M_vec),         length(alpha_rad_cases));



% --------- CL_tot_alpha per M_cases ------------------------------------
for i = 1:length(M_cases)

    M = M_cases(i);
    vlauncher = M * vsound;
    q = 0.5 * rho * vlauncher^2;

    for j = 1:length(alpha_rad_vec)
        alpha = alpha_rad_vec(j);

        % --- Body: CN e CA
        aero.Cdn = Cdn;
        flow.q   = q;

        % 1) contributo slender-body + crossflow SEMPRE
        CN_body_val = CN_body(M, alpha, bodyGeom, aero);

        % 2) se M > 5 aggiungo contributo Newtoniano del blunt nose
        if M > 5
            % area windward efficace del nose (~ Aref o Anose)
            S_nose_wind = Aref;
            CN_nose_newt = CN_newton_blunt(alpha, M, gamma, S_nose_wind, Aref);
            CN_body_val  = CN_body_val + CN_nose_newt;
        end

        CA_body_val = CA_body(M, alpha, bodyGeom, flow, isPowered, a_sub, b_sub);

        % --- Fins: CN e CA per UNA ALETTTA
        [CN_fin_single, CD0_fric_single, CD0_wave_single] = ...
            AerodynCoefFins_new(alpha, vlauncher, vsound, be, Se, q, Aref, ...
                                cmac, delta_le, lambda_le, b, tmac);

        CN_fins_tot   = Nfins * CN_fin_single;
        CD0_fins_tot  = Nfins * (CD0_fric_single + CD0_wave_single);
        CA_fins_tot   = CD0_fins_tot * cos(alpha)^2;

        % Somma totali
        CN_tot = CN_body_val + CN_fins_tot;
        CA_tot = CA_body_val + CA_fins_tot;

        % CL / CD totali
        [CL_val, ~] = calculate_CL_CD(CN_tot, CA_tot, rad2deg(alpha));

        CL_tot_alpha(j, i) = CL_val;
    end
end



% --------- CD_tot_mach per alpha_cases --------------------------------
for k = 1:length(alpha_deg_cases)

    alpha = alpha_rad_cases(k);

    for z = 1:length(M_vec)

        M = M_vec(z);
        vlauncher = M * vsound;
        q = 0.5 * rho * vlauncher^2;

        % --- Body
        aero.Cdn = Cdn;
        flow.q   = q;

        % 1) contributo slender-body + crossflow SEMPRE
        CN_body_val = CN_body(M, alpha, bodyGeom, aero);

        % 2) contributo Newtoniano del blunt nose solo per M > 5
        if M > 5
            S_nose_wind = Aref;  
            CN_nose_newt = CN_newton_blunt(alpha, M, gamma, S_nose_wind, Aref);
            CN_body_val  = CN_body_val + CN_nose_newt;
        end

        CA_body_val = CA_body(M, alpha, bodyGeom, flow, isPowered, a_sub, b_sub);

        % --- Fins
        [CN_fin_single, CD0_fric_single, CD0_wave_single] = ...
            AerodynCoefFins_new(alpha, vlauncher, vsound, be, Se, q, Aref, ...
                                cmac, delta_le, lambda_le, b, tmac);

        CN_fins_tot  = Nfins * CN_fin_single;
        CD0_fins_tot = Nfins * (CD0_fric_single + CD0_wave_single);
        CA_fins_tot  = CD0_fins_tot * cos(alpha)^2;

        CN_tot = CN_body_val + CN_fins_tot;
        CA_tot = CA_body_val + CA_fins_tot;

        [~, CD_val] = calculate_CL_CD(CN_tot, CA_tot, rad2deg(alpha));

        CD_tot_mach(z, k) = CD_val;
    end
end



% --------- PLOT --------------------------------------------------------

% CL
figure;
plot(alpha_deg_vec, CL_tot_alpha, 'LineWidth', 1.5);
grid on;
hold on;
xlabel('\alpha [deg]');
ylabel('C_L^{tot}');
title('Total lift coefficient vs \alpha (body + fins): slender body (M <= 5) + blunt nose Newtonian (M > 5)');
leg_str = arrayfun(@(Mval) sprintf('M = %.1f', Mval), M_cases, 'UniformOutput', false);
legend(leg_str, 'Location', 'best');

% CD
figure;
plot(M_vec, CD_tot_mach, 'LineWidth', 1.5);
grid on;
xlabel('Mach number M [-]');
ylabel('C_D^{tot}');
title('Total drag coefficient vs Mach (body + fins)');
leg_str2 = arrayfun(@(aval) sprintf('\\alpha = %.1f^\\circ', aval), alpha_deg_cases, 'UniformOutput', false);
legend(leg_str2, 'Location', 'best');

end



% --------- COMMENTI -----------------------------------------------------
%
% Si combina il contributo "slender–body + crossflow" (Allen–Jørgensen / Fleeman) 
% con un contributo aggiuntivo di tipo Newtoniano per il nose blunt
% ("Newtonian / Modified Newtonian Model for Blunt-Nose Bodies" sul latex in
% appendice)
%
% - Modello slender–body sempre utilizzato come base, in tutti i
%   regimi di Mach (sia M<=5 che M>5), perché garantisce un comportamento 
%   fisicamente corretto per piccoli angoli d’attacco (C_N ~ α) ed è 
%   validato per regimi subsonici, transonici e supersonici moderati (M<=5).
%
% - Per Mach ipersonici (M > 5) si aggiunge un contributo Newtoniano
%   legato al nose blunt, ottenuto tramite il coefficiente di pressione
%   massimo C_{p,max} ricavato dalle relazioni di Rayleigh–Pitot
%   (onda normale + ristagno) --> termine che rappresenta il carico normale
%   aggiuntivo dovuto allo shock forte sul nose a Mach elevati.
%
% - Il modello Newtoniano è sommato al modello slender–bod --> così si 
%   mantiene una pendenza C_{N,α} realistica vicino ad α = 0 e allo stesso 
%   tempo si ha l’incremento di pressione tipico del regime ipersonico.
%
% - I contributi delle alette (C_N e C_A) vengono sempre aggiunti al totale,
%   indipendentemente dal regime di Mach.
%
%--------------------------------------------------------------------------

