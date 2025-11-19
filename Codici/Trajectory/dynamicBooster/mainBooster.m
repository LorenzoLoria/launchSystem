%% Code to Recover the Booster

addpath(genpath("..\..\"))
mission = dataStruct;

%%

% Stato iniziale
r0  = [0; 0; 6500000];        % posizione iniziale
v0  = [5000; 0; 0];      % velocità iniziale

% Orientazione iniziale (identità): nessuna rotazione rispetto all'inerziale
q0 = [1; 0; 0; 0];

% Velocità angolare iniziale (ad es. spin attorno all'asse y del corpo)
w0 = [0; 0 ; 0];        % rad/s

% Stato iniziale completo
y0 = [r0; v0; q0; w0];

% Intervallo temporale
tspan = [0 10000];

% Integrazione
windDirection =[1/sqrt(2) 1/sqrt(2) 0]';
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@groundEvent);
f = @(t,x) dynBooster(t, x,mission, windDirection);
[t_sol, y_sol] = ode113(f, tspan, y0,options);

%%

plot(t_sol, y_sol(:,11:13))
legend(["OmegaX" "OmegaY" "OmegaZ"])


%%

t = t_sol;
q_hist = y_sol(:,7:10);      % quaternioni
r_hist = y_sol(:,1:3);       % posizioni

N = length(t);

% ... (Il tuo codice precedente per t, q_hist, r_hist, N rimane uguale) ...

% Ricalcolo assi inerziali per sicurezza (se non lo hai già fatto sopra)
ex_I = zeros(N,3); ey_I = zeros(N,3); ez_I = zeros(N,3);
for k = 1:N
    q = q_hist(k,:).'; q = q / norm(q);
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    Rbi = [ 1-2*(q2^2+q3^2),   2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2);
            2*(q1*q2 + q0*q3), 1-2*(q1^2+q3^2),   2*(q2*q3 - q0*q1);
            2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1-2*(q1^2+q2^2) ];
    ex_I(k,:) = Rbi(:,1).'; ey_I(k,:) = Rbi(:,2).'; ez_I(k,:) = Rbi(:,3).';
end

figure; hold on; grid on; axis equal;
xlabel('X_I [km]'); ylabel('Y_I [km]'); zlabel('Z_I [km]');
title('Animazione Booster: Traiettoria e Assetto');
%EarthPlot(6371)
view(45,15)
%xlim([-500 1500])
%ylim([-700 700])
%zlim([6000 6700])

% --- 1. SETUP DEL CILINDRO (SYSTEM BODY) ---
% Definisci dimensioni visibili sul grafico (in km se la traiettoria è in km)
r_cyl = 5;   % Raggio cilindro (es. 1 km per vederlo bene)
h_cyl = 25;   % Altezza cilindro
n_res = 20;    % Risoluzione (bassa per rendere l'animazione fluida)

% Genera cilindro standard (allineato su Z per default)
[Xc_base, Yc_base, Zc_base] = cylinder(r_cyl, n_res);
Zc_base = Zc_base * h_cyl; 
% Centro il cilindro sull'origine (opzionale, dipende da dove è il CM)
Zc_base = Zc_base - h_cyl/2; 

% IMPORTANTE: Ruoto i punti "base" per allineare il cilindro all'asse X (Roll axis)
% Se vuoi che il cilindro segua l'asse Z (blu), commenta queste 3 righe e usa X=Xc, Y=Yc, Z=Zc
Xc_body = Xc_base;  % La lunghezza ora è su X
Yc_body = Yc_base;
Zc_body = Zc_base;  % Il raggio è su Y e Z

% Creiamo l'oggetto grafico 'mesh' iniziale (vuoto o all'origine)
% Usiamo 'EdgeColor' definito e 'FaceColor' none come richiesto (outline)
hCylMesh = mesh(Xc_body, Yc_body, Zc_body, 'FaceColor', 'none', 'EdgeColor', 'y', 'FaceAlpha', 0);

% --- SETUP ASSI E TRAIETTORIA ---
plot3(r_hist(:,1)/1000, r_hist(:,2)/1000, r_hist(:,3)/1000, 'r-', 'LineWidth', 1.0);
L = 10; % Lunghezza frecce assi (scalata per visibilità)
hX = quiver3(0,0,0, 0,0,0, 'r', 'LineWidth', 2); 
hY = quiver3(0,0,0, 0,0,0, 'g', 'LineWidth', 2); 
hZ = quiver3(0,0,0, 0,0,0, 'b', 'LineWidth', 2); 



% --- 2. LOOP DI ANIMAZIONE ---
step = max(1, floor(N/2000)); % Step per circa 500 frames
for k = 1:step:N
    p = r_hist(k,:)/1000; % Posizione attuale in km

    % Ricalcolo Rbi locale (o lo prendi da un array salvato)
    q = q_hist(k,:).'; q = q / norm(q);
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    Rbi = [ 1-2*(q2^2+q3^2),   2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2);
            2*(q1*q2 + q0*q3), 1-2*(q1^2+q3^2),   2*(q2*q3 - q0*q1);
            2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1-2*(q1^2+q2^2) ];

    % --- AGGIORNAMENTO CILINDRO ---
    % Per ruotare la mesh, dobbiamo operare su tutti i punti come vettori
    % Le matrici Xc_body sono 2x(n_res+1). Le linearizziamo per la rotazione.
    dims = size(Xc_body);
    P_body = [Xc_body(:)'; Yc_body(:)'; Zc_body(:)']; % Matrice 3xM

    % Rotazione + Traslazione
    % P_inerziale = Rbi * P_body + Posizione
    P_inerz = Rbi * P_body + p.';

    % Rimettiamo i dati nella forma di matrice per la funzione mesh
    X_new = reshape(P_inerz(1,:), dims);
    Y_new = reshape(P_inerz(2,:), dims);
    Z_new = reshape(P_inerz(3,:), dims);

    % Aggiorna la grafica del cilindro
    set(hCylMesh, 'XData', X_new, 'YData', Y_new, 'ZData', Z_new);

    % --- AGGIORNAMENTO ASSI (QUIVER) ---
    ex = Rbi(:,1); ey = Rbi(:,2); ez = Rbi(:,3);
    set(hX, 'XData', p(1), 'YData', p(2), 'ZData', p(3), 'UData', L*ex(1), 'VData', L*ex(2), 'WData', L*ex(3));
    set(hY, 'XData', p(1), 'YData', p(2), 'ZData', p(3), 'UData', L*ey(1), 'VData', L*ey(2), 'WData', L*ey(3));
    set(hZ, 'XData', p(1), 'YData', p(2), 'ZData', p(3), 'UData', L*ez(1), 'VData', L*ez(2), 'WData', L*ez(3));

    % Aggiorna limiti assi per seguire il booster (Camera Follow) - Opzionale
    % xlim([p(1)-20 p(1)+20]); ylim([p(2)-20 p(2)+20]); zlim([p(3)-20 p(3)+20]);

    drawnow limitrate; % limitrate rende l'animazione più fluida scartando frame se necessario
end
legend([hX hY hZ], 'X_b', 'Y_b', 'Z_b');
