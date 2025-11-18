% Parametri fisici del cilindro e dell'ambiente
params.m   = 1.0;       % kg
params.Cd  = 0.8;       % coeff. drag traslazionale
params.A   = 0.01;      % m^2 (area di riferimento)
params.rho = 1.225;     % kg/m^3
params.g   = 9.81;      % m/s^2

% Geometria del cilindro
R = 0.05;   % raggio [m]
L = 0.30;   % lunghezza [m]

% Momenti di inerzia di un cilindro pieno, asse z del corpo allineato all'asse del cilindro
Jx = 1/12 * params.m * (3*R^2 + L^2);
Jy = Jx;
Jz = 0.5  * params.m * R^2;
params.J = diag([Jx Jy Jz]);

% Smorzamento rotazionale (puoi regolare o mettere 0)
params.C_omega = 0.01;  % N m s

% Stato iniziale
r0  = [0; 0; 6500000];        % posizione iniziale
v0  = [3000; 0; 0];      % velocità iniziale

% Orientazione iniziale (identità): nessuna rotazione rispetto all'inerziale
q0 = [1; 0; 0; 0];

% Velocità angolare iniziale (ad es. spin attorno all'asse y del corpo)
w0 = [0; 0 ; 0];        % rad/s

% Stato iniziale completo
y0 = [r0; v0; q0;w0];

% Intervallo temporale
tspan = [0 10000];

% Integrazione
windDirection =[0 0 1]';
options = odeset('RelTol',1e-6,'AbsTol',1e-2,'Events',@groundEvent);
f = @(t,x) dynBooster(t, x,mission, windDirection);
tic
for i= 1:100
[t_sol, y_sol] = ode113(f, tspan, y0,options);
end
toc
%%
plot(t_sol, y_sol(:,11:13))
legend(["OmegaX" "OmegaY" "OmegaZ"])






%%

t = t_sol;
q_hist = y_sol(:,7:10);      % quaternioni
r_hist = y_sol(:,1:3);       % posizioni

N = length(t);

ex_I = zeros(N,3);  % asse X_body in inerziale
ey_I = zeros(N,3);  % asse Y_body in inerziale
ez_I = zeros(N,3);  % asse Z_body in inerziale

for k = 1:N
    q = q_hist(k,:).';
    q = q / norm(q);

    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);

    % Rbi: rotazione body -> inerziale (stessa formula usata in dynBooster)
    Rbi = [ 1-2*(q2^2+q3^2),   2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2);
            2*(q1*q2 + q0*q3), 1-2*(q1^2+q3^2),   2*(q2*q3 - q0*q1);
            2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1-2*(q1^2+q2^2) ];

    % Le colonne di Rbi sono gli assi body espressi nell'inerziale
    ex_I(k,:) = Rbi(:,1).';   % X_body
    ey_I(k,:) = Rbi(:,2).';   % Y_body
    ez_I(k,:) = Rbi(:,3).';   % Z_body
end




figure; hold on; grid on;
axis equal;

% Traiettoria
plot3(r_hist(:,1)/1000, r_hist(:,2)/1000, r_hist(:,3)/1000, 'r-', 'LineWidth', 2.0);
hold on
EarthPlot(6371)
xlabel('X_I [km]');
ylabel('Y_I [km]');
zlabel('Z_I [km]');
title('Terna del booster nel sistema inerziale');
view(0,25)
% Lunghezza degli assi da plottare
L = 100;   % scegli tu la scala (m)

% Inizializza le tre frecce per Xb, Yb, Zb
hX = quiver3(0,0,0, 0,0,0, 'r', 'LineWidth', 2); % X_body (rosso)
hY = quiver3(0,0,0, 0,0,0, 'g', 'LineWidth', 2); % Y_body (verde)
hZ = quiver3(0,0,0, 0,0,0, 'b', 'LineWidth', 2); % Z_body (blu)

% Limiti del grafico (ad es. dal minimo al massimo della traiettoria)
%xlim([min(r_hist(:,1)) max(r_hist(:,1))]);
%ylim([min(r_hist(:,2)) max(r_hist(:,2))]);
%zlim([min(r_hist(:,3)) max(r_hist(:,3))]);

% Loop di animazione (puoi cambiare step per avere più/meno frame)
step = max(1, floor(N/2000));  % al massimo ~200 frame

for k = 1:step:N
    p  = r_hist(k,:)/1000;       % posizione del booster
    ex = ex_I(k,:);         % asse X_body in inerziale
    ey = ey_I(k,:);         % asse Y_body in inerziale
    ez = ez_I(k,:);         % asse Z_body in inerziale

    % Aggiorna le frecce:
    set(hX, 'XData', p(1), 'YData', p(2), 'ZData', p(3), ...
            'UData', L*ex(1), 'VData', L*ex(2), 'WData', L*ex(3));

    set(hY, 'XData', p(1), 'YData', p(2), 'ZData', p(3), ...
            'UData', L*ey(1), 'VData', L*ey(2), 'WData', L*ey(3));

    set(hZ, 'XData', p(1), 'YData', p(2), 'ZData', p(3), ...
            'UData', L*ez(1), 'VData', L*ez(2), 'WData', L*ez(3));

    drawnow;
end

legend('Traiettoria','X_b','Y_b','Z_b');
