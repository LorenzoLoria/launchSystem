clear
clc

x = [30, 2, 10, 3];
nComponents = 4;
launcher.mass = ones(nComponents,1) * 50e3;
launcher.drag = ones(nComponents,1) * 50000;
launcher.acceleration = 4 * 9.81;
launcher.alpha = deg2rad(45);
launcher.theta = deg2rad(4);
launcher.thrust = 4*845e3;
launcher.lift= ones(nComponents,1) * 50000;
launcher.g0 = 9.81;

loads = loadsFinder(4, x, launcher);
N = loads(1:3:end);
T = loads(2:3:end);
M = loads(3:3:end);

x = x(:);

x1 = linspace(0, x(1), 100);
x2 = linspace(x(1), x(1)+x(2), 100);
x3 = linspace(x(1)+x(2), x(1)+x(2)+x(3), 100);
x4 = linspace(x(1)+x(2)+x(3), x(1)+x(2)+x(3)+x(4), 100);

N1 = ones(100,1) * N(1);
N2 = ones(100,1) * N(2);
N3 = ones(100,1) * N(3);
N4 = ones(100,1) * N(4);

T1 = ones(100,1) * T(1);
T2 = ones(100,1) * T(2);
T3 = ones(100,1) * T(3);
T4 = ones(100,1) * T(4);

M1 = ones(100,1) * M(1);
M2 = ones(100,1) * M(2);
M3 = ones(100,1) * M(3);
M4 = ones(100,1) * M(4);

x_all = [x1, x2, x3, x4];
N_all = [N1.', N2.', N3.', N4.'];   % qui metto il .' per avere riga
T_all = [T1.', T2.', T3.', T4.'];
M_all = [M1.', M2.', M3.', M4.'];

% ============================== PLOTS ====================================
figure; plot(x_all, N_all, 'LineWidth', 1.5);
grid on;
xlabel('x');
ylabel('N');
title('Andamento N(x)');

figure; plot(x_all, T_all, 'LineWidth', 1.5);
grid on;
xlabel('x');
ylabel('T');
title('Andamento T(x)');

figure; plot(x_all, M_all, 'LineWidth', 1.5);
grid on;
xlabel('x');
ylabel('M');
title('Andamento M(x)');

% s_nodes = [0; cumsum(x)];                     % nodi globali
% s_elem  = (s_nodes(1:end-1) + s_nodes(2:end))/2;
% 
% figure; plot(s_elem, N, '-o'); xlabel('s [m]'); ylabel('N');
% figure; plot(s_elem, T, '-o'); xlabel('s [m]'); ylabel('T');
% figure; plot(s_elem, M, '-o'); xlabel('s [m]'); ylabel('M');
%%
nPoints = 20;  % punti interni per rappresentare l’elemento
for i = 1:nComponents
    xi_local = linspace(0, x(i), nPoints);  % da 0 alla lunghezza dell’elemento i

    % Se per ora consideri N e T costanti dentro l’elemento:
    Ni = N(i) * ones(size(xi_local));
    Ti = T(i) * ones(size(xi_local));

    % M_i(x) puoi ricostruirlo, ad esempio lineare, se hai T costante
    % Esempio: momento lineare lungo l’elemento con M(i) al nodo di sinistra
    % e M(i+1) a destra (qui ci serve M(i+1), con un check sull’ultimo elemento).
    if i < nComponents
        Mi_local = M(i) + (M(i+1) - M(i)) * (xi_local / x(i));
    else
        Mi_local = M(i) * ones(size(xi_local));  % oppure una BC tua
    end

    figure(1); hold on; plot(xi_local, Ni); xlabel('x_{local} [m]'); ylabel('N');
    figure(2); hold on; plot(xi_local, Ti); xlabel('x_{local} [m]'); ylabel('T');
    figure(3); hold on; plot(xi_local, Mi_local); xlabel('x_{local} [m]'); ylabel('M');
end
