clear; clc; close all;

nComponents = 4;
nPointsPerComponent = 100;

launcher.mass = [8e3, 20e3, 1000, 180e3];
launcher.drag = ones(nComponents,1) * 5000;
launcher.alpha = deg2rad(2);
launcher.theta = deg2rad(2);
launcher.thrust = 845e3 * 4;
launcher.acceleration = launcher.thrust / sum(launcher.mass(:));
launcher.lift = ones(nComponents,1) * 0;
launcher.g0 = 9.81;

launcher.stagesDimensions = [10, 15, 2, 30];

[b, loadsResults] = loadsFinder_freefree(nComponents, launcher);

nNodes = nComponents + 1;
N = loadsResults(1:3:end);
T = loadsResults(2:3:end);
M = loadsResults(3:3:end);

% ----- Coordinates -----
x_coordinates = cumsum([0, launcher.stagesDimensions]);
nPointsPerComponent = 100;
x_all = []; N_all = []; T_all = []; M_all = [];

for i = 1:nComponents
    x_start = x_coordinates(i);
    x_end   = x_coordinates(i+1);
    x_current = linspace(x_start, x_end, nPointsPerComponent);
    % Axial Load
    N_current = ones(1, nPointsPerComponent) * N(i);
    
    % Shear Load
    T_current = ones(1, nPointsPerComponent) * T(i);
    
    % Bending Moment (changes linearly)
    M_current = linspace(M(i), M(i+1), nPointsPerComponent);
    
    x_all = [x_all, x_current];
    N_all = [N_all, N_current];
    T_all = [T_all, T_current];
    M_all = [M_all, M_current];
end

% ============================== PLOTS ====================================

% Costruisco versioni "spezzate" con NaN tra un elemento e l'altro
xN_plot = []; N_plot = [];
xT_plot = []; T_plot = [];
xM_plot = []; M_plot = [];

for i = 1:nComponents
    i_start = (i-1)*nPointsPerComponent + 1;
    i_end   = i*nPointsPerComponent;

    % Segmento i-esimo
    x_segN = x_all(i_start:i_end);
    N_seg  = N_all(i_start:i_end);

    x_segT = x_all(i_start:i_end);
    T_seg  = T_all(i_start:i_end);

    x_segM = x_all(i_start:i_end);
    M_seg  = M_all(i_start:i_end);

    % Aggiungo il segmento + un NaN per spezzare la linea
    xN_plot = [xN_plot, x_segN, NaN];
    N_plot  = [N_plot,  N_seg,  NaN];

    xT_plot = [xT_plot, x_segT, NaN];
    T_plot  = [T_plot,  T_seg,  NaN];

    xM_plot = [xM_plot, x_segM, NaN];
    M_plot  = [M_plot,  M_seg,  NaN];
end

% ----- Axial Load -----
figure;
hN = area(xN_plot, N_plot*1e-6);
hN.EdgeColor = 'b';           % bordo blu
hN.LineWidth = 1.5;
hN.FaceColor = [0.8 0.9 1];   % azzurrino
hold on;
yline(0,'k','LineWidth',1.0);
grid on;
xlabel('x [m]');
ylabel('Axial Load N [MN]');
title('Axial Load');
xlim([min(x_all), max(x_all)]);

% ----- Shear Load -----
figure;
hT = area(xT_plot, T_plot*1e-6);
hT.EdgeColor = 'b';
hT.LineWidth = 1.5;
hT.FaceColor = [0.8 0.9 1];
hold on;
yline(0,'k','LineWidth',1.0);
grid on;
xlabel('x [m]');
ylabel('Shear Load T [MN]');
title('Shear Load');
xlim([min(x_all), max(x_all)]);

% ----- Bending Moment -----
figure;
hM = area(xM_plot, M_plot*1e-6);
hM.EdgeColor = 'b';
hM.LineWidth = 1.5;
hM.FaceColor = [0.8 0.9 1];
hold on;
yline(0,'k','LineWidth',1.0);
grid on;
xlabel('x [m]');
ylabel('Bending Moment M [MNm]');
title('Bending Moment');
xlim([min(x_all), max(x_all)]);

