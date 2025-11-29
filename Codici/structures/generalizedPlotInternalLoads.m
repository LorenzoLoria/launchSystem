clear
clc
close all

% ============================== DATA =====================================
nComponents = 8;
nPointsPerComponent = 100;

m1 = 180e3;
m2 = 40e3;
mp = 8e3;

launcher.mass = [mp, m2, m1];
launcher.drag = 350000;
launcher.gamma = deg2rad(2);
launcher.alpha = deg2rad(3.5);
launcher.accelerationAxial = 1.6 * 9.81;
launcher.lift = 0;
launcher.g0 = 9.81;
launcher.accelerationNormal = 2 * 9.81;
launcher.elementLength = [4, 1, 2, 4, 4, 15, 2, 13]; % sono le lunghezze degli elementi

% ========================== SOLUTION =====================================

[N, T, M, A] =  loadsFinder(nComponents, launcher);

%% ============================== PLOTS ====================================

x_all = [];
N_all = [];
T_all = [];
M_all = [];

x_coordinates = cumsum([0, launcher.elementLength]); % defines the coordinates of
% start and finish of each component

% Required for interpolation
M_end_values = [M(2:end); 0];

for i = 1:nComponents
    
    x_start = x_coordinates(i);
    x_end   = x_coordinates(i+1);
    x_current = linspace(x_start, x_end, nPointsPerComponent);
    
    % Axial Load
    N_current = ones(1, nPointsPerComponent) * N(i);
    
    % Shear Load
    T_current = ones(1, nPointsPerComponent) * T(i);
    
    % Bending Moment (changes linearly)
    M_start = M(i);
    M_end   = M_end_values(i);
    M_current = linspace(M_start, M_end, nPointsPerComponent);
    
    x_all = [x_all, x_current];
    N_all = [N_all, N_current];
    T_all = [T_all, T_current];
    M_all = [M_all, M_current];
end

figure(1); 
plot(x_all, N_all, 'LineWidth', 1.5, 'Color', 'b');
hold on
yline(0, 'LineWidth', 1.5)
grid on;
xlabel('x [m]');
ylabel('Axial Load [N]');
xlim([0,x_all(end)])

figure(2); 
plot(x_all, T_all, 'LineWidth', 1.5, 'Color', 'b');
hold on
yline(0, 'LineWidth', 1.5)
grid on;
xlabel('x [m]');
ylabel('Shear Load [N]');
xlim([0,x_all(end)])

figure; 
plot(x_all, M_all, 'LineWidth', 1.5, 'Color', 'b');
hold on
yline(0, 'LineWidth', 1)
grid on;
xlabel('x [m]');
ylabel('Bending Moment [Nm]');
xlim([0,x_all(end)])