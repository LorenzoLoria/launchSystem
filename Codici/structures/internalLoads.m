clear; clc; close all;

nComponents = 4;
nPointsPerComponent = 100;

launcher.mass = [180e3, 1000, 20e3, 8e3];
launcher.drag = ones(nComponents,1) * 5000;
launcher.alpha = deg2rad(2);
launcher.theta = deg2rad(4);
launcher.thrust = 4*845e3;
launcher.acceleration = launcher.thrust / launcher.mass(1);
launcher.lift = ones(nComponents,1) * 0;
launcher.g0 = 9.81;

launcher.stagesDimensions = [30, 2, 15, 10];

loadsResults = loadsFinder_freefree(nComponents, launcher.stagesDimensions(:), launcher);

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
    N_current = linspace(N(i), N(i+1), nPointsPerComponent);
    T_current = linspace(T(i), T(i+1), nPointsPerComponent);
    M_current = linspace(M(i), M(i+1), nPointsPerComponent);
    x_all = [x_all, x_current];
    N_all = [N_all, N_current];
    T_all = [T_all, T_current];
    M_all = [M_all, M_current];
end

% ----- Plots -----
figure; plot(x_all, N_all); title('Axial Load'); xlabel('x [m]'); ylabel('N [N]'); grid on;
figure; plot(x_all, T_all); title('Shear Load'); xlabel('x [m]'); ylabel('T [N]'); grid on;
figure; plot(x_all, M_all); title('Bending Moment'); xlabel('x [m]'); ylabel('M [Nm]'); grid on;
