clear;
clc;

%% ============================= SCRIPT ====================================

% Cd sfera, cubo, Crew Dragon + limite Crew Dragon hypersonic
M_range = 0 : 0.01 : 10;

Cd_s = Cd_sphere(M_range);
Cd_c = Cd_cube(M_range);
Cd_CD = Cd_crewDragon(M_range);

Cd_s_asympt = Cd_s(end); % 0.92
Cd_c_asympt = Cd_c(end); % 1.67
% Crew dragon ~ 1.23

% Plot sphere and cube
figure;
plot(M_range, Cd_s, '-b', 'LineWidth', 1.5);
grid on;
hold on;
plot(M_range, Cd_c, '-r', 'LineWidth', 1.5);
plot(M_range, Cd_CD, '-g', 'LineWidth', 1.5);
yline(1.23, '--k', 'LineWidth', 1.5);
xlabel('Mach number [-]');
ylabel('Drag coefficient [-]');
title('Drag coefficient for a sphere and a cube as a function of M');
legend('Cd sphere', 'Cd cube', 'Cd Crew Dragon', 'Approx Cd_{hypersonic} Crew Dragon');




%% ============================= FUNCTIONS =================================
function Cd_s = Cd_sphere(M_range)

Cd_s = zeros(length(M_range), 1);
i = 0;
for M = M_range
    i = i + 1;
    if M <= 0.722
        Cd_s(i) = 0.45 * M^2 + 0.424;
    else
        Cd_s(i) = 2.1 * exp(-1.2 * (M + 0.35)) - 8.9 * exp(-2.2 * (M + 0.35)) + 0.92;
    end
end

end

% -------------------------------------------------------------------------
function Cd_c = Cd_cube(M_range)

Cd_c = zeros(length(M_range), 1);
i = 0;
for M = M_range
    i = i + 1;
    if M <= 1.15
        Cd_c(i) = 0.60 * M^2 + 1.04;
    else
        Cd_c(i) = 2.1 * exp(-1.16 * (M + 0.35)) - 6.5 * exp(-2.23 * (M + 0.35)) + 1.67;
    end
end

end

% -------------------------------------------------------------------------
function Cd_s = Cd_crewDragon(M_range)

Cd_s = zeros(length(M_range), 1);
i = 0;
for M = M_range
    i = i + 1;
    if M <= 0.722
        Cd_s(i) = 0.68 * M^2 + 0.424;
    else
        Cd_s(i) = 1.73 * exp(-1.02 * (M + 0.25)) - 8.9 * exp(-1.95 * (M + 0.35)) + 1.23;
    end
end

end





% REFERTENCES
% 1- https://ntrs.nasa.gov/api/citations/20110016614/downloads/20110016614.pdf (andamento Cd sfera e cubo)