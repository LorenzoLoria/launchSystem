function validate_Fleeman_xCP
% VALIDATE_FLEEMAN_XCP
% Funzione per la validazione della formula di Fleeman per la posizione del centro di pressione
% con valori di L_body/L_nose = 1, 2, 5, 10 e angolo di attacco alpha da 0° a 90°.

clc;
clear;
close all;

%% 1) Parametri di riferimento
L_nose = 1; % Lunghezza del naso (unità arbitrarie)
L_body_ratios = [1, 2, 5, 10]; % Rapporto L_body / L_nose

alpha_deg = linspace(0, 90, 91);  % Angolo di attacco da 0° a 90°
alpha_rad = deg2rad(alpha_deg);   % Conversione in radianti

%% 2) Calcolo della posizione del centro di pressione (x_CP) usando la formula di Fleeman
x_CP = zeros(length(alpha_deg), length(L_body_ratios));

for i = 1:length(L_body_ratios)
    L_body = L_nose * L_body_ratios(i);  % Calcolare L_body
    x_CP(:, i) = 0.63 * (1 - sin(alpha_rad).^2) + 0.5 * (L_body / L_nose) .* sin(alpha_rad).^2;
end

%% 3) Grafico della posizione del centro di pressione (x_CP)
figure;
hold on;
grid on;
box on;

colors = lines(length(L_body_ratios)); % Color palette for different L_body/L_nose ratios

% Plotting delle curve
for i = 1:length(L_body_ratios)
    plot(alpha_deg, x_CP(:, i), 'LineWidth', 2, 'Color', colors(i, :), ...
        'DisplayName', sprintf('L_{body}/L_{nose} = %.1f', L_body_ratios(i)));
end

% Etichette e titolo
xlabel('\alpha [deg]', 'FontSize', 15);
ylabel('$\frac{x_{CP,\,body}}{\ell_{nose}}$', 'Interpreter', 'latex', 'FontSize', 30);
ylim([0 6])
%title('Validazione della formula di Fleeman per la posizione del centro di pressione');
lgd = legend('Location','best');  
lgd.FontSize = 15;   


end
