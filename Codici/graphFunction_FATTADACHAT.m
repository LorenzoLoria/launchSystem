% Plot h_apogee vs eps_p including MER_solid feasibility
clear; close all; clc;

% --- rocket data (come nel tuo esempio) ---
rocket.Isp  = 280;      % s
rocket.m0   = 12300;    % kg
rocket.T    = 450000;   % N
rocket.mpay = 800;      % kg
rocket.v0   = 0;        % m/s
rocket.h0   = 0;        % m

g = 9.81;

% --- evaluation grid for eps_p ---
eps_min = 0.6;   % estendo un po' il range per vedere il comportamento
eps_max = 0.95;
N = 1000;
eps_vec = linspace(eps_min, eps_max, N)';

% Preallocazioni
h_apogee = NaN(size(eps_vec));
feasible = false(size(eps_vec));
error_m_s = NaN(size(eps_vec));
mp_vec = NaN(size(eps_vec));

for k = 1:length(eps_vec)
    eps_p = eps_vec(k);
    mp = eps_p * rocket.m0;
    mp_vec(k) = mp;
    
    % protezioni numeriche minime
    if mp <= 0
        continue
    end
    
    tb = mp * rocket.Isp * g / rocket.T;
    c = rocket.Isp * g;
    
    MR = (rocket.m0 - mp) / rocket.m0;   % m_f / m0
    if MR <= 0
        continue
    end
    
    % velocità e altezza al cutoff (formule analitiche)
    v_ct = rocket.v0 - c * log(MR) - g * tb;
    h_ct = rocket.h0 + c * rocket.m0 / mp * tb * (MR * log(MR) - MR + 1) - 0.5 * g * tb^2;
    
    % delta h in coast = v^2/(2*g)
    h_ap = h_ct + 0.5 * v_ct.^2 / g;
    h_apogee(k) = h_ap;
    
    % massa secca / struttura e vincolo MER_solid
    mi = rocket.m0 - mp;            % massa immediatamente dopo burnout
    m_s = mi - rocket.mpay;         % massa strutturale stimata dal modello semplice
    ms_emp = MER_solid(mp);         % massa inerte dalla correlazione empirica
    if ms_emp <= 0
        % se la correlazione produce valore non sensato, consideriamo non ammissibile
        error_m_s(k) = NaN;
        feasible(k) = false;
    else
        error_m_s(k) = (m_s - ms_emp) / ms_emp;
        feasible(k) = (h_ap >= 0) && (abs(error_m_s(k)) <= 0.05);
    end
end

% --- Plot ---
figure('Color','w','Position',[200 200 900 500]);
hold on; grid on;
% linea dell'apogeo (anche dove NaN)
plot(eps_vec, h_apogee/1000, 'LineWidth', 1.5);   % convert to km on y axis

% punti Feasible / Infeasible
idx_feas = find(feasible & ~isnan(h_apogee));
idx_infeas = find(~feasible & ~isnan(h_apogee));

% sovrapponi punti: non forzare colori particolari, ma rendere leggibile
plot(eps_vec(idx_feas), h_apogee(idx_feas)/1000, 'o', 'MarkerSize',5, 'MarkerFaceColor', [0 .6 0], 'MarkerEdgeColor', [0 .3 0]);
plot(eps_vec(idx_infeas), h_apogee(idx_infeas)/1000, 'o', 'MarkerSize',5, 'MarkerFaceColor', [0.8 0 0], 'MarkerEdgeColor', [0.5 0 0]);

% bounds originali (se vuoi evidenziarli)
xline(0.7, '--', 'Lower bound = 0.7', 'LabelHorizontalAlignment','left');
xline(0.9, '--', 'Upper bound = 0.9', 'LabelHorizontalAlignment','right');

xlabel('\epsilon_p = m_p / m_0');
ylabel('h_{apogee} [km]');
title('h_{apogee} vs \epsilon_p  (punti verdi: ammissibili; rossi: non ammissibili)');
legend({'h_{apogee} (line)', 'ammissibile', 'non ammissibile'}, 'Location','best');
xlim([eps_min eps_max]);

% --- Optional: mostra la curva di errore relativo sulla secondaria asse Y ---
ax1 = gca;
ax2 = axes('Position', ax1.Position, 'Color','none', 'YAxisLocation','right', 'XTick',[], 'Box','off');
hold(ax2, 'on');
plot(ax2, eps_vec, error_m_s, 'LineStyle','-.', 'LineWidth', 1);
ylabel(ax2, 'error\_m\_s (relative)');
ylim(ax2, [-0.5 0.5]);    % scala utile per vedere +-50%
set(ax2, 'YColor', [0.15 0.15 0.15]);

% --- stampa alcuni risultati numerici chiave ---
[~, idx_max] = max(h_apogee);
fprintf('Max h_apogee (ignoring feasibility) = %.2f km at eps_p = %.4f\n', h_apogee(idx_max)/1000, eps_vec(idx_max));
feas_eps = eps_vec(idx_feas);
if ~isempty(feas_eps)
    [~, idx_best_feas] = max(h_apogee(idx_feas));
    fprintf('Best feasible h_apogee = %.2f km at eps_p = %.4f\n', h_apogee(idx_feas(idx_best_feas))/1000, eps_vec(idx_feas(idx_best_feas)));
else
    fprintf('Nessuna soluzione ammissibile trovata nel range considerato.\n');
end

%% --- sottofunzione MER_solid ---
function ms = MER_solid(mp)
    pc = 7; % assumption (come nel tuo codice)
    params_Me = [800.340063456389, 0.0970643753730160, ...
                 -162.431518885189, -0.0180992880437423, -8.48759341267327, -6441.47167561289];
    Minert = @(params_Me, mp, pc) params_Me(1) + params_Me(2).*mp + params_Me(3).*pc ...
        + params_Me(4).*(pc + params_Me(5)).*(mp + params_Me(6));
    if mp > 5000
        ms = Minert(params_Me, mp, pc);
    else
        % per mp <= 5000 usi la legge logaritmica
        ms = 140.94.*log(mp) - 823.29;
    end
end
