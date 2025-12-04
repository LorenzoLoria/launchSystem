function coeff = Xcpfinscurve(A)
% plotFinCPvsMach Plots fin Xcp vs Mach and returns transonic parabola coefficients
%
% Input:
%   A : fin slenderness factor (be^2 / Se)
%
% Output:
%   coeff : coefficients [a, b, c] of the transonic parabola (Xcp = a*M^2 + b*M + c)

% Define Mach numbers for plotting
M_sub  = linspace(0, 0.7, 50);
M_trans = linspace(0.7, 2, 100);
M_sup  = linspace(2, 5, 50);

% Subsonic Xcp
Xcp_sub = 0.25 * ones(size(M_sub));

% Supersonic Xcp
Xcp_sup = (A*sqrt(M_sup.^2 - 1) - 0.67) ./ (2*A*sqrt(M_sup.^2 - 1) - 1);

% Transonic parabola coefficients
M0 = 0.7;    X0 = 0.25;
M1 = 2.0;    X1 = (A*sqrt(M1^2 - 1) - 0.67)/(2*A*sqrt(M1^2 - 1) - 1);
Mmid = 1.3; Xmid = 0.44; % midpoint for smooth curve

Aeq = [M0^2, M0, 1;
       M1^2, M1, 1;
       Mmid^2, Mmid, 1];
Beq = [X0; X1; Xmid];

coeff = Aeq\Beq;  % [a; b; c]

% Compute transonic Xcp
Xcp_trans = coeff(1)*M_trans.^2 + coeff(2)*M_trans + coeff(3);

% Plot the full curve
% figure;
% plot(M_sub, Xcp_sub, 'b-', 'LineWidth', 2); hold on;
% plot(M_trans, Xcp_trans, 'r-', 'LineWidth', 2);
% plot(M_sup, Xcp_sup, 'g-', 'LineWidth', 2);
% xlabel('Mach number');
% ylabel('X_{cp} / cmac');
% title('Fin Center of Pressure vs Mach');
% grid on;
% legend('Subsonic (M<0.7)', 'Transonic (0.7 ≤ M ≤ 2)', 'Supersonic (M>2)', 'Location', 'best');

end
