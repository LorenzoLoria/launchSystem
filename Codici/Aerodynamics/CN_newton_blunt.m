function CN = CN_newton_blunt(alpha, M, gamma, S_nose_wind, S_ref)
% CN_newton_blunt
%   Normal force coefficient of a blunt nose in hypersonic flow
%   using the modified Newtonian model.
%
%   CN = CN_newton_blunt(alpha, M, gamma, S_nose_wind, S_ref)
%
%   INPUTS:
%     alpha        : angle of attack [rad] (can be scalar or vector)
%     M            : freestream Mach number (hypersonic: M > 5)
%     gamma        : ratio of specific heats (typically 1.4)
%     S_nose_wind  : effective windward area of the blunt nose [m^2]
%     S_ref        : aerodynamic reference area [m^2]
%
%   OUTPUT:
%     CN           : normal force coefficient of the blunt nose
%
%   Model:
%     C_p(theta)   = C_p,max * sin^2(theta)
%     CN(alpha)    ~ C_p,max * sin^2(alpha) * (S_nose_wind / S_ref)

    % Stagnation pressure coefficient from modified Newtonian model
    Cp_max = cpmax_modified_newtonian(M, gamma);

    % Global normal force coefficient for blunt nose
    CN = Cp_max .* (sin(alpha).^2) .* (S_nose_wind ./ S_ref);
end
