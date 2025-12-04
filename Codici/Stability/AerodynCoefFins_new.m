function [CN_surf, CD0_surf_friction, CD0_surf_wave] = AerodynCoefFins_new(alpha_p,A, vlauncher, vsound, Se, q, Sref, cmac, delta_le, lambda_le,b, tmac)
% Calculates aerodynamic coefficients for fins
% Inputs:
%   alpha_p    : local angle of attack [rad]
%   vlauncher  : speed of launcher at time t [m/s]
%   vsound     : speed of sound at time t [m/s]
%   b          : span of the fins [m]
%   q          : dynamic pressure at time t [Pa]
%   Se         : surface of the fin [m^2]
%   Sref       : reference surface [m^2]
%   cmac       : mean aerodynamic chord of fin [m]
%   delta_le   : leading edge sweep [deg]
%   lambda_le  : fin base angle [deg]
%   tmac       : max thickness of MAC [m]
%
% Outputs:
%   CN_surf           : normal force coefficient
%   CD0_surf_friction : surface friction drag coefficient
%   CD0_surf_wave     : wave drag coefficient

    M = vlauncher/vsound;
    M_ale = M * cosd(lambda_le);

    %Data proportional SaturnV 
    % A = 5.3;
    % Se = be^2*0.1875;
    % delta_le = 20.6;
    % lambda_le = 36.9;


    % --- Normal force coefficient
    if M > sqrt(1 + (8/(pi*A))^2)
        CN_surf = ((4*abs(sin(alpha_p)*cos(alpha_p)) / sqrt(M^2 - 1)) + 2*sin(alpha_p)^2) * Se / Sref;
    else
        CN_surf = ((pi*A/2*abs(sin(alpha_p)*cos(alpha_p)) + 2*sin(alpha_p)^2) * Se / Sref);
    end

    % --- CD0 surface friction
    CD0_surf_friction = 0.0133 * (M / (q*cmac))^0.2 * 2 * Se / Sref;

    % --- CD0 surface wave

    if M_ale < 1
        CD0_surf_wave = 0;
    else
        CD0_surf_wave = (1.429 / M_ale^2) * ((1.2*M_ale^2)^3.5 * (2.4/(2.8*M_ale^2 - 0.4))^2.5 - 1) * (sin(deg2rad(delta_le))^2 * cos(deg2rad(lambda_le)) * tmac * b) / Sref;
    end
    
end
