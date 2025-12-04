function [CN_flare] = AerodynCoefFlare(S1, S2, d)
% Calculates aerodynamic coefficients for fins
% Inputs:
%   S1         : upstream section [m]
%   S2         : downstream section of the flare [m]
%   d          : diameter of the launcher [m]
%
% Outputs:
%   CN_flare          : normal force coefficient of the flare


    Sref = pi*d^2/4;
    CN_flare = 8*(S2-S1)/Sref;
 
end
