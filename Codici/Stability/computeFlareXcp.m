function Xcp = computeFlareXcp(N, lc1, lc2, lco,d, hf, db)
% Calculates the center of pressure of the launcher
% Inputs:
%   d   : reference diameter of launcher, [m]
%   h   : cone length, [m]
%   hf  : flare length, [m]
%   db  : base radius of the flare cone, [m]
%
% Output:
%   Xcp : center of pressure location (from top), [m]

    S = pi*d^2/4; %reference surface area [m^2]
    Sb = pi*db^2/4; %base surface area for volume normalization [m^2]
    dm = (db+d)/2; % mean radius of the flare cone, [m]
  
    % --- Compute nondimensional slender-body volume term v/(Sb*d)
    Xcp_over_d = (2/3)*(h/d)*(S/Sb) + (1 - S/Sb)*(Lb/d) - (hf/d)*((dm^2)/(d^2) - 1)*(S/Sb);

    % --- Convert to dimensional Xcp
    Xcp = Xcp_over_d * d;
end
