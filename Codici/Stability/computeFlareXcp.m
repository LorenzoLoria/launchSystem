function Xcp = computeFlareXcp( lc1, lc2, lco, d, hf, db)
% Calculates the center of pressure of the launcher with a single flare at
% the end
% Inputs:
%   N    : Number of stages, [-]
%   lc1, lc2 : length of stages 1 and 2, [m]
%   lco  : length of cone, [m]
%   d    : reference diameter of launcher, [m]
%   hf   : flare length, [m]
%   db   : base radius of the flare cone, [m]
%
% Output:
%   Xcp  : center of pressure location (from top), [m]

    % Reference areas
    S  = pi*d^2/4;   % reference surface area [m^2]
    Sb = pi*db^2/4;  % base surface area for volume normalization [m^2]
    dm = (db + d)/2; % mean radius of the flare cone, [m]

    % Compute center of pressure
    
    l = lco + lc2 + lc1;
    Xcp_over_d = (2/3)*(hf/d)*(S/Sb) + (1 - S/Sb)*(l/d) - (hf/d)*((dm^2)/(d^2) - 1)*(S/Sb);

    
    Xcp = Xcp_over_d * d;

end

