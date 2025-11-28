function Xcp = computeFlareandFinsXcp( hf, d, S, Sb, Kf, l, dm)
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


    Xcp_over_d = ( (2/3) * (hf/d) * (S/Sb) + (Kf + 1 - S/Sb) * (l/d) - (hf/d) * ( (dm^2)/(d^2) - 1 ) * (S/Sb) ) / (1 + Kf);

    Xcp = Xcp_over_d * d;   % dimensional CP location

end
