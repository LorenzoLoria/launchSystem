function xCP = computeFinXcp(l, d, h, hf, dm, S, Sb, Kf)
% computeFinXcp   Computes the new center of pressure (xCP) for a 
%                 fin-stabilized launch vehicle.
%
% Inputs:
%   l   : total length of vehicle [m]
%   d   : reference body diameter [m]
%   h   : cone length [m]
%   hf  : fin axial extent (flare length) [m]
%   dm  : outer diameter at fin base (flare diameter) [m]
%   S   : reference area (usually pi*d^2/4) [m^2]
%   Sb  : base area for normalization [m^2]
%   Kf  : fin correction factor (default 1/2)
%
% Output:
%   xCP : center of pressure from nose tip [m]

    % Default Kf if not provided (typical missile assumption)
    if nargin < 8 || isempty(Kf)
        Kf = 1/2;
    end

    % --- Non-dimensional xCP/d expression ---
    xCP_over_d = ((2/3)*(h/d)*(S/Sb)+ (Kf + 1 - S/Sb)*(l/d) - (hf/d)*( (dm^2)/(d^2) - 1 )*(S/Sb)) / (1 + Kf);

    % Convert nondimensional → dimensional
    xCP = xCP_over_d * d;
end
