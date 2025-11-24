function Xcp = computeXcp(alpha, Lb, h)
% Calculates the center of pressure of the launcher
% Inputs:
%   alpha   : angle of attack [deg]
%   Lb   : length of vehicle, [m]
%   h   : length of cone, [m]
% Output:
%   Xcp : center of pressure location (from top), [m]


    % --- Compute nondimensional center of pressure Xcp/h
    Xcp_over_h = 0.63*(1-(sin(alpha))^2)+ 0.5*Lb/h*(sin(alpha))^2;

    % --- Convert to dimensional Xcp
    Xcp = Xcp_over_h * h;
end
