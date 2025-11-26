function Xcp = computeXcp(N, alpha, lc1, lc2, lco)
% Calculates the center of pressure of the launcher
% Inputs:
%   alpha   : angle of attack [deg]
%   lc1, lc2   : length of stages 1 and 2, [m]
%   lco   : length of cone, [m]
% Output:
%   Xcp : center of pressure location (from top), [m]



    if N == 2
        l = lco + lc2 + lc1;
    elseif N == 1
        l = lco + lc2;
    else
        l = lco;
    end
    
    Xcp_over_l = 0.63*(1-(sin(alpha))^2)+ 0.5*l/lco*(sin(alpha))^2;
    
end
