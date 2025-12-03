function Xcp = computeXcp(mission, opt)
% Calculates the center of pressure of the launcher
% Inputs:
%   N : number of stage (variable, depends on flight condition)
%   alpha   : angle of attack [deg]
%   lc1, lc2   : length of stages 1 and 2, [m]
%   lco   : length of cone, [m]
% Output:
%   Xcp : center of pressure location (from top), [m]


% ========================== DATA CONVERSION ==============================
N = opt.stage

if N == 3
    l = lco + li3 + lc3 + li2 + lc2  + li1 + lc1; 
elseif N == 2
    l = lco + li2 + lc2 + li1 + lc1;
elseif N == 1
    l = lco + + li2 + lc2;
else
    l = lco;
end

Xcp_over_l = 0.63*(1-(sin(alpha))^2)+ 0.5*l/lco*(sin(alpha))^2;
Xcp = Xcp_over_l*l;

end
