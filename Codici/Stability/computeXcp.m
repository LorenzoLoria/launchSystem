function Xcp = computeXcp(mission, opt,launcher)
% Calculates the center of pressure of the launcher
% Inputs:
%   N : number of stage (variable, depends on flight condition)
%   alpha   : angle of attack [deg]
%   lc1, lc2, lc3   : length of stages 1, 2 and 3, [m]
%   li1, li2, li3   : length of interstages 1, 2 and 3 [m]
%   lco   : length of cone, [m]
% Output:
%   Xcp : center of pressure location (from top), [m]


% ========================== DATA CONVERSION ==============================
N = launcher(1);
alpha = 0;
lco = mission.capsule.height;

% ============================ SOLUTION ===================================
if N == 3
    lc1 = opt.geometry.stage{1}.length;
    lc2 = opt.geometry.stage{2}.length;
    lc3 = opt.geometry.stage{3}.length;
    li1 = opt.geometry.stage{1}.interstage.length;
    li2 = opt.geometry.stage{2}.interstage.length;
    li3 = opt.geometry.stage{3}.interstage.length;
    l   = opt.geometry.totalLength; 
elseif N == 2
    lc1 = opt.geometry.stage{1}.length;
    lc2 = opt.geometry.stage{2}.length;
    li1 = opt.geometry.stage{1}.interstage.length;
    li2 = opt.geometry.stage{2}.interstage.length;
    l   = opt.geometry.totalLength;
elseif N == 1
    lc1 = opt.geometry.stage{1}.length;
    li1 = opt.geometry.stage{1}.interstage.length;
    l   = opt.geometry.totalLength;
else
    l = lco;
end

Xcp_over_l = 0.63*(1-(sin(alpha))^2)+ 0.5*l/lco*(sin(alpha))^2;
Xcp = Xcp_over_l*lco;

end