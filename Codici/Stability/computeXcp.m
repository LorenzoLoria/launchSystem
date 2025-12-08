function Xcp = computeXcp(mission, opt)
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
N = opt.nStages;
alpha = mission.structure.alphaQmax;
lco = mission.capsule.height;
lc1 = opt.stage{1}.length;
lc2 = opt.stage{2}.length;
lc3 = opt.stage{3}.length;
li1 = mission.structures{1}.lengthInterstage;
li2 = mission.structures{2}.lengthInterstage;
li3 = mission.structures{3}.lengthInterstage;


% ============================ SOLUTION ===================================
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
Xcp = Xcp_over_l*lco;

end