function Xcp = computeXcpGlobal(mission, configuration,launcher, alpha, stageNumber)
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
N = launcher(1) - (stageNumber - 1);
lco = mission.capsule.height;
lc = zeros(1, N);
% ============================ SOLUTION ===================================

for i = (launcher(1)-stageNumber+1) : -1 : 1
    lc(i) = configuration.geometry.stage{i}.tanksLength;
    li(i) = configuration.geometry.stage{i}.interstage.length;
end

l = lco + sum(lc) ;

Xcp_over_l = 0.63*(1-(sin(alpha))^2)+ 0.5*l/lco*(sin(alpha))^2;
Xcp = Xcp_over_l*lco;


end