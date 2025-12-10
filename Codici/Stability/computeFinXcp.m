function Xcp = computeFinXcp(mission)
% computeFinXcp Computes the center of pressure (xCP) for a fin-stabilized launch vehicle
%
% Inputs:
%   vlauncher : launcher speed [m/s]
%   vsound    : speed of sound [m/s]
%   cmac      : fin mean chord [m]
%   be        : 2*fin axial base [m]
%   Se        : 2*fin surface [m^2]
%
% Output:
%   Xcp : center of pressure of fin [m]     



% ========================== DATA CONVERSION ==============================
be   = mission.aerodynamics.finsGeom.be;
cmac = mission.aerodynamics.finsGeom.cmac;
Se   = mission.aerodynamics.finsGeom.Se;
soundSpeed = mission.aerodynamics.soundspeedFun(mission.structure.hMaxQ);
launcherSpeed = mission.structure.vMaxQ' * mission.structure.vMaxQ;

% ========================= SOLUTION ======================================

Mach = launcherSpeed / soundSpeed;  

A = be^2 / Se;            

Mach0 = 0.7;    Xcp_sub = 0.25;
Mach1 = 2.0;    Xcp_sup = (A*sqrt(Mach1^2 - 1) - 0.67)/(2*A*sqrt(Mach1^2 - 1) - 1);
Machmid = 1.3;  Xmid = 0.44; % midpoint for smooth curve

Aeq = [Mach0^2, Mach0, 1;
       Mach1^2, Mach1, 1;
       Machmid^2, Machmid, 1];
Beq = [Xcp_sub; Xcp_sup; Xmid];

coeff = Aeq\Beq;

a = coeff(1); b = coeff(2); c = coeff(3);


if Mach < 0.7
    Xcp_on_cmac = Xcp_sub;
elseif Mach > 2
    Xcp_on_cmac = Xcp_sup;
else
    Xcp_on_cmac = a*Mach^2 + b*Mach + c;
end


Xcp = Xcp_on_cmac * cmac;

end