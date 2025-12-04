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
launcherSpeed = norm(mission.structure.vMaxQ);

Mach = launcherSpeed / soundSpeed;  

A = be^2 / Se;            


Xcp_sub = 0.25;

Xcp_sup = (A*sqrt(Mach^2-1) - 0.67) / (2*A*sqrt(Mach^2-1) - 1);


coeff = Xcpfinscurve(A);  
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