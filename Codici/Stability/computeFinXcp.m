function Xcp = computeFinXcp(vlauncher, vsound, cmac, be, Se)
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

    M = vlauncher / vsound;  
    A = be^2 / Se;            


    Xcp_sub = 0.25;

    Xcp_sup = (A*sqrt(M^2-1) - 0.67) / (2*A*sqrt(M^2-1) - 1);

  
    coeff = Xcpfinscurve(A);  
    a = coeff(1); b = coeff(2); c = coeff(3);

  
    if M < 0.7
        Xcp_on_cmac = Xcp_sub;
    elseif M > 2
        Xcp_on_cmac = Xcp_sup;
    else
        Xcp_on_cmac = a*M^2 + b*M + c;
    end

    
    Xcp = Xcp_on_cmac * cmac;

end
