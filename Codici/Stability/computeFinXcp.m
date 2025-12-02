function Xcp = computeFinXcp(vlauncher, vsound, cmac, be, Se)

% computeFinXcp   Computes the new center of pressure (xCP) for a 
%                 fin-stabilized launch vehicle.
%
% Inputs:
%   vlauncher   : speed of launcher [m/s]
%   vsound   : speed of sound [m/s]
%   cmac   : fin mean chord [m]
%   be   : 2*fin axial base [m]
%   Se   : 2*fin surface [m^2]
%
% Output:
%   xCP : center of pressure of fin [m]
    
    M = vlauncher/vsound;
    A = be^2/Se;


    if M < 0.7
        Xcp_on_cmac = 0.25;
 

    elseif M > 2
        Xcp_on_cmac = (A*(M^2-1)^1/2-0.67)/(2*A*(M^2-1)^1/2-1);

    else 
        Xcp_on_cmac = (A*(M^2-1)^1/2-0.67)/(2*A*(M^2-1)^1/2-1);
    end

    Xcp = Xcp_on_cmac*cmac;

     
end
