%function Xcp = computeFinXcp(vlauncher, vsound, cmac, be, Se)
function Xcp = computeFinXcp(N, lc2, lc1, lco, alpha, d, hf, db)
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
    
    % M = vlauncher/vsound;
    % A = be^2/Se;
    % 
    % 
    % if M < 0.7
    %     Xcp_on_cmac = 0.25;
    % end
    % 
    % if M > 2
    %     Xcp_on_cmac = (A*(M^2-1)^1/2-0.67)/(2*A*(M^2-1)^1/2-1);
    % end
    % 
    % Xcp = Xcp_on_cmac*cmac;

    %Sforza combo fins and flare 

    S  = pi*d^2/4;   % reference surface area [m^2]
    Sb = pi*db^2/4;  % base surface area for volume normalization [m^2]
    dm = (db + d)/2; 
   
    Kf = 1/2;

    if N == 2
        l = lco + lc1 + lc2;
        xCP_over_d = ((2/3)*(lc2/d)*(S/Sb)+ (Kf + 1 - S/Sb)*(l/d) - (hf/d)*( (dm^2)/(d^2) - 1 )*(S/Sb)) / (1 + Kf);
        Xcp = xCP_over_d * d;
    elseif N == 1 
        Xcp_over_l = 0.63*(1-(sin(alpha))^2)+ 0.5*(lco+lc2)/lco*(sin(alpha))^2;
        Xcp = Xcp_over_l*(lco+lc2);
    else 
        Xcp = lco * 2/3; 
    end


   
   
end
