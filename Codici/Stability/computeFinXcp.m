function Xcp = computeFinXcp(vlauncher, vsound, cmac, be, Se)
% computeFinXcp Computes the center of pressure (xCP) for a fin-stabilized launch vehicle.
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

    M = vlauncher / vsound;  % Mach number
    A = be^2 / Se;            % fin slenderness factor

    % Subsonic (M < 0.7)
    Xcp_sub = 0.25;

    % Supersonic (M > 2)
    Xcp_sup = (A*sqrt(M^2-1) - 0.67) / (2*A*sqrt(M^2-1) - 1);

    % Parabolic fit for transonic region (0.7 <= M <= 2)
    if M < 0.7
        Xcp_on_cmac = Xcp_sub;
    elseif M > 2
        Xcp_on_cmac = Xcp_sup;
    else
        % Fit a parabola: Xcp = a*M^2 + b*M + c
        % Pass through the endpoints: (0.7, Xcp_sub) and (2, Xcp_sup)
        % Also ensure derivative continuity (optional)
        
        % Simplest: quadratic through endpoints and midpoint at M=1.35
        M0 = 0.7;   X0 = Xcp_sub;
        M1 = 2.0;   X1 = Xcp_sup;
        Mmid = 1.35; Xmid = (X0 + X1)/2;  % midpoint value for smooth curve

        % Solve for coefficients a, b, c
        Aeq = [M0^2, M0, 1;
               M1^2, M1, 1;
               Mmid^2, Mmid, 1];
        Beq = [X0; X1; Xmid];

        coeff = Aeq\Beq;  % [a; b; c]

        Xcp_on_cmac = coeff(1)*M^2 + coeff(2)*M + coeff(3);
    end

    % Final CP position
    Xcp = Xcp_on_cmac * cmac;

end
