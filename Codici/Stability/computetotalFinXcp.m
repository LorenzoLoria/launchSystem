function Xcp_total = computetotalFinXcp(opt, mission, vlauncher, vsound, be, Se, Sref)
% computeTotalFinXcp Computes the total center of pressure of the vehicle including fins
%
% Inputs:
%   N         : current stage number (2=stage1, 1=stage2, 0=no stage)
%   alpha     : angle of attack [rad]
%   lc1, lc2  : stage 1 and 2 lengths [m]
%   lco       : cone length [m]
%   vlauncher : launcher speed [m/s]
%   vsound    : speed of sound [m/s]
%   be        : fin axial base [m]
%   Se        : fin surface area [m^2]
%   Sref      : vehicle reference area [m^2]
%
% Output:
%   Xcp_total : combined center of pressure [m]

    N     = opt.nStages;
    cmac = mission.aerodynamics.finsGeom.cmac ;

    % --- Body CP
    Xcp_body = computeXcpcomputeXcp(mission, opt);

    % --- If stage 1, include fins
    if N == 2
        % Fin CP
        Xcp_fin = computeFinXcp(vlauncher, vsound, cmac, be, Se);

        % Estimate fin normal force coefficient
        alpha_p = mission.alpha;    % small angle, in radians
        b       = be/2;             % fin span
      % b = mission.aerodynamics.finsGeom.b ;
        delta_le = 20;              % example leading edge sweep [deg]
      % delta_le = mission.aerodynamics.finsGeom.delta_le ;
        lambda_le = 36.9;           % example base angle [deg]
      % lambda_le = mission.aerodynamics.finsGeom.lambda_le
      
        tmac = 0.1*cmac;            % max thickness
      % tmac = mission.aerodynamics.finsGeom.tmac ;
      
        q = 0.5*1.225*vlauncher^2;  % dynamic pressure, rho=1.225 kg/m3

        [CN_fin, ~, ~] = AerodynCoefFins_new(alpha_p, be^2/Se, vlauncher, vsound, Se, q, Sref, cmac, delta_le, lambda_le, b, tmac);

        % Body area
        S_body = Sref - Se;

        % Weighted average to get total Xcp
        Xcp_total = (Xcp_body*S_body + Xcp_fin*Se*CN_fin) / (S_body + Se*CN_fin);
    else
        % No fins: total CP = body CP
        Xcp_total = Xcp_body;
    end
end
