function Xcp_total = computetotalFinXcp(opt, mission)
% Computes the total center of pressure of the vehicle including fins for both stages
%
% Inputs:
%   opt       : options structure (contains nStages)
%   mission   : mission structure (contains geometry, aerodynamics, etc.)
%
% Output:
%   Xcp_total : combined center of pressure [m]

    N = opt.nStages; % Current stage number (2=both stages, 1=stage 2 only, 0=no stages)

    % --- Body CP
    Xcp_body = computeXcpcomputeXcp(mission, opt);
    % --- Weighted average to get total Xcp
    S_body = mission.aerodynamics.bodyGeom.Aref;

    % --- Data fins stage 

    delta_le = mission.aerodynamics.finsGeomStage.delta_le;
    lambda_le = mission.aerodynamics.finsGeomStage.lambda_le;
    bfin   = mission.aerodynamics.finsGeomStage.bfin;
    be = 2*bfin;
    Sfin   = mission.aerodynamics.finsGeomStage.Sfin;
    Se = 2*Sfin;
    cmac = mission.aerodynamics.finsGeomStage.cmac;
    tmac = mission.aerodynamics.finsGeomStage.tmac;

    % --- Center of pressure fins stage


    Xcp_fin = computeFinXcp(mission, be, cmac, Se, maxQData);
        [CN_fin, ~, ~] = AerodynCoefFins_new(mission, cmac, lambda_le, Se, delta_le, tmac);


    % --- Fin CP and normal force coefficient
    
    % --- 2 stages 
    if N == 2
        
       Xcp_total = (Xcp_body*S_body + 2*Xcp_fin*Se*CN_fin) / ...
                    (S_body + 2*Se*CN_fin );

        
    end

    % --- 1 stages
    if N == 1

        Xcp_total = (Xcp_body*S_body + Xcp_fin*Se*CN_fin) / ...
                    (S_body + Se*CN_fin);

    % --- 0 stages
    else
        
        Xcp_total = Xcp_body;
    end
end
