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

    % --- Data fins 1st stage 

    delta_le = mission.aerodynamics.finsGeomStage.delta_le;
    lambda_le = mission.aerodynamics.finsGeomStage.lambda_le;
    bfin1   = mission.aerodynamics.finsGeomStage.bfin;
    be1 = 2*bfin1;
    Sfin1   = mission.aerodynamics.finsGeomStage.Sfin;
    Se1 = 2*Sfin1;
    cmac1 = mission.aerodynamics.finsGeomStage1.cmac;
    tmac1 = mission.aerodynamics.finsGeomStage1.tmac;

    % --- Center of pressure fins 1st stage


    Xcp_fin1 = computeFinXcp(mission, be1, cmac1, Se1, maxQData);
        [CN_fin1, ~, ~] = AerodynCoefFins_new(mission, cmac1, lambda_le, Se1, delta_le, tmac1);

    % --- Data fins 2nd stage 

    bfin2   = d2/d1*bfin1;
    be2 = 2*bfin2;
    Sfin2   = (d2/d1)^2*Sfin1;
    Se2 = 2*Sfin2;
    cmac2 = d2/d1*cmac1; 
    tmac2 = d2/d1*tmac1;

    % --- Center of pressure fins 2nd stage

    Xcp_fin2 = computeFinXcp(mission, be2, cmac2, Se2, maxQData);
        [CN_fin2, ~, ~] = AerodynCoefFins_new(mission, cmac2, lambda_le, Se2, delta_le, tmac2);


    
    % --- Fin CP and normal force coefficient
    
    % --- 2 stages 
    if N == 2
        
       Xcp_total = (Xcp_body*S_body + Xcp_fin1*Se1*CN_fin1 + Xcp_fin2*Se2*CN_fin2) / ...
                    (S_body + Se1*CN_fin1 + Se2*CN_fin2);

        
    end

    % --- 1 stages
    if N == 1

        Xcp_total = (Xcp_body*S_body + Xcp_fin1*Se1*CN_fin1) / ...
                    (S_body + Se1*CN_fin1);

    % --- 0 stages
    else
        
        Xcp_total = Xcp_body;
    end
end
