function[CL,CD,CN,CA, mainbodyCL, mainbodyCD, finsCL, finsCD] = CLCDcomputation(Mach,alpha,dynamicPressure,isPoweredFlag,mission,currentStage,dimensions,engineVec, finsVec)



%Last stage, coupler to the capsule is not considered here, but rather in
%the nosecone
% lCylinder = geoStages{nStages}.tanksLength;                                            % lunghezza totale [m] lunghezza cilindro mediato
% dCylinder = lCylinder *stageRadius*2;     % diametro [m] diametro cilindro mediato
% 
% %Other stages
% for i=nStages-1:-1:currentStage
%     lCylinder = lCylinder + geoStages{i}.tanksLength + geoStages{i}.interstage.length ;
%     dCylinder = dCylinder + geoStages{i}.tanksLength * geoStages{i}.radius*2 +...
%                 geoStages {i}.interstage.length *...
%                 (geoStages {i}.radius + geoStages {i+1}.radius);
% 
% end

lCylinder = dimensions(2);
dCylinder = dimensions(3);

dCylinder = dCylinder / lCylinder ;
lambda = lCylinder / dCylinder;                                                                 % Fineness ratio

lNose = dimensions(4);                % lunghezza nose [m] lunghezza capsula + interstage
Aref = pi/4 * dCylinder^2 ;                                                                     % area di riferimento [m^2] cross section cilindro mediato
Sref = Aref ;
Anose = dimensions(5);              % area nose [m^2] max tra cross section capsula e ultimo stadio
phi = pi/4 ;                                                                                    % angolo di giunzione [rad] angolo tra ultimo stadio e interstage di connessione alla capsula
nEngines = engineVec(1);

if currentStage == 1
    AexitTot = engineVec(2) * nEngines;   % aree uscita motori
else
    AexitTot = engineVec(3) * nEngines;
end

stage1Radius = dimensions(6);
boatTailRadius = dimensions(7) ;      % boat tail not present r=radius of the current stagev
Abase = pi * boatTailRadius^2 ;       % area base [m^2] se è presente una boattail non coincide con l'area dello stadio corrente
Ab = pi/4 * dCylinder^2;              % area max cross section
Ap = dCylinder * lCylinder ;          % area razzo vista da lato

% Fins Data
cr = finsVec(1);                      % root chord [m]
ct = finsVec(2);                       % tip chord [m] (triangolo puro)
bfin = finsVec(3);                      % semispan [m]
Nfins = finsVec(4);              
Sfin = finsVec(5);     % area della fin [m^2]
cmac = finsVec(6);     % mean aerodynamic chord [m]
deltaLE = finsVec(7);  % [deg]
lambdaLE = finsVec(8); % [deg]
tmac = finsVec(9);    
aSub = 0;
bSub = 1;


%--------------------------------------------------------------------------
%-------------------------CA COMPUTATION-----------------------------------
%--------------------------------------------------------------------------

if dynamicPressure < 1e-16
    
    CL = 0; CD = 0; CN = 0; CA = 0; mainbodyCL = 0; mainbodyCD = 0; finsCL = 0; finsCD = 0;

else

    % -----------------------------
    % 1) Wave drag (C_A)_W
    % -----------------------------
    
    % 1.1) Contributi supersonici (M >= 1.3)
    CAW_sharp = (1.586 + 1.834 ./ (Mach^2)) * ...
                (atan(0.5 / (lNose / dCylinder)))^1.69;
    
    CAW_hemi  = 0.665 * (1.586 + 1.834 / (Mach^2));
    
    CAW_sup   = CAW_sharp .* ((Aref - Anose) / Aref) + ...
                CAW_hemi  .* (Anose / Aref);
    
    % 1.2) Limite fortemente subsonico (M ≈ 0)
    CAW_M0 = 0.8 * sin(phi).^2;   % scalare
    
    % 1.3) Modello subsonico fino a M <= 0.8
    CAW_sub = aSub * (Mach^bSub) + CAW_M0;
    
    % regioni:
    %   M <= 0       : poniamo CAW = CAW_M0
    %   0 < M <= 0.8 : CAW_sub
    %   0.8 < M < 1.3: interpolazione lineare tra CAW_sub(0.8) e CAW_sup(1.3)
    %   M >= 1.3     : CAW_sup
    
    if Mach <= 0.1
        CAW = CAW_M0;
        
    elseif Mach > 0.1 && Mach<=0.8    % 0 < M <= 0.8
        CAW = CAW_sub;
    elseif Mach > 0.8 && Mach<=1.3    % 0.8 < M < 1.3  -> interpolazione lineare
        % valore a M=0.8 (subsonico)
        M1   = 0.8;
        CA1  = aSub * (M1^bSub) + CAW_M0;
        % valore a M=1.3 (supersonico)
        M2   = 1.3;
        CAW_sharp_13 = (1.586+1.834/M2^2) * ...
                       (atan(0.5/(lNose / dCylinder)))^1.69;
        CAW_hemi_13  = 0.665*(1.586+1.834/M2^2);
        CA2  = CAW_sharp_13*((Aref-Anose)/Aref) + ...
               CAW_hemi_13*(Anose/Aref);
        % interp lineare in funzione di M
        CAW = CA1 + (CA2 - CA1)*(Mach - M1)/(M2 - M1);
    elseif Mach > 1.3                 % M >= 1.3
        CAW = CAW_sup;
    end
    
    
    % -----------------------------
    % 2) Base drag (C_A)_B
    % -----------------------------
    
    % Coefficiente base drag riferito all'area di base
    
    
    CD0B = 0.12 + 0.13 .* (Mach.^2)*(Mach<1) + (0.25/Mach)*(Mach>=1);
    
    
    % Area di base efficace (motore acceso)
    Abase_eff = Abase - AexitTot;
    
    
    CA_B = (CD0B .* (Abase_eff ./ Aref))*isPoweredFlag + CD0B .* (Abase ./ Aref)*(not(isPoweredFlag)) ;
    
    
    % -----------------------------
    % 3) Friction drag (C_A)_f  (Jerger/Fleeman in SI)
    % -----------------------------
    % (C_A)_f ≈ 0.091 * (L/d) * ( M / (q * L) )^0.2
    % q = pressione dinamica [Pa], L in [m]
    
    CA_f = 0 + (0.091*lambda*(Mach/(dynamicPressure*lCylinder))^0.2)*(Mach~=0);
    
    
    
    % -----------------------------
    % 4) C_A,alpha=0 e dipendenza in alpha
    % -----------------------------
    
    CA0 = CAW + CA_B + CA_f;
    
    % C_A,body(M, alpha) = C_A,alpha=0(M) * cos^2(alpha)
    CAbody = CA0*cos(alpha)^2;
    
    
    %--------------------------------------------------------------------------
    %-------------------------CN COMPUTATION-----------------------------------
    %--------------------------------------------------------------------------
    
    Ab_over_A = Ab / Aref;
    Ap_over_A = Ap / Aref;
    Cdn = 1.2; %mission.aerodynamics.bodyInfo.Cdn;  %lascia come è costante
    
    eta = (0.05 * lambda + 0.52).*(Mach<1) + 1.*(Mach>1);
    
    
    % -----------------------------
    % Calcolo dei termini slender-body e crossflow
    % -----------------------------
    CN_slender   = Ab_over_A .* sin(2 .* alpha) .* cos(alpha ./ 2);
    CN_crossflow = eta .* Cdn .* Ap_over_A .* (sin(alpha).^2);
    
    % Somma totale
    CNbody = CN_slender + CN_crossflow;
    
    % se M > 5 aggiungo contributo Newtoniano del blunt nose
    
    if Mach > 5
        % area windward efficace del nose (~ Aref o Anose)
        S_nose_wind = Aref;
        CN_nose_newt = CN_newton_blunt(alpha, Mach, 1.4, S_nose_wind, Aref);
        CNbody  = CNbody + CN_nose_newt;
    end
    
    
    %--------------------------------------------------------------------------
    %------------------------------FINS----------------------------------------
    %--------------------------------------------------------------------------
    
    A = (bfin*2)^2/(Sfin*2);
    M_ale = Mach * cosd(lambdaLE);
    alphafins = 10 * pi / 180;
    
    % --- Normal force coefficient
    mRef = sqrt(1 + (8/(pi*A))^2);
    CN_surf =( ((4*abs(sin(alphafins)*cos(alphafins)) / sqrt(Mach^2 - 1)) + 2*sin(alphafins)^2) * Sfin / Sref)*(Mach > mRef) + (((pi*A/2*abs(sin(alphafins)*cos(alphafins)) + 2*sin(alphafins)^2) * Sfin / Sref))*(Mach <= mRef);
    
    
    % --- CD0 surface friction
    
    CD0_surf_friction = 0.0133 * (Mach / (dynamicPressure*cmac))^0.2 * 2 * Sfin / Sref;
    
    
    
    % --- CD0 surface wave
    
    CD0_surf_wave = 0 + ((1.429 / M_ale^2) * ((1.2*M_ale^2)^3.5 * (2.4/(2.8*M_ale^2 - 0.4))^2.5 - 1) * (sin((deltaLE))^2 * cos((lambdaLE)) * tmac * bfin) / Sref).*(M_ale>=1);
    
    
    CN_fins_tot = Nfins * CN_surf;
    CD0_fins_tot = Nfins * (CD0_surf_friction + CD0_surf_wave);
    CA_fins_tot = CD0_fins_tot * cos(alphafins)^2;
    
    % Somma totali
    CN = CNbody + CN_fins_tot;
    CA = CAbody + CA_fins_tot;
    
    % Calcola il coefficiente di lift (C_L)
    CL = CN * cos(alpha) - CA * sin(alpha);
    
    % Calcola il coefficiente di drag (C_D)
    CD = CA * cos(alpha) + CN * sin(alpha);
    
    % CL, CD Main Body
    mainbodyCL = CNbody * cos(alpha) - CAbody * sin(alpha);
    mainbodyCD = CAbody * cos(alpha) + CNbody * sin(alpha);
    
    % CL, CD Fins
    finsCL = CN_fins_tot * cos(alphafins) - CA_fins_tot * sin(alphafins);
    finsCD = CA_fins_tot * cos(alphafins) + CN_fins_tot * sin(alphafins);

end


end