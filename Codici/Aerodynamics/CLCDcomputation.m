function[CL,CD,CN,CA, mainbodyCL, mainbodyCD, finsCL, finsCD] = CLCDcomputation(Mach,alpha,dynamicPressure,isPoweredFlag,mission,currentStage,opt)

nStages = length(opt.stage);
geoStages = opt.geometry.stage;
optStage = opt.stage;
q = dynamicPressure;                       % pressione dinamica [Pa]


%Last stage, coupler to the capsule is not considered here, but rather in
%the nosecone
stageRadius = geoStages{nStages}.radius;
lCylinder = geoStages{nStages}.tanksLength;                                            % lunghezza totale [m] lunghezza cilindro mediato
dCylinder = lCylinder *stageRadius*2;     % diametro [m] diametro cilindro mediato

%Other stages
for i=nStages-1:-1:currentStage
    lCylinder = lCylinder + geoStages{i}.tanksLength + geoStages{i}.interstage.length ;
    dCylinder = dCylinder + geoStages{i}.tanksLength * geoStages{i}.radius*2 +...
                geoStages {i}.interstage.length *...
                (geoStages {i}.radius + geoStages {i+1}.radius);
    
end

dCylinder = dCylinder / lCylinder ;
lambda = lCylinder / dCylinder;                                                                 % Fineness ratio

lNose = mission.capsule.height + geoStages{nStages}.interstage.length ;                % lunghezza nose [m] lunghezza capsula + interstage
Aref = pi/4 * dCylinder^2 ;                                                                     % area di riferimento [m^2] cross section cilindro mediato
Sref = Aref ;
Anose = max(mission.capsule.Area , pi*stageRadius^2) ;                   % area nose [m^2] max tra cross section capsula e ultimo stadio
phi = pi/4 ;                                                                                    % angolo di giunzione [rad] angolo tra ultimo stadio e interstage di connessione alla capsula

nEngines = optStage{currentStage}.nEngines;

if currentStage == 1
    AexitTot = optStage{currentStage}.engine.effAreaZero * nEngines;   % aree uscita motori
else
    AexitTot = optStage{currentStage}.engine.effAreaVac * nEngines;
end

stage1Radius = geoStages{1}.radius;
boatTailRadius = geoStages{currentStage}.radius ;                                      %boat tail not present r=radius of the current stagev
Abase = pi * boatTailRadius^2 ;                                                                 % area base [m^2] se è presente una boattail non coincide con l'area dello stadio corrente
Ab = pi/4 * dCylinder^2;                                                                        %area max cross section
Ap = dCylinder * lCylinder ;                                                                    %area razzo vista da lato
ce = 4.525 * stage1Radius*2 / 10;                                               %mean chord fins??? mediato rispetto a saturn5?
be = 4.525 * stage1Radius*2 / 10;                                               %semi span delle fin mediato rispetto a saturn5?
Se = 0.5 * be^2;                                                                                %superficie di una fin
cmac = 2/3 * ce;                                                                                %mean aerochord fin
deltaLE = pi/4;                                                                                   %angolo rombo ???
lambdaLE = 0;                                                                                   % sweep leading edge
b = 2 * be ;                                                                                    % 2*be
tmac = 0.08 * cmac;                                                                             %spessore massimo sezione
Nfins = 4;                                                                                      %numero fins
aSub = 0;                                                                                       %parametri wave drag
bSub = 1;                                                                                       %parametri wave drag





%--------------------------------------------------------------------------
%-------------------------CA COMPUTATION-----------------------------------
%--------------------------------------------------------------------------

% -----------------------------
% 1) Wave drag (C_A)_W
% -----------------------------

% 1.1) Contributi supersonici (M >= 1.3)
CAW_sharp = (1.586 + 1.834 ./ (Mach.^2)) .* ...
            (atan(0.5 ./ (lNose ./ dCylinder))).^1.69;

CAW_hemi  = 0.665 * (1.586 + 1.834 ./ (Mach.^2));

CAW_sup   = CAW_sharp .* ((Aref - Anose) ./ Aref) + ...
            CAW_hemi  .* (Anose ./ Aref);

% 1.2) Limite fortemente subsonico (M ≈ 0)
CAW_M0 = 0.8 * sin(phi).^2;   % scalare

% 1.3) Modello subsonico fino a M <= 0.8
CAW_sub = aSub .* (Mach.^bSub) + CAW_M0;

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
    CAW_sharp_13 = (1.586 + 1.834 / M2^2) * ...
                   (atan(0.5 / (lNose / dCylinder)))^1.69;
    CAW_hemi_13  = 0.665 * (1.586 + 1.834 / M2^2);
    CA2  = CAW_sharp_13 * ((Aref - Anose) / Aref) + ...
           CAW_hemi_13  * (Anose / Aref);
    % interp lineare in funzione di M
    CAW = CA1 + (CA2 - CA1) .* (Mach - M1) / (M2 - M1);
elseif Mach > 1.3                 % M >= 1.3
    CAW = CAW_sup;
end


% -----------------------------
% 2) Base drag (C_A)_B
% -----------------------------

% Coefficiente base drag riferito all'area di base

if Mach < 1
    CD0B = 0.12 + 0.13 .* (Mach.^2);
else
    CD0B = 0.25 ./ Mach;
end

% Area di base efficace (motore acceso)
Abase_eff = Abase - AexitTot;

if isPoweredFlag
    CA_B = CD0B .* (Abase_eff ./ Aref);
else
    CA_B = CD0B .* (Abase ./ Aref);
end

% -----------------------------
% 3) Friction drag (C_A)_f  (Jerger/Fleeman in SI)
% -----------------------------
% (C_A)_f ≈ 0.091 * (L/d) * ( M / (q * L) )^0.2
% q = pressione dinamica [Pa], L in [m]
if Mach ==0
    CA_f = 0;
else
    CA_f = 0.091 * lambda .* ( Mach ./ (q * lCylinder) ).^0.2;
end
% -----------------------------
% 4) C_A,alpha=0 e dipendenza in alpha
% -----------------------------

CA0 = CAW + CA_B + CA_f;

% C_A,body(M, alpha) = C_A,alpha=0(M) * cos^2(alpha)
CAbody = CA0 .* cos(alpha).^2;


%--------------------------------------------------------------------------
%-------------------------CN COMPUTATION-----------------------------------
%--------------------------------------------------------------------------

Ab_over_A = Ab / Aref;
Ap_over_A = Ap / Aref;
Cdn = mission.aerodynamics.bodyInfo.Cdn;  %lascia come è costante

if Mach < 1
    eta = 0.05 * lambda + 0.52;
else
    eta = 1;
end

% -----------------------------
% Calcolo dei termini slender-body e crossflow
% -----------------------------
CN_slender   = Ab_over_A .* sin(2 .* alpha) .* cos(alpha ./ 2);
CN_crossflow = eta .* Cdn .* Ap_over_A .* (sin(alpha).^2);

% Somma totale
CNbody = CN_slender + CN_crossflow;

%--------------------------------------------------------------------------
%------------------------------FINS----------------------------------------
%--------------------------------------------------------------------------

A = be^2/Se;
M_ale = Mach * cosd(lambdaLE);

% --- Normal force coefficient
if Mach > sqrt(1 + (8/(pi*A))^2)
    CN_surf = ((4*abs(sin(alpha)*cos(alpha)) / sqrt(Mach^2 - 1)) + 2*sin(alpha)^2) * Se / Sref;
else
    CN_surf = ((pi*A/2*abs(sin(alpha)*cos(alpha)) + 2*sin(alpha)^2) * Se / Sref);
end

% --- CD0 surface friction
CD0_surf_friction = 0.0133 * (Mach / (q*cmac))^0.2 * 2 * Se / Sref;

% --- CD0 surface wave

if M_ale < 1
    CD0_surf_wave = 0;
else
    CD0_surf_wave = (1.429 / M_ale^2) * ((1.2*M_ale^2)^3.5 * (2.4/(2.8*M_ale^2 - 0.4))^2.5 - 1) * (sin((deltaLE))^2 * cos((lambdaLE)) * tmac * b) / Sref;
end

CN_fins_tot = Nfins * CN_surf;
CD0_fins_tot = Nfins * (CD0_surf_friction + CD0_surf_wave);
CA_fins_tot = CD0_fins_tot * cos(alpha)^2;

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
finsCL = CN_fins_tot * cos(alpha) - CA_fins_tot * sin(alpha);
finsCD = CA_fins_tot * cos(alpha) + CN_fins_tot * sin(alpha);
end
