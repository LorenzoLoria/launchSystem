% data script
function [mission,optimisation,settings] = dataStruct

mission = struct();

mission.launcher.engines{1}.name = 'Merlin1D' ;
mission.launcher.engines{2}.name = 'Raptor' ;
mission.launcher.engines{3}.name = 'Vinci' ;
mission.launcher.engines{4}.name = 'RS-25' ;
mission.launcher.engines{5}.name   = 'HM7B';
mission.launcher.engines{6}.name   = 'Aestus2'; % !! TRL = 6 !!
mission.launcher.engines{7}.name   = 'AJ10_118K';
mission.launcher.engines{8}.name   = 'RL10C-1';
mission.launcher.engines{9}.name   = 'BE-3U';

mission.launcher.engines{1}.couple = 'RP1-LOX' ;
mission.launcher.engines{2}.couple = 'CH4-LOX' ;
mission.launcher.engines{3}.couple = 'LH2-LOX' ;
mission.launcher.engines{4}.couple = 'LH2-LOX' ;
mission.launcher.engines{5}.couple = 'LH2-LOX';
mission.launcher.engines{6}.couple = 'MMH-N2O4';
mission.launcher.engines{7}.couple = 'A50-N2O4'; 
mission.launcher.engines{8}.couple = 'LH2-LOX';
mission.launcher.engines{9}.couple = 'LH2-LOX';

mission.launcher.engines{1}.ispZero = 283;
mission.launcher.engines{1}.ispVac = 311;

mission.launcher.engines{1}.thrustZero = 854*1e3;
mission.launcher.engines{1}.thrustVacum = 981*1e3;
mission.launcher.engines{1}.weight = 470;
mission.launcher.engines{1}.OF = 2.36;
mission.launcher.engines{1}.oxDens = 1143;
mission.launcher.engines{1}.fuelDens = 835;   
mission.launcher.engines{1}.effAreaZero = 0.9;
mission.launcher.engines{1}.effAreaVac = 4.9;

% Questi valori dovrebbero essere dati da GA:
mission.launcher.engines{1}.mPropellant1 =  200e3;
mission.launcher.engines{1}.ms1 = 12e3;
mission.launcher.engines{1}.mPropellant2 = 37.6e3;
mission.launcher.engines{1}.ms2 = 3e3;
mission.launcher.engines{1}.m0 =208e3;% mission.launcher.engines{1}.mPropellant1 + mission.launcher.engines{1}.ms1 + mission.launcher.engines{1}.mPropellant2 + mission.launcher.engines{1}.ms2 + 8600;

mission.launcher.engines{2}.ispZero = 350;
mission.launcher.engines{2}.ispVac = 380;
mission.launcher.engines{2}.thrustZero = 2750*1e3;
mission.launcher.engines{2}.thrustVacum = 2986 *1e3;
mission.launcher.engines{2}.weight = 1525;
mission.launcher.engines{2}.OF = 3.6;
mission.launcher.engines{2}.oxDens = 1143;
mission.launcher.engines{2}.fuelDens = 423;
mission.launcher.engines{2}.effAreaZero = 1.3^2/4*pi;
mission.launcher.engines{2}.effAreaVac = 3.09^2/4*pi;

mission.launcher.engines{3}.ispZero = 270;
mission.launcher.engines{3}.ispVac = 457.2;
mission.launcher.engines{3}.thrustZero = 270/458.2 * 180 *1e3;
mission.launcher.engines{3}.thrustVacum = 180 *1e3;
mission.launcher.engines{3}.weight = 550;
mission.launcher.engines{3}.OF = 6.1;
mission.launcher.engines{3}.oxDens = 1143;
mission.launcher.engines{3}.fuelDens = 71;
mission.launcher.engines{3}.effAreaZero = 1.84^2/4*pi;
mission.launcher.engines{3}.effAreaVac = 1.84^2/4*pi;

mission.launcher.engines{4}.ispZero = 366;
mission.launcher.engines{4}.ispVac = 452;
mission.launcher.engines{4}.thrustZero = 1860*1e3;
mission.launcher.engines{4}.thrustVacum = 2279*1e3;
mission.launcher.engines{4}.weight = 3526;
mission.launcher.engines{4}.OF = 6;
mission.launcher.engines{4}.oxDens = 1143;
mission.launcher.engines{4}.fuelDens = 71;
mission.launcher.engines{4}.effAreaZero = 2.3^2/4*pi;
mission.launcher.engines{4}.effAreaVac = 2.3^2/4*pi;

% =========================
% 8) HM7B (Ariane upper stage)
% =========================


% Motore solo-vuoto: ispZero/thrustZero sono valori indicativi
mission.launcher.engines{5}.ispZero = 360;        % [s] stima
mission.launcher.engines{5}.ispVac  = 445;        % [s] ~444–446

mission.launcher.engines{5}.thrustZero  = 60e3;   % [N] circa
mission.launcher.engines{5}.thrustVacum = 62.7e3; % [N] dato Ariane 5

mission.launcher.engines{5}.weight = 155;         % [kg] ordine di grandezza 

mission.launcher.engines{5}.OF      = 5.0;        % rapporto O/F ~5:1
mission.launcher.engines{5}.oxDens  = 1143;       % [kg/m^3] LOX (come gli altri)
mission.launcher.engines{5}.fuelDens= 71;         % [kg/m^3] LH2

% Nozzle exit diameter ≈ 2.1 m → area efficace
d_HM7B = 2.1;                                     % [m]
mission.launcher.engines{5}.effAreaZero = (d_HM7B^2)*pi/4;
mission.launcher.engines{5}.effAreaVac  = (d_HM7B^2)*pi/4;

% =========================
% 9) Aestus II / RS-72 (upper stage ipergolico, ~55 kN)
% =========================


mission.launcher.engines{6}.ispZero = 320;        % [s] ipotesi (motore solo-vuoto)
mission.launcher.engines{6}.ispVac  = 338;        % [s] ~336–340

mission.launcher.engines{6}.thrustZero  = 55.4e3; % [N]
mission.launcher.engines{6}.thrustVacum = 55.4e3; % [N]

mission.launcher.engines{6}.weight = 138;         % [kg] dry mass

mission.launcher.engines{6}.OF      = 1.9;        % mixture ratio (massa) tipico
% Densità tipiche ipergolici (ordine di grandezza)
mission.launcher.engines{6}.oxDens  = 1440;       % [kg/m^3] N2O4
mission.launcher.engines{6}.fuelDens= 880;        % [kg/m^3] MMH

% Nozzle exit diameter ≈ 1.36 m (Aestus II)
d_Aestus2 = 1.36;                                 % [m]
mission.launcher.engines{6}.effAreaZero = (d_Aestus2^2)*pi/4;
mission.launcher.engines{6}.effAreaVac  = (d_Aestus2^2)*pi/4;

% =========================
% 10) AJ10-118K (Delta-K, ~44 kN)
% =========================

mission.launcher.engines{7}.ispZero = 300;       % [s] stima (altitude engine)
mission.launcher.engines{7}.ispVac  = 320.5;     % [s] ~319–321

mission.launcher.engines{7}.thrustZero  = 43.6e3; % [N] ~9.8 klbf
mission.launcher.engines{7}.thrustVacum = 43.6e3; % [N]

mission.launcher.engines{7}.weight = 125;        % [kg] ≈ 275 lb

mission.launcher.engines{7}.OF      = 1.9;       % mixture ratio tipico AJ10
% Densità da letteratura (ordine di grandezza)
mission.launcher.engines{7}.oxDens  = 1448;      % [kg/m^3] N2O4
mission.launcher.engines{7}.fuelDens= 903;       % [kg/m^3] Aerozine-50

% Diametro ugello ≈ 1.7 m
d_AJ10 = 1.7;                                     % [m]
mission.launcher.engines{7}.effAreaZero = (d_AJ10^2)*pi/4;
mission.launcher.engines{7}.effAreaVac  = (d_AJ10^2)*pi/4;

% =========================
% 7) RL10C-1 (Centaur / Vulcan)
% =========================

mission.launcher.engines{8}.ispZero = 380;        % ipotesi (motore solo-vuoto)
mission.launcher.engines{8}.ispVac  = 450;        % ~449.7 s arrotondato 
mission.launcher.engines{8}.thrustZero  = 1.01e5; % N, ≈ vac
mission.launcher.engines{8}.thrustVacum = 1.01e5; % N, ~101 kN in vuoto
mission.launcher.engines{8}.weight      = 190;    % kg, ~420 lb 

mission.launcher.engines{8}.OF      = 5.5;   % nominal mixture ratio RL10C-1
mission.launcher.engines{8}.oxDens  = 1143;
mission.launcher.engines{8}.fuelDens= 71;

% RL10C-1 nozzle diameter ~57" ≈ 1.45 m
mission.launcher.engines{8}.effAreaZero = (1.45^2)*pi/4;
mission.launcher.engines{8}.effAreaVac  = (1.45^2)*pi/4;

% =========================
% 6) BE-3U (New Glenn upper)
% =========================

mission.launcher.engines{9}.ispZero = 390;        % ipotesi (non opera davvero a SL)
mission.launcher.engines{9}.ispVac  = 445;        % s, dichiarato da Blue Origin 
mission.launcher.engines{9}.thrustZero  = 7.78e5; % N, assumo ≈ thrust vac
mission.launcher.engines{9}.thrustVacum = 7.78e5; % N, ~778 kN in vuoto
mission.launcher.engines{9}.weight      = 1000;   % kg, stima ordine grandezza (<~1.1 t) 

mission.launcher.engines{9}.OF      = 5.5;   % ≈ rapporto tipico BE-3 LH2/LOX
mission.launcher.engines{9}.oxDens  = 1143;
mission.launcher.engines{9}.fuelDens= 71;

% diametro ugello non pubblicato in modo chiaro: assumo simile a RS-25
mission.launcher.engines{9}.effAreaZero = (2.3^2)*pi/4; % m^2 (ipotesi)
mission.launcher.engines{9}.effAreaVac  = (2.3^2)*pi/4;

mission.capsule.weigth = 8600;
mission.capsule.Area = 4^2*pi/4;
mission.capsule.supersonicCD = 1.23;
mission.capsule.subsonicCD = 0.45;
mission.capsule.height = 2.9;

mission.environment.altRange = (-1000:100:1000000);
rhoVal = zeros(length(mission.environment.altRange),1);
Tval= zeros(length(mission.environment.altRange),1);
%warning('off', 'all'); 

%for i=1:length(mission.environment.altRange)
%    [T,rho] = atmosnrlmsise00(mission.environment.altRange(i), 45, 120, 2020,50,1300);
%    rhoVal(i) = rho(6);
%    Tval(i) = T(1,2) ; 
%end

%warning('on', 'all');   
%mission.environment.rho = rhoVal;
%mission.environment.T = Tval;

%mission.environment.gridInterp = griddedInterpolant( ...
%    mission.environment.altRange, ...   % grid points
%    mission.environment.rho, ...        % values
%    'linear', ...                       % interpolation method
%    'linear');  

%mission.environment.gridInterpTemp = griddedInterpolant( ...
%    mission.environment.altRange, ...   % grid points
%    mission.environment.T, ...        % values
%    'linear', ...                       % interpolation method
%    'linear');  


mission.environment.rEarth = 6371e3;
mission.environment.g0 = 9.81;
mission.environment.GM = 398600.433e9;


mission.options.fmincon = optimoptions("fmincon","Display","iter","MaxIterations",200,'MaxFunctionEvaluations',10000,'StepTolerance',1e-19,'OptimalityTolerance',1e-6,'FunctionTolerance',1e-15,'ConstraintTolerance',1e-10,'EnableFeasibilityMode',true,"UseParallel",false);
mission.options.gaOptions = optimoptions('ga', 'PlotFcn',{'gaplotbestf', 'gaplotbestindiv'}, 'display', 'iter','MaxStallGenerations', 10, ...
        'FunctionTolerance', 1e-6, 'EliteCount',  6,...
        'MaxGenerations', 100, 'PopulationSize', 40, ...
        'NonlinearConstraintAlgorithm', 'penalty',"UseParallel",false);


% Robe di prova con il booster

mission.launcher.booster.mass = 5600; 
mission.launcher.booster.cp = [0 0 -2]';

mission.environment.rhoFun = @(h) 1.29*exp(-h/8433);


mission.launchBase.latInitial = deg2rad(-39.261515237910196) ;
mission.launchBase.lonInitial = deg2rad(+177.86521705520965) ;
mission.launchBase.initialPointECI = 6371000*[cos(mission.launchBase.latInitial)*cos(mission.launchBase.lonInitial ); cos(mission.launchBase.latInitial)*sin(mission.launchBase.lonInitial ); sin(mission.launchBase.latInitial)]' ;
mission.target.latInitial = - mission.launchBase.latInitial;
mission.target.lonInitial = mission.launchBase.lonInitial - pi ;
mission.target.initialPointECI = 6371000*[cos(mission.target.latInitial)*cos(mission.target.lonInitial); cos(mission.target.latInitial)*sin(mission.target.lonInitial); sin(mission.target.latInitial) ];
mission.target.omega = 2*pi/86164.1 ;


n = cross(mission.launchBase.initialPointECI,[0,1,0])/(norm(cross(mission.launchBase.initialPointECI,[0,1,0])));


ex = mission.launchBase.initialPointECI' / norm(mission.launchBase.initialPointECI);
ez = n' / norm(n) ;
ey = cross(ez,ex)/norm(cross(ez,ex));

rot = [ex,ey,ez]';

%EarthPlot(mission.environment.rEarth)
%hold on
%plot3(mission.target(1),mission.target(2),mission.target(3),'ro')
%plot3(mission.initialPoint(1),mission.initialPoint(2),mission.initialPoint(3),'bo')
%hold off


    mission.Rfinal = rot ;

    mission.optimisation.GA.variables = 5;    


%% differentStages Simulation 
optimisation = struct();
optimisation.nStages = 1;
optimisation.stage{1}.engine = mission.launcher.engines{2};
optimisation.stage{2}.engine = mission.launcher.engines{3};
optimisation.stage{3}.engine = mission.launcher.engines{3};

optimisation.stage{1}.nEngines = 2;
optimisation.stage{2}.nEngines = 1.25;
optimisation.stage{3}.nEngines = 0.6;

optimisation.stage{1}.percentage = 0.75;
optimisation.stage{2}.percentage= 0.5;
optimisation.stage{3}.percentage= 0.4;

mProp1 = optimisation.stage{1}.nEngines*optimisation.stage{1}.engine.thrustZero/9.81*optimisation.stage{1}.percentage;
epsS1 = 0.05;

mProp2 = optimisation.stage{2}.nEngines*optimisation.stage{2}.engine.thrustVacum/9.81*optimisation.stage{2}.percentage;
epsS2 = 0.06;

mProp3 = optimisation.stage{3}.nEngines*optimisation.stage{3}.engine.thrustVacum/9.81*optimisation.stage{3}.percentage;
epsS3 = 0.06;

mS1 = epsS1/(1-epsS1)*mProp1;
mS2 = epsS2/(1-epsS2)*(mProp2);
mS3 = epsS3/(1-epsS3)*(mProp3);

mStage1 = mProp1+mS1;
mStage2 = mProp2+mS2;
mStage3 = mProp3+mS3;



optimisation.m0Tot = mStage1 + mission.capsule.weigth;
optimisation.totalMass = mStage1 + mission.capsule.weigth;

optimisation.stage{1}.mStage = mStage1;
optimisation.stage{1}.mProp = mProp1;
optimisation.stage{1}.structuralMass = mS1;
%optimisation.stage{1}.ispZero = mission.launcher.engines{1}.ispZero;
%optimisation.stage{1}.ispVac = mission.launcher.engines{1}.ispVac;

%optimisation.stage{1}.ispZero = mission.launcher.engines{1}.ispZero;
%optimisation.stage{1}.ispVac = mission.launcher.engines{1}.ispVac;


optimisation.stage{2}.mStage = mStage2;
optimisation.stage{2}.mProp = mProp2;
optimisation.stage{2}.structuralMass = mS2;

optimisation.stage{3}.mStage = mStage3;
optimisation.stage{3}.mProp = mProp3;
optimisation.stage{3}.structuralMass = mS3;
%optimisation.stage{2}.ispZero = mission.launcher.engines{2}.ispZero;
%optimisation.stage{2}.ispVac = mission.launcher.engines{2}.ispVac;

%optimisation.stage{2}.ispZero = mission.launcher.engines{2}.ispZero;
%optimisation.stage{2}.ispVac = mission.launcher.engines{2}.ispVac;


optimisation.stage{1}.length = 30;
optimisation.stage{2}.length = 15;
optimisation.stage{3}.length = 10;


optimisation.nDeval = 100;


% ============================ Aerodynamics ===============================

mission.aerodynamics.bodyGeom.l = 54.9;
mission.aerodynamics.bodyGeom.d = 4;
mission.aerodynamics.bodyGeom.Lnose = 2.9;
mission.aerodynamics.bodyGeom.Aref = pi*(mission.aerodynamics.bodyGeom.d/2)^2; 
mission.aerodynamics.bodyGeom.Anose = mission.aerodynamics.bodyGeom.Aref;
mission.aerodynamics.bodyGeom.Abase = mission.aerodynamics.bodyGeom.Aref;
mission.aerodynamics.bodyGeom.Aexit = 0;  
mission.aerodynamics.bodyGeom.phi = deg2rad(0);
mission.aerodynamics.bodyGeom.Ab = mission.aerodynamics.bodyGeom.Aref;
mission.aerodynamics.bodyGeom.Ap = mission.aerodynamics.bodyGeom.l * mission.aerodynamics.bodyGeom.d; % sezione verticale del cilindro



mission.aerodynamics.rootChord   = 1.81;                      % root chord [m]
mission.aerodynamics.tipChord   = 0.45;                       % tip chord [m] (triangolo puro)
mission.aerodynamics.semispan    = 1.81;                      % semispan [m]

mission.aerodynamics.finsGeom.Nfins = 1;
mission.aerodynamics.finsGeom.be = mission.aerodynamics.semispan;                  % span equivalente ~ semispan reale
mission.aerodynamics.finsGeom.Se = 2 * 0.5 * (mission.aerodynamics.rootChord + mission.aerodynamics.tipChord) * mission.aerodynamics.semispan;       % area della fin [m^2]
mission.aerodynamics.finsGeom.cmac = (2/3) * mission.aerodynamics.rootChord;       % mean aerodynamic chord [m]
mission.aerodynamics.finsGeom.delta_le = 45;
mission.aerodynamics.finsGeom.lambda_le = 0;
mission.aerodynamics.finsGeom.b = 2 * mission.aerodynamics.rootChord;
mission.aerodynamics.finsGeom.tmac = 0.08 * mission.aerodynamics.finsGeom.cmac;    



mission.aerodynamics.bodyInfo.isPowered = false;
mission.aerodynamics.bodyInfo.a_sub = 0;
mission.aerodynamics.bodyInfo.b_sub = 1;
mission.aerodynamics.bodyInfo.Cdn = 1.2;



mission.aerodynamics.finsInfo.rho = 1.225;
mission.aerodynamics.finsInfo.vsound = 340;

alt = linspace(0,500000,100);
soundSpeedVec = load("soundspeed.mat");
soundSpeedVec = linspace(soundSpeedVec.soundspeed(1),soundSpeedVec.soundspeed(end),100);
%mission.aerodynamics.soundspeedFun = @(h) interp1(alt,soundSpeedVec,h,"linear","extrap");


mission.aerodynamics.soundspeedFun= griddedInterpolant( ...
   alt, ...   % grid points
   soundSpeedVec, ...        % values
   'linear', ...                       % interpolation method
   'linear');  


[~,~,pressure] = atmosisa(alt);
pressure = pressure.*(pressure>0.5);

mission.environment.pressure= griddedInterpolant( ...
   alt, ...   % grid points
   pressure, ...        % values
   'linear', ...                       % interpolation method
   'linear');  

% ========================= Structures ====================================

mission.structure.SF = 1.5; % da verificare

% Al2219 - cryogenic tanks and primary structures for 1st/2nd stage
mission.structure.rho      = 2840;
mission.structure.E        = 72e9;
mission.structure.yield    = 390e6;
mission.structure.ultimate = 480e6;

% Pressurizzazione serbatoi
mission.structure.tankPressure       = [0, 0, 0, 3.4e5, 3.4e5, 3.4e5, 3.4e5, 0];

% Dimensions 
mission.structure.diameter           = ones(8,1) * 4;


% ========================= Structures ====================================

%----------------------------------------------------------------------
% Al 7075-T6 - highly loaded fittings, secondary structures
%----------------------------------------------------------------------

mission.structure.rho = [2840, 2720, 2780, 2810];
mission.structure.E   = [72e9, 75e9, 73e9, 72e9];
mission.structure.ultimate = [480e6, 560e6, 483e6, 572e6];
mission.structure.yield = [390e6, 500e6, 345e6, 503e6];

% Pressurizzazione serbatoi
mission.structure.tankPressure = 3.4e5;

% Dimensions 
mission.structure.diameter = 4;

% Interstages
mission.structures{1}.mInterstage = 150;
mission.structures{2}.mInterstage = 150;
mission.structures{3}.mInterstage = 150;

mission.structures{1}.lengthInterstage = 3.5;
mission.structures{2}.lengthInterstage = 3.5;
mission.structures{3}.lengthInterstage = 3.5;









%% opt vere questa volta

settings = struct();

settings.lowerBoundsFMC(:,:,1) = [0.1*ones(mission.optimisation.GA.variables,1);...
                        0*ones(mission.optimisation.GA.variables,1)];
settings.lowerBoundsFMC(:,:,2) = [0.1*ones(mission.optimisation.GA.variables,1);...
                        0*ones(mission.optimisation.GA.variables,1)];
settings.lowerBoundsFMC(:,:,3) = [0.1*ones(mission.optimisation.GA.variables,1);...
                        0*ones(mission.optimisation.GA.variables,1)];


settings.upperBoundsFMC(:,:,1) = [ones(mission.optimisation.GA.variables,1);...
                        90*ones(mission.optimisation.GA.variables,1)];
settings.upperBoundsFMC(:,:,2) = [ones(mission.optimisation.GA.variables,1);...
                        120*ones(mission.optimisation.GA.variables,1)];
settings.upperBoundsFMC(:,:,3) = [ones(mission.optimisation.GA.variables,1);...
                        150*ones(mission.optimisation.GA.variables,1)];


settings.lowerBoundsGA = settings.lowerBoundsFMC(:);
settings.upperBoundsGA = settings.upperBoundsFMC(:);


settings.gaTrajOptions = optimoptions("ga", ...
                        "Display","iter", ...
                        "MaxGenerations",20, ...
                        "PopulationSize",50,...
                        "UseParallel",false,...
                        "FunctionTolerance", 1e-4);


settings.fminconTrajOptions = optimoptions("fmincon",...
                                "Display","iter",...
                                "MaxIterations",200,...
                                'MaxFunctionEvaluations',10000,...
                                'StepTolerance',1e-19,...
                                'OptimalityTolerance',1e-6,...
                                'FunctionTolerance',1e-15,...
                                'ConstraintTolerance',1e-10, ...
                                'EnableFeasibilityMode',true,...
                                "UseParallel",false);

settings.TrajOptimisationPoints = 5;

settings.fsolveOptions = optimoptions('fsolve','Display','none');





end




