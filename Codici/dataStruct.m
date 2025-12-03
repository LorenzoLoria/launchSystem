% data script
function [mission,optimisation] = dataStruct

mission = struct();

mission.launcher.engines{1}.name = 'Merlin1D' ;
mission.launcher.engines{2}.name = 'Raptor' ;
mission.launcher.engines{3}.name = 'RS-25' ;

mission.launcher.engines{1}.couple = 'RP1-LOX' ;
mission.launcher.engines{2}.couple = 'CH4-LOX' ;
mission.launcher.engines{3}.couple = 'LH2-LOX' ;

mission.launcher.engines{1}.ispZero = 283;
mission.launcher.engines{1}.ispVac = 311;

mission.launcher.engines{1}.thrust = 854*1e3;
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
mission.launcher.engines{2}.thrust = 2750*1e3;
mission.launcher.engines{2}.weight = 1525;
mission.launcher.engines{2}.OF = 3.6;
mission.launcher.engines{2}.oxDens = 1143;
mission.launcher.engines{2}.fuelDens = 423;
mission.launcher.engines{2}.effAreaZero = 1.3^2/4*pi;
mission.launcher.engines{2}.effAreaVac = 3.09^2/4*pi;

mission.launcher.engines{3}.ispZero = 366;
mission.launcher.engines{3}.ispVac = 452;
mission.launcher.engines{3}.thrust = 1860*1e3;
mission.launcher.engines{3}.weight = 3526;
mission.launcher.engines{3}.OF = 6;
mission.launcher.engines{3}.oxDens = 1143;
mission.launcher.engines{3}.fuelDens = 71;
mission.launcher.engines{3}.effAreaZero = 2.30^2/4*pi;
mission.launcher.engines{3}.effAreaVac = 2.30^2/4*pi;

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
optimisation.nStages = 2;
optimisation.stage{1}.engine = mission.launcher.engines{1};
optimisation.stage{2}.engine = mission.launcher.engines{3};

optimisation.stage{1}.nEngines = 4;
optimisation.stage{2}.nEngines = 1;

mProp1 = optimisation.stage{1}.nEngines*optimisation.stage{1}.engine.thrust/9.81*0.5;
epsS1 = 0.05;

mProp2 = optimisation.stage{2}.nEngines*optimisation.stage{2}.engine.thrust/9.81*0.3;
epsS2 = 0.06;

mS1 = epsS1/(1-epsS1)*mProp1;
mS2 = epsS2/(1-epsS2)*(mProp2);

mStage1 = mProp1+mS1;
mStage2 = mProp2+mS2;


optimisation.m0Tot = mStage1 + mStage2 + mission.capsule.weigth;

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


%optimisation.stage{2}.ispZero = mission.launcher.engines{2}.ispZero;
%optimisation.stage{2}.ispVac = mission.launcher.engines{2}.ispVac;

%optimisation.stage{2}.ispZero = mission.launcher.engines{2}.ispZero;
%optimisation.stage{2}.ispVac = mission.launcher.engines{2}.ispVac;


optimisation.stage{1}.length = 30;
optimisation.stage{2}.length = 15;
optimisation.stage{3}.length = 10;


optimisation.nDeval = 100;


% Aerodynamics

mission.aerodynamics.bodyGeom.l = 3;
mission.aerodynamics.bodyGeom.d = 0.3;
mission.aerodynamics.bodyGeom.Lnose = 0.6;
mission.aerodynamics.bodyGeom.Aref = pi*(mission.aerodynamics.bodyGeom.d/2)^2; 
mission.aerodynamics.bodyGeom.Anose = mission.aerodynamics.bodyGeom.Aref;
mission.aerodynamics.bodyGeom.Abase = mission.aerodynamics.bodyGeom.Aref;
mission.aerodynamics.bodyGeom.Aexit = 0;  
mission.aerodynamics.bodyGeom.phi = deg2rad(0);
mission.aerodynamics.bodyGeom.Ab = mission.aerodynamics.bodyGeom.Aref;
mission.aerodynamics.bodyGeom.Ap = mission.aerodynamics.bodyGeom.l * mission.aerodynamics.bodyGeom.d;



mission.aerodynamics.rootChord   = 1.81;                      % root chord [m]
mission.aerodynamics.tipChord   = 0.45;                       % tip chord [m] (triangolo puro)
mission.aerodynamics.semispan    = 1.81;                      % semispan [m]

mission.aerodynamics.finsGeom.Nfins = 1;
mission.aerodynamics.finsGeom.be = mission.aerodynamics.semispan;                  % span equivalente ~ semispan reale
mission.aerodynamics.finsGeom.Se = 0.5 * mission.aerodynamics.rootChord * mission.aerodynamics.semispan;       % area della fin [m^2]
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

mission.structure.nComponents = 8;
mission.structure.elementLength = [4,2,3,10,11,10,3,12]; % da modificare con funzione constance

mission.structure.SF = 2.5; % da verificare

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

mission.structure.SF = 1; % da verificare

%----------------------------------------------------------------------
% Al 7075-T6 - highly loaded fittings, secondary structures
%----------------------------------------------------------------------
mission.structure.rho      = 2720;
mission.structure.E        = 75e9;
mission.structure.yield    = 500e6;
mission.structure.ultimate = 560e6;

% Pressurizzazione serbatoi
mission.structure.tankPressure = [0, 0, 3.4e5, 0, 3.4e5];

% Dimensions 
mission.structure.diameter = 4;

% Interstages
mission.structures{1}.mInterstage = 480;
mission.structures{2}.mInterstage = 480;
mission.structures{3}.mInterstage = 480;

mission.structures{1}.lengthInterstage = 3.5;
mission.structures{2}.lengthInterstage = 3.5;
mission.structures{3}.lengthInterstage = 3.5;