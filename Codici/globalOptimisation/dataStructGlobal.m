function [mission,settings] = dataStructGlobal


%% Mission constants

mission = struct();

% ============================ Engines ===============================
mission.engines{1}.name = 'Merlin1D' ;
mission.engines{2}.name = 'Raptor' ;
mission.engines{3}.name = 'Vinci' ;
mission.engines{4}.name = 'RS-25' ;

mission.engines{1}.couple = 'RP1-LOX' ;
mission.engines{2}.couple = 'CH4-LOX' ;
mission.engines{3}.couple = 'LH2-LOX' ;
mission.engines{4}.couple = 'LH2-LOX' ;


mission.engines{1}.ispZero = 283;
mission.engines{1}.ispVac = 311;
mission.engines{1}.thrustZero = 854*1e3;
mission.engines{1}.thrustVacum = 981*1e3;
mission.engines{1}.weight = 470;
mission.engines{1}.OF = 2.36;
mission.engines{1}.oxDens = 1143;
mission.engines{1}.fuelDens = 835;   
mission.engines{1}.effAreaZero = 0.9;
mission.engines{1}.effAreaVac = 4.9;
mission.engines{1}.chamberRadius = 0.609/2;
mission.engines{1}.length = 2.69;


mission.engines{2}.ispZero = 350;
mission.engines{2}.ispVac = 380;
mission.engines{2}.thrustZero = 2750*1e3;
mission.engines{2}.thrustVacum = 2986 *1e3;
mission.engines{2}.weight = 1525;
mission.engines{2}.OF = 3.6;
mission.engines{2}.oxDens = 1143;
mission.engines{2}.fuelDens = 423;
mission.engines{2}.effAreaZero = 1.3^2/4*pi;
mission.engines{2}.effAreaVac = 3.09^2/4*pi;
mission.engines{2}.chamberRadius = 1.3/2;
mission.engines{2}.length = 3.1;


mission.engines{3}.ispZero = 270;
mission.engines{3}.ispVac = 457.2;
mission.engines{3}.thrustZero = 270/458.2 * 180 *1e3;
mission.engines{3}.thrustVacum = 180 *1e3;
mission.engines{3}.weight = 550;
mission.engines{3}.OF = 6.1;
mission.engines{3}.oxDens = 1143;
mission.engines{3}.fuelDens = 71;
mission.engines{3}.effAreaZero = 1.84^2/4*pi;
mission.engines{3}.effAreaVac = 1.84^2/4*pi;
mission.engines{3}.chamberRadius = 0.7/2;
mission.engines{3}.length = 3.22;


mission.engines{4}.ispZero = 366;
mission.engines{4}.ispVac = 452;
mission.engines{4}.thrustZero = 1860*1e3;
mission.engines{4}.thrustVacum = 2279*1e3;
mission.engines{4}.weight = 3526;
mission.engines{4}.OF = 6;
mission.engines{4}.oxDens = 1143;
mission.engines{4}.fuelDens = 71;
mission.engines{4}.effAreaZero = 2.3^2/4*pi;
mission.engines{4}.effAreaVac = 2.3^2/4*pi;
mission.engines{4}.chamberRadius = 0.45/2;
mission.engines{4}.length = 4.3;
% ============================ Capsule ===============================

mission.capsule.weight = 8600;
mission.capsule.Area = 4^2*pi/4;
mission.capsule.radius = 2;
mission.capsule.supersonicCD = 1.23;
mission.capsule.subsonicCD = 0.45;
mission.capsule.height = 2.9;


% ============================ Launch base and target ===============================

mission.launchBase.latInitial = deg2rad(-39.261515237910196) ;
mission.launchBase.lonInitial = deg2rad(+177.86521705520965) ;
mission.launchBase.initialPointECI = 6371000*[cos(mission.launchBase.latInitial)*cos(mission.launchBase.lonInitial ); cos(mission.launchBase.latInitial)*sin(mission.launchBase.lonInitial ); sin(mission.launchBase.latInitial)]' ;

mission.target.latInitial = - mission.launchBase.latInitial;
mission.target.lonInitial = mission.launchBase.lonInitial - pi ;
mission.target.initialPointECI = 6371000*[cos(mission.target.latInitial)*cos(mission.target.lonInitial); cos(mission.target.latInitial)*sin(mission.target.lonInitial); sin(mission.target.latInitial) ];


n = cross(mission.launchBase.initialPointECI,[0,1,0])/(norm(cross(mission.launchBase.initialPointECI,[0,1,0])));
ex = mission.launchBase.initialPointECI' / norm(mission.launchBase.initialPointECI);
ez = n' / norm(n) ;
ey = cross(ez,ex)/norm(cross(ez,ex));
rot = [ex,ey,ez]';

mission.target.Rfinal = rot ;

% ============================ Environment ===============================

mission.environment.rEarth = 6371e3;
mission.environment.g0 = 9.81;
mission.environment.GM = 398600.433e9;
mission.environment.rhoFun = @(h) 1.29*exp(-h/8433);
mission.environment.omega = 2*pi/86164.1 ;

hVec = 0:100;
meanWind = GRAM07_HWM07_annual(hVec);
windAngVel = meanWind ./ (mission.environment.rEarth + hVec);
vxWind = - windAngVel .* (mission.environment.rEarth + hVec) .* sin(mission.launchBase.lonInitial) ;
vyWind = windAngVel .* (mission.environment.rEarth + hVec) .* cos(mission.launchBase.lonInitial) ;
mission.environment.windXFun = griddedInterpolant(hVec,vxWind,'linear','linear');
mission.environment.windYFun = griddedInterpolant(hVec,vyWind,'linear','linear');

soundSpeedVec = load("soundspeed.mat");

soundSpeedVec = soundSpeedVec.soundspeed(1:70:end);
alt = linspace(0,500000,length(soundSpeedVec));


mission.aerodynamics.soundspeedFun= griddedInterpolant( ...
                                    alt, ...   % grid points
                                    soundSpeedVec, ...        % values
                                    'linear', ...                       % interpolation method
                                    'linear');  


[~,~,pressure] = atmosisa(alt);
pressure = pressure.*(pressure>0.5);

mission.environment.pressure=   griddedInterpolant( ...
                                alt, ...   % grid points
                                pressure, ...        % values
                                'linear', ...                       % interpolation method
                                'linear');  

% ============================ Aerodynamics ===============================



mission.aerodynamics.finsGeom.rootChord   = 1.81;                      % root chord [m]
mission.aerodynamics.finsGeom.tipChord   = 0.45;                       % tip chord [m] (triangolo puro)
mission.aerodynamics.finsGeom.bfin    = 1.81;                      % semispan [m]

mission.aerodynamics.finsGeom.Nfins = 2;     
mission.aerodynamics.finsGeom.Sfin = 0.5 * (mission.aerodynamics.finsGeom.rootChord+mission.aerodynamics.finsGeom.tipChord) * mission.aerodynamics.finsGeom.bfin;       % area della fin [m^2]
mission.aerodynamics.finsGeom.lambda = mission.aerodynamics.finsGeom.tipChord / mission.aerodynamics.finsGeom.rootChord;
mission.aerodynamics.finsGeom.cmac = (2/3) * mission.aerodynamics.finsGeom.rootChord * (1 + mission.aerodynamics.finsGeom.lambda + mission.aerodynamics.finsGeom.lambda^2) / (1 + mission.aerodynamics.finsGeom.lambda);       % mean aerodynamic chord [m]

mission.aerodynamics.finsGeom.delta_le = 30; % [deg]
mission.aerodynamics.finsGeom.Lambda_le = 36.9; % [deg]

mission.aerodynamics.finsGeom.tmac = 0.1;    



mission.aerodynamics.bodyInfo.isPowered = true;
mission.aerodynamics.bodyInfo.a_sub = 0;
mission.aerodynamics.bodyInfo.b_sub = 1;
mission.aerodynamics.bodyInfo.Cdn = 1.2;

% ============================ Materials and structures ===============================

% Al2219 - cryogenic tanks and primary structures for 1st/2nd stage
% mission.materials.Al2219.rho      = 2840;
% mission.materials.Al2219.E        = 72e9;
% mission.materials.Al2219.yield    = 390e6;
% mission.materials.Al2219.ultimate = 480e6;

% {
% ----------------------------------------------------------------------
% Al-Li 2195 - lightweight cryogenic tanks for high-performance stages 
% ----------------------------------------------------------------------
mission.structure.rho      = 2710; %https://www.matweb.com/search/datasheet_print.aspx?matguid=4363dafc7f5545688506d8b4af1e9468&utm_source
mission.structure.E        = 75e9; %https://it.scribd.com/document/922715788/Al-Li-2195-Mechanical-Properties-Research?utm_source
mission.structure.yield    = 500e6; %https://it.scribd.com/document/922715788/Al-Li-2195-Mechanical-Properties-Research?utm_source
mission.structure.ultimate = 560e6; %https://it.scribd.com/document/922715788/Al-Li-2195-Mechanical-Properties-Research?utm_source
% }

%{
%----------------------------------------------------------------------
% Al 2024-T3 - structural panels, frames, non-cryogenic structures
%----------------------------------------------------------------------
mission.structure.rho      = 2780;
mission.structure.E        = 73e9;
mission.structure.yield    = 345e6;
mission.structure.ultimate = 483e6;
%}

%{
%----------------------------------------------------------------------
% Al 7075-T6 - highly loaded fittings, secondary structures
%----------------------------------------------------------------------
mission.structure.rho      = 2810;
mission.structure.E        = 72e9;
mission.structure.yield    = 503e6;
mission.structure.ultimate = 572e6;
%}


mission.structure.rho = [2840, 2710, 2780, 2810];
mission.structure.E   = [72e9, 75e9, 73e9, 72e9];
mission.structure.ultimate = [480e6, 560e6, 483e6, 572e6];
mission.structure.yield = [390e6, 500e6, 345e6, 503e6];
mission.structure.minThickness = 1e-3 ;

mission.structure.tankPressure       = 3.4e5;
mission.structure.safetyFactor       = 1.4;



%% settings for functions


settings = struct();

settings.globalGAoptVariables = 7 ;
settings.nOptPointsTraj = 5 ;

% ============================ Initial population internal GA ===============================

settings.initialPopulationGATraj.oneStage = load("InitialPop1stages.mat") ;
settings.initialPopulationGATraj.oneStage = settings.initialPopulationGATraj.oneStage.X ;
clear settings.initialPopulationGATraj.oneStage.X ;


settings.initialPopulationGATraj.twoStages = load("InitialPop2stages.mat") ;
settings.initialPopulationGATraj.twoStages = settings.initialPopulationGATraj.twoStages.X ;
clear settings.initialPopulationGATraj.twoStages.X ;


settings.initialPopulationGATraj.threeStages = load("InitialPop3stages.mat") ;
settings.initialPopulationGATraj.threeStages = settings.initialPopulationGATraj.threeStages.X ;
clear settings.initialPopulationGATraj.threeStages.X ;




% ============================ Lower and upper bounds ===============================

settings.lowerBoundsFMC(:,:,1) = [0.9*ones(settings.nOptPointsTraj,1);...
                        0*ones(settings.nOptPointsTraj,1)];
settings.lowerBoundsFMC(:,:,2) = [0.4*ones(settings.nOptPointsTraj,1);...
                        0*ones(settings.nOptPointsTraj,1)];
settings.lowerBoundsFMC(:,:,3) = [0.4*ones(settings.nOptPointsTraj,1);...
                        0*ones(settings.nOptPointsTraj,1)];


settings.upperBoundsFMC(:,:,1) = [ones(settings.nOptPointsTraj,1);...
                        90*ones(settings.nOptPointsTraj,1)];
settings.upperBoundsFMC(:,:,2) = [ones(settings.nOptPointsTraj,1);...
                        120*ones(settings.nOptPointsTraj,1)];
settings.upperBoundsFMC(:,:,3) = [ones(settings.nOptPointsTraj,1);...
                        150*ones(settings.nOptPointsTraj,1)];


settings.lowerBoundsGA = settings.lowerBoundsFMC(:);
settings.upperBoundsGA = settings.upperBoundsFMC(:);


settings.lowerBoundsGlobalGA = [2,1,3,3,0.35,0.3,0.3];
settings.upperBoundsGlobalGA = [3,4,3,3,0.65,0.9,0.9];


% ============================ Function options ===============================

settings.gaTrajOptions = optimoptions("ga", ...
                        "Display","iter", ...
                        "MaxGenerations",20, ...
                        "PopulationSize",500,...
                        "UseParallel",false,...
                        "FunctionTolerance", 1e-4);


settings.fminconTrajOptions = optimoptions("fmincon",...
                                "Display","iter",...
                                "MaxIterations",100,...
                                'MaxFunctionEvaluations',10000,...
                                'StepTolerance',1e-15,...
                                'OptimalityTolerance',1e-6,...
                                'FunctionTolerance',1e-15,...
                                'ConstraintTolerance',1e-10, ...
                                'EnableFeasibilityMode',true,...
                                "UseParallel",false);


settings.globalGAOptions = optimoptions("ga", ...
                                "Display","iter", ...
                                "MaxGenerations",20, ...
                                "PopulationSize",48,...
                                "EliteCount", 1,...
                                "UseParallel",true,...
                                "FunctionTolerance", 1e-4, ...
                                "NonlinearConstraintAlgorithm", "penalty", ...
                                "PlotFcn",{"gaplotbestf", "gaplotbestindiv"});


settings.fsolveOptions = optimoptions('fsolve','Display','none');



% ============================ Other optimization parameters ===============================

settings.nEvalPointsTraj = 100 ;

settings.intconGlobalGA = [1 2 3 4];
settings.trajectoryOption2D = 1;


%============================ Plots ======================================

settings.color.green = '#777754';
settings.color.orange = '#cf9143';
settings.color.blu = '#2f5d62';
settings.color.gray = '#3a3a3a';
settings.color.terracotta = '#8f4e3a';

end