function [mission,settings] = dataStructGlobal

mission = struct();

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

mission.engines{4}.ispZero = 366;
mission.engines{4}.ispVac = 452;
mission.engines{4}.thrustZero = 1860*1e3;
mission.engines{4}.thrustVacum = 2279*1e3;
mission.engines{4}.weight = 3526;
mission.engines{3}.OF = 6;
mission.engines{4}.oxDens = 1143;
mission.engines{4}.fuelDens = 71;
mission.engines{4}.effAreaZero = 2.3^2/4*pi;
mission.engines{4}.effAreaVac = 2.3^2/4*pi;


mission.capsule.weight = 8600;
mission.capsule.Area = 4^2*pi/4;
mission.capsule.supersonicCD = 1.23;
mission.capsule.subsonicCD = 0.45;
mission.capsule.height = 2.9;


mission.environment.rEarth = 6371e3;
mission.environment.g0 = 9.81;
mission.environment.GM = 398600.433e9;
mission.environment.rhoFun = @(h) 1.29*exp(-h/8433);
mission.environment.omega = 2*pi/86164.1 ;



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


% ============================ Aerodynamics ===============================

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



% Al2219 - cryogenic tanks and primary structures for 1st/2nd stage
mission.materials.Al2219.rho      = 2840;
mission.materials.Al2219.E        = 72e9;
mission.materials.Al2219.yield    = 390e6;
mission.materials.Al2219.ultimate = 480e6;

% Pressurizzazione serbatoi
mission.structure.tankPressure       = [0, 0, 0, 3.4e5, 3.4e5, 3.4e5, 3.4e5, 0];


%----------------------------------------------------------------------
% Al 7075-T6 - highly loaded fittings, secondary structures
%----------------------------------------------------------------------
mission.materials.Al7075.rho      = 2720;
mission.materials.Al7075.E        = 75e9;
mission.materials.Al7075.yield    = 500e6;
mission.materials.Al7075.ultimate = 560e6;



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

settings.fsolveOptions = optimoptions('fsolve','Display','none');

settings.nOptPointsTraj = 5 ;
settings.nEvalPointsTraj = 100 ;



end