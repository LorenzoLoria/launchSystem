% data script
function [mission,optimisation] = dataStruct

mission = struct();

mission.launcher.engines{1}.name = 'Merlin1D' ;
mission.launcher.engines{2}.name = 'Raptor' ;
mission.launcher.engines{3}.name = 'RS-25' ;

mission.launcher.engines{1}.isp = 283;
mission.launcher.engines{1}.thrust = 854*1e3;
mission.launcher.engines{1}.weight = 470;
mission.launcher.engines{1}.OF = 2.36;
mission.launcher.engines{1}.oxDens = 1143;
mission.launcher.engines{1}.fuelDens = 835;   

% Questi valori dovrebbero essere dati da GA:
mission.launcher.engines{1}.mPropellant1 =  200e3;
mission.launcher.engines{1}.ms1 = 12e3;
mission.launcher.engines{1}.mPropellant2 = 37.6e3;
mission.launcher.engines{1}.ms2 = 3e3;
mission.launcher.engines{1}.m0 =208e3;% mission.launcher.engines{1}.mPropellant1 + mission.launcher.engines{1}.ms1 + mission.launcher.engines{1}.mPropellant2 + mission.launcher.engines{1}.ms2 + 8600;

mission.launcher.engines{2}.isp = 327;
mission.launcher.engines{2}.thrust = 2750*1e3;
mission.launcher.engines{2}.weight = 1525;
mission.launcher.engines{2}.OF = 3.6;
mission.launcher.engines{2}.oxDens = 1143;
mission.launcher.engines{2}.fuelDens = 423;

mission.launcher.engines{3}.isp = 366;
mission.launcher.engines{3}.thrust = 1860*1e3;
mission.launcher.engines{3}.weight = 3526;
mission.launcher.engines{3}.OF = 6;
mission.launcher.engines{3}.oxDens = 1143;
mission.launcher.engines{3}.fuelDens = 71;

mission.capsule.weigth = 8600;
mission.capsule.Area = 4^2*pi/4;
mission.capsule.supersonicCD = 1.23;
mission.capsule.subsonicCD = 0.45;

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


mission.options.fmincon = optimoptions("fmincon","Display","iter","MaxIterations",200,'MaxFunctionEvaluations',10000,'StepTolerance',1e-19,'OptimalityTolerance',1e-6,'FunctionTolerance',1e-15,'ConstraintTolerance',1e-10,'EnableFeasibilityMode',true,"UseParallel",true);
mission.options.gaOptions = optimoptions('ga', 'PlotFcn',{'gaplotbestf', 'gaplotbestindiv'}, 'display', 'iter','MaxStallGenerations', 10, ...
        'FunctionTolerance', 1e-6, 'EliteCount',  6,...
        'MaxGenerations', 100, 'PopulationSize', 40, ...
        'NonlinearConstraintAlgorithm', 'penalty',"UseParallel",true);


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

optimisation.stage{1}.nEngines = 5;
optimisation.stage{2}.nEngines = 1;

mProp1 = optimisation.stage{1}.nEngines*optimisation.stage{1}.engine.thrust/9.81*0.6;
epsS1 = 0.05;

mProp2 = optimisation.stage{2}.nEngines*optimisation.stage{2}.engine.thrust/9.81*0.4;
epsS2 = 0.06;

mS1 = epsS1/(1-epsS1)*mProp1;
mS2 = epsS2/(1-epsS2)*(mProp2);

mStage1 = mProp1+mS1;
mStage2 = mProp2+mS2;


optimisation.m0Tot = mStage1 + mStage2 + mission.capsule.weigth;

optimisation.stage{1}.mStage = mStage1;
optimisation.stage{1}.mProp = mProp1;
optimisation.stage{1}.structuralMass = mS1;
optimisation.stage{1}.Isp = mission.launcher.engines{1}.isp;

optimisation.stage{2}.mStage = mStage2;
optimisation.stage{2}.mProp = mProp2;
optimisation.stage{2}.structuralMass = mS2;
optimisation.stage{2}.Isp = mission.launcher.engines{3}.isp;







% Aerodynamics

mission.aerodynamics.bodyGeom.l = 3;
mission.aerodynamics.bodyGeom.d = 0.3;
mission.aerodynamics.bodyGeom.Lnose = 0.6;
mission.aerodynamics.bodyGeom.Aref = pi*(bodyGeom.d/2)^2; 
mission.aerodynamics.bodyGeom.Anose = mission.aerodynamics.bodyGeom.Aref;
mission.aerodynamics.bodyGeom.Abase = mission.aerodynamics.bodyGeom.Aref;
mission.aerodynamics.bodyGeom.Aexit = 0;  
mission.aerodynamics.bodyGeom.phi = deg2rad(0);
mission.aerodynamics.bodyGeom.Ab = mission.aerodynamics.bodyGeom.Aref;
mission.aerodynamics.bodyGeom.Ap = mission.aerodynamics.bodyGeom.l * mission.aerodynamics.bodyGeom.d;



cr   = 0.35;                      % root chord [m]
ct   = 0.0;                       % tip chord [m] (triangolo puro)
s    = 0.20;                      % semispan [m]

mission.aerodynamics.finsGeom.Nfins = 1;
mission.aerodynamics.finsGeom.be = s;                  % span equivalente ~ semispan reale
mission.aerodynamics.finsGeom.Se = 0.5 * cr * s;       % area della fin [m^2]
mission.aerodynamics.finsGeom.cmac = (2/3) * cr;       % mean aerodynamic chord [m]
mission.aerodynamics.finsGeom.delta_le = 45;
mission.aerodynamics.finsGeom.lambda_le = 0;
mission.aerodynamics.finsGeom.b = 2 * cr;
mission.aerodynamics.finsGeom.tmac = 0.08 * finsGeom.cmac;    



mission.aerodynamics.bodyInfo.isPowered = false;
mission.aerodynamics.bodyInfo.a_sub = 0;
mission.aerodynamics.bodyInfo.b_sub = 1;
mission.aerodynamics.bodyInfo.Cdn = 1.2;



mission.aerodynamics.finsInfo.rho = 1.225;
mission.aerodynamics.finsInfo.vsound = 340;


end
