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
mission.capsule.Area = 3.7^2*pi;
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


mission.options.fmincon = optimoptions("fmincon","Display","iter","MaxIterations",200,'MaxFunctionEvaluations',10000,'StepTolerance',1e-15,'OptimalityTolerance',1e-8,'FunctionTolerance',1e-19,'ConstraintTolerance',1e-10,'EnableFeasibilityMode',true,"UseParallel",true);
mission.options.gaOptions = optimoptions('ga', 'PlotFcn',{'gaplotbestf', 'gaplotbestindiv'}, 'display', 'iter','MaxStallGenerations', 10, ...
        'FunctionTolerance', 1e-6, 'EliteCount',  6,...
        'MaxGenerations', 100, 'PopulationSize', 40, ...
        'NonlinearConstraintAlgorithm', 'penalty',"UseParallel",true);


% Robe di prova con il booster

mission.launcher.booster.mass = 5600; 
mission.launcher.booster.cp = [0 0 -2]';

mission.environment.rhoFun = @(h) 1.29*exp(-h/8433);


lat = deg2rad(-39.261515237910196);
lon =deg2rad(+177.86521705520965-180);
mission.initialPoint = 6371000*[cos(lat)*cos(lon); cos(lat)*sin(lon); sin(lat)]';

%mission.initialPoint = [0 0 6371000];

latFinal = deg2rad(+39.261515237910196);
lonFinal =deg2rad(+177.86521705520965);
mission.target = 6371000*[cos(latFinal)*cos(lonFinal); cos(latFinal)*sin(lonFinal); sin(latFinal) ];

n = cross(mission.initialPoint,[1,0,0])/(norm(cross(mission.initialPoint,[1,0,0])));


ex = mission.initialPoint' / norm(mission.initialPoint);
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
optimisation.stage{1}.Isp = mission.launcher.engines{1}.isp;

optimisation.stage{2}.mStage = mStage2;
optimisation.stage{2}.mProp = mProp2;
optimisation.stage{2}.Isp = mission.launcher.engines{3}.isp;









end
