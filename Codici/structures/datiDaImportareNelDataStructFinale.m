% ======================= STRUCTURES DATA =================================
% Density [kg/m^3], Young Modulus [Pa], Yield and Ultimate Stresses
% of different materials 

% Al2219 - cryogenic tanks and primary structures for 1st/2nd stage (LOX/LH2, LOX/RP-1)
mission.launcher.structures{1}.rho      = 2840;
mission.launcher.structures{1}.E        = 72e9;
mission.launcher.structures{1}.yield    = 390e6;
mission.launcher.structures{1}.ultimate = 480e6;

% Al-Li 2195 - lightweight cryogenic tanks for high-performance stages 
% (mainly 1st stage and upper stages)
mission.launcher.structures{2}.rho      = 2720;
mission.launcher.structures{2}.E        = 75e9;
mission.launcher.structures{2}.yield    = 500e6;
mission.launcher.structures{2}.ultimate = 560e6;

% Al 2024-T3 - structural panels, frames, non-cryogenic structures (1st/2nd stage)
mission.launcher.structures{3}.rho      = 2780;
mission.launcher.structures{3}.E        = 73e9;
mission.launcher.structures{3}.yield    = 345e6;
mission.launcher.structures{3}.ultimate = 483e6;

% Al 7075-T6 - highly loaded fittings, secondary structures (1st/2nd stage)
mission.launcher.structures{4}.rho      = 2810;
mission.launcher.structures{4}.E        = 72e9;
mission.launcher.structures{4}.yield    = 503e6;
mission.launcher.structures{4}.ultimate = 572e6;

% Ti-6Al-4V (Grade 5) - thrust rings, engine mounts, highly loaded nodes (mostly 1st stage)
mission.launcher.structures{5}.rho      = 4430;
mission.launcher.structures{5}.E        = 112e9;
mission.launcher.structures{5}.yield    = 880e6;
mission.launcher.structures{5}.ultimate = 950e6;

% CFRP UD carbon/epoxy (fiber direction) - interstages, fairings, cylindrical 
% structures, and solid motor casings (mainly upper stages)
mission.launcher.structures{6}.rho      = 1600;
mission.launcher.structures{6}.E        = 130e9;
mission.launcher.structures{6}.yield    = 1000e6;
mission.launcher.structures{6}.ultimate = 1500e6;

% Maraging steel C250 - solid rocket motor casings for 1st stage boosters,
% ultra-high-strength structural elements
mission.launcher.structures{7}.rho      = 8000;
mission.launcher.structures{7}.E        = 195e9;
mission.launcher.structures{7}.yield    = 1724e6;
mission.launcher.structures{7}.ultimate = 1758e6;

% D6AC steel - large solid booster casings and critical structural components (1st stage)
mission.launcher.structures{8}.rho      = 7850;
mission.launcher.structures{8}.E        = 210e9;
mission.launcher.structures{8}.yield    = 1724e6;
mission.launcher.structures{8}.ultimate = 1930e6;

% Inconel 718 - hot structures near the nozzle, turbopump support regions (1st stage)
mission.launcher.structures{9}.rho      = 8250;
mission.launcher.structures{9}.E        = 200e9;
mission.launcher.structures{9}.yield    = 1035e6;
mission.launcher.structures{9}.ultimate = 1240e6;
