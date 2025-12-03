%% === SCEGLI IL MATERIALE PER LA STRUTTURA ===
% Decommenta SOLO uno dei blocchi qui sotto

%----------------------------------------------------------------------
% Al2219 - cryogenic tanks and primary structures for 1st/2nd stage (LOX/LH2, LOX/RP-1)
%----------------------------------------------------------------------
mission.structure.rho      = 2840;
mission.structure.E        = 72e9;
mission.structure.yield    = 390e6;
mission.structure.ultimate = 480e6;

%{
%----------------------------------------------------------------------
% Al-Li 2195 - lightweight cryogenic tanks for high-performance stages 
%----------------------------------------------------------------------
mission.structure.rho      = 2720;
mission.structure.E        = 75e9;
mission.structure.yield    = 500e6;
mission.structure.ultimate = 560e6;
%}

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

%{
%----------------------------------------------------------------------
% Ti-6Al-4V (Grade 5) - thrust rings, engine mounts, highly loaded nodes
%----------------------------------------------------------------------
mission.structure.rho      = 4430;
mission.structure.E        = 112e9;
mission.structure.yield    = 880e6;
mission.structure.ultimate = 950e6;
%}

%{
%----------------------------------------------------------------------
% CFRP UD carbon/epoxy (fibra longitudinale)
%----------------------------------------------------------------------
mission.structure.rho      = 1600;
mission.structure.E        = 130e9;
mission.structure.yield    = 1000e6;
mission.structure.ultimate = 1500e6;
%}

%{
%----------------------------------------------------------------------
% Maraging steel C250 - solid rocket motor casings, high-strength parts
%----------------------------------------------------------------------
mission.structure.rho      = 8000;
mission.structure.E        = 195e9;
mission.structure.yield    = 1724e6;
mission.structure.ultimate = 1758e6;
%}

%{
%----------------------------------------------------------------------
% D6AC steel - large solid booster casings, critical structures
%----------------------------------------------------------------------
mission.structure.rho      = 7850;
mission.structure.E        = 210e9;
mission.structure.yield    = 1724e6;
mission.structure.ultimate = 1930e6;
%}

%{
%----------------------------------------------------------------------
% Inconel 718 - hot structures near nozzle, turbopump support
%----------------------------------------------------------------------
mission.structure.rho      = 8250;
mission.structure.E        = 200e9;
mission.structure.yield    = 1035e6;
mission.structure.ultimate = 1240e6;
%}
