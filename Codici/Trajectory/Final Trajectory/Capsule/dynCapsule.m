function [dxdt] = dynCapsule(~, x,mission, windDirection)
%   KEPLERIAN_RHS  Evaluates the right-hand-side of a 2-body (keplerian) propagator
%   Evaluates the right-hand-side of a newtonian 2-body propagator.
%   x is the state


% Vettori Stato
rr = x(1:3);
v  = x(4:6); 
    
% Calcoli della Norma (Fatti SOLO UNA volta)
rMag = sqrt(rr'*rr); 
vMag = sqrt(v'*v);

GM = mission.environment.GM;
h =rMag-mission.environment.rEarth; 
%rho = mission.environment.rhoFun(h);
rho = 1.29*exp(-h/8433);

soundspeed = mission.aerodynamics.soundspeedFun(h);
Mach = vMag/soundspeed;
%Cd = Cd_CrewDragon(Mach);

%Evaluate Cd of crewDragon
Cd = (0.68*Mach^2+0.424)*(Mach<= 0.722) + (1.73 * exp(-1.02 * (Mach + 0.25)) - 8.9 * exp(-1.95 * (Mach + 0.35)) + 1.23)*(Mach> 0.722);

aeroForce = - 0.5 .* rho .* vMag .* v .* mission.capsule.Area .* Cd ;

% Initialize right-hand-side
dxdt = zeros(6,1);

% Position detivative is object's velocity
dxdt(1:3) = x(4:6);   

% Compute the gravitational acceleration using Newton's law + air drag
dxdt(4:6) = - GM * rr /(rMag^3) + aeroForce./mission.capsule.weight;

end






