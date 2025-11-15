function [dxdt] = dynCapsule(~, x,mission, windDirection)

%   KEPLERIAN_RHS  Evaluates the right-hand-side of a 2-body (keplerian) propagator
%   Evaluates the right-hand-side of a newtonian 2-body propagator.
%   x is the state


GM = mission.environment.GM;
rho = interp1(mission.environment.altRange,mission.environment.rho,norm(x(1:3))-mission.environment.rEarth);

%Modello vento preso da MIMP
h = norm(x(1:3)) / 1e3 - mission.environment.rEarth / 1e3; %[km]
windIntensity = (6.9288 * h + 9.144).*(h<9.6) + 76.2 .* (h>=9.6 && h<14) +...
                (76.2-8.9474 * (h-14)) .* (h>=14 && h<20) + 24.384 .* (h>=20) ; 


windVelocity = x(4:6) - windIntensity * windDirection ;

aeroForce = - 0.5 .* rho .* norm((windVelocity)) .* windVelocity .* mission.capsule.Area .* mission.capsule.supersonicCD ;

% Initialize right-hand-side
dxdt = zeros(6,1);

% Extract positions
rr = x(1:3);

% Compute square distance and distance
dist2 = dot(rr, rr);
dist = sqrt(dist2);

% Position detivative is object's velocity
dxdt(1:3) = x(4:6);   

% Compute the gravitational acceleration using Newton's law + air drag
dxdt(4:6) = - GM * rr /(dist*dist2) + aeroForce./mission.capsule.weigth;

end






