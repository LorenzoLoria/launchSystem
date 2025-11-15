function [dxdt] = dynCapsule(~, x,mission)

%   KEPLERIAN_RHS  Evaluates the right-hand-side of a 2-body (keplerian) propagator
%   Evaluates the right-hand-side of a newtonian 2-body propagator.
%   x is the state


GM = mission.environment.GM;

rho = interp1(mission.envirnoment.altRange,mission.envirnoment.rho,norm(x(1:3))-mission.envirnoment.rEarth);

v = norm(x(4:6));

D = - 0.5 * rho * v^2 * mission.capsule.Area * mission.capsule.Cd .* x(4:6)./v;

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
dxdt(4:6) = - GM * rr /(dist*dist2) + D./mission.capsule.weigth;

end






