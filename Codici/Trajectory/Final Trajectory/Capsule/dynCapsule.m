function [dxdt] = dynCapsule(~, x,mission, windDirection)
%   KEPLERIAN_RHS  Evaluates the right-hand-side of a 2-body (keplerian) propagator
%   Evaluates the right-hand-side of a newtonian 2-body propagator.
%   x is the state


GM = mission.environment.GM;
h =norm(x(1:3))-mission.environment.rEarth; 
rho = mission.environment.rhoFun(h);

v = x(4:6) ;
%dynamicPressure = 0.5 * rho * norm(v)^2;
%[soundspeed] = soundSpeedFun(h);
%Mach = norm(v)/soundspeed;
%[CL,Cd,~,~] = CLCDcomputation(Mach,0,dynamicPressure,0,mission);
Cd = mission.capsule.supersonicCD;
aeroForce = - 0.5 .* rho .* norm((v)) .* v .* mission.capsule.Area .* Cd ;

% Initialize right-hand-side
dxdt = zeros(6,1);

% Extract positions
rr = x(1:3);

% Position detivative is object's velocity
dxdt(1:3) = x(4:6);   

% Compute the gravitational acceleration using Newton's law + air drag
dxdt(4:6) = - GM * rr /(dot(rr, rr)^(3/2)) + aeroForce./mission.capsule.weigth;

end






