function [tt,xx] = ballisticTrajectory(x0,mission,windDirection,t0,nDeval)

persistent x0Copy copyttCaps copyxxCaps

if isempty(x0Copy)
x0Copy = zeros(6,1);
end

if x0 == x0Copy 
tt = copyttCaps;
xx = copyxxCaps;
else 

tSpan = [t0 3*24*3600];

options = odeset('RelTol',1e-6,'AbsTol',1e-2,'Events',@groundEvent);

solution = ode113(@(t,x) dynCapsule(t,x,mission,windDirection), tSpan, x0,options);

tt = linspace(solution.x(1),solution.x(end),nDeval);

xx = deval(solution,tt);

copyttCaps = tt;
copyxxCaps = xx;
x0Copy = x0;
end
end