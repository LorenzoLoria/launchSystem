function [tt,xx] = ballisticTrajectory(x0,mission,windDirection,t0)

tSpan = [t0 3*24*3600];

options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@groundEvent);

solution = ode113(@(t,x) dynCapsule(t,x,mission,windDirection), tSpan, x0,options);

tt = solution.x;

xx = solution.y;

end