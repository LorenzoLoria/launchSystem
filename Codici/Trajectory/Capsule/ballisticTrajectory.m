function [tt,xx] = ballisticTrajectory(x0,mission,windDirection,t0)

tSpan = [t0 3*3600];

options = odeset('RelTol',1e-6,'AbsTol',1e-2,'Events',@groundEvent);

solution = ode45(@(t,x) dynCapsule(t,x,mission,windDirection), tSpan, x0,options);

tt = linspace(solution.x(1),solution.x(end),100);

xx = deval(solution,tt);

end