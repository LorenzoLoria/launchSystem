function [tt,xx] = ballisticTrajectory(x0,mission)

tSpan = [0 3*24*3600];

options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@groundEvent);

solution = ode45(@(t,x) dynCapsule(t,x,mission), tSpan, x0,options);

tt = linspace(solution.x(1), solution.x(end), 1000);

xx = deval(solution,tt);

end