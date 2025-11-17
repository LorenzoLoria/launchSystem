function [tt,xx] = launcherTrajectory(x0,mission,thrustData)

% ThrustDataVec deve essere una matrice n*3

tSpan = [0 3*3600];

% Manca la propEvent che deve riferirsi a quando finisce il carburante
% dello stadio o qualcosa di simile.

options = odeset('RelTol',1e-6,'AbsTol',1e-2 ,'Events',@(t,x) propEvent(t,x, mission));

solution = ode113(@(t,x) launcherDynamicsECI(t, x,thrustData, mission), tSpan, x0,options);

tt = linspace(solution.x(1),solution.x(end),100);

xx = deval(solution,tt);
end