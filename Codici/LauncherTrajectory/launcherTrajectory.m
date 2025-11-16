function [tt,xx] = launcherTrajectory(x0,mission,thrustDataVec)

% ThrustDataVec deve essere una matrice n*3

tSpan = [0 3*24*3600];

% Manca la propEvent che deve riferirsi a quando finisce il carburante
% dello stadio o qualcosa di simile.

options = odeset('RelTol',1e-12,'AbsTol',1e-12); % ,'Events',@propEvent);

thrustData = @(t) [interp1(tVec,thrustDataVec(:,1),t);interp1(tVec,thrustDataVec(:,2),t);interp1(tVec,thrustDataVec(:,3),t)];

solution = ode45(@(t,x) launcherDynamicsECI(t, x,@(t) thrustData, mission), tSpan, x0,options);

tt = linspace(solution.x(1), solution.x(end), 1000);

xx = deval(solution,tt);