function [tt,xx] = launcherTrajectory(x0,mission,thrustData,t0,nDeval,stageNumber,opt)

%persistent copyThrustData copytt copyxx

%if isempty(copyThrustData)
%copyThrustData = zeros(mission.optimisation.GA.variables,2);
%end

%if thrustDataVec == copyThrustData 
%tt = copytt;
%xx = copyxx;
%else 

% ThrustDataVec deve essere una matrice n*3

tSpan = [t0 3*24*3600];

% Manca la propEvent che deve riferirsi a quando finisce il carburante
% dello stadio o qualcosa di simile.

options = odeset('RelTol',1e-6,'AbsTol',1e-2 ,'Events',@(t,x) propEvent(t,x, mission,opt.stage{stageNumber}.mProp,opt));

solution = ode113(@(t,x) launcherDynamicsECI(t, x,thrustData, mission,stageNumber,opt), tSpan, x0,options);

tt = linspace(solution.x(1),solution.x(end),nDeval);

xx = deval(solution,tt);

%copytt = tt;
%copyxx = xx;
%copyThrustData = thrustDataVec;
%end
end