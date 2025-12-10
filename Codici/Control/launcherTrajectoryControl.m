function [tt,xx] = launcherTrajectoryControl(x0,mission, configuration, launcher, thrustData,tSpan,nDeval, guidancePoints, guidanceTime, gains)
    

options = odeset('RelTol',1e-6,'AbsTol',1e-4 ,'Events');

solution = ode113(@(t,x) auncherDynamicsAndControlECI(t, x,thrustData, mission, configuration, launcher, stageNumber, guidancePoints, guidanceTime, gains), tSpan, x0,options);

tt = linspace(solution.x(1),solution.x(end), nDeval);

xx = deval(solution,tt);

end