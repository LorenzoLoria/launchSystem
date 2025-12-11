function [tt,xx] = launcherTrajectoryControl(x0,mission, mer, configuration, launcher, ~,tSpan,nDeval, stageNumber, guidancePoints, guidanceTime, gains)
    
options = odeset('RelTol',1e-6,'AbsTol',1e-6 ,'Events',@(t,x) propEvent(t,x, mission,configuration.stage{stageNumber}.mProp,configuration,stageNumber));


solution = ode113(@(t,x) launcherDynamicsAndControlECI(t, x, mission, mer, configuration, launcher, stageNumber, guidancePoints, guidanceTime, gains), tSpan, x0, options);

tt = linspace(solution.x(1),solution.x(end), nDeval);

xx = deval(solution,tt);

end