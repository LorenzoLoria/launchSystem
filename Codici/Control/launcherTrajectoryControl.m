function [tt,xx] = launcherTrajectoryControl(x0,mission,guidancePoints,thrustData,tSpan,nDeval,stageNumber,opt)


options = odeset('RelTol',1e-6,'AbsTol',1e-4 ,'Events');

solution = ode113(@(t,x) launcherDynamicsAndControlECI(t, x, mission,configuration, thrustData, guidancePoints,gains), tSpan, x0,options);

tt = linspace(solution.x(1),solution.x(end),nDeval);

xx = deval(solution,tt);

end