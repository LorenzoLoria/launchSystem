function [value,isterminal,direction] = propEvent(t,x, mission)

mP = mission.launcher.engines{1}.mPropellant;
m0 = mission.launcher.engines{1}.m0;

value = x(end) - (m0 - mP);

isterminal = 1;       
direction  = -1;     


end