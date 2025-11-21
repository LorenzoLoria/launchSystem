function [value,isterminal,direction] = propEvent(t,x, mission,mP)

m0 = mission.launcher.engines{1}.m0;
h = norm(x(1:3)) - mission.environment.rEarth + 0.2;


value(1) = x(end) - (m0 - mP);

isterminal(1) = 1;       
direction(1)  = -1;     

value(2) = h;

isterminal(2) = 1;       
direction(2)  = -1;  

end