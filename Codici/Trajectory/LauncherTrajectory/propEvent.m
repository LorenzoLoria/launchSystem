function [value,isterminal,direction] = propEvent(t,x, mission,mP,opt,stageNumber)

m0 = opt.m0Tot;
h = norm(x(1:3)) - mission.environment.rEarth + 0.1*stageNumber;


value(1) = x(end) - (m0 - mP);

isterminal(1) = 1;       
direction(1)  = -1;     

value(2) = h;

isterminal(2) = 1;       
direction(2)  = -1;  
if isnan(value)
    keyboard
end
if isinf(value)
    keyboard
end

end