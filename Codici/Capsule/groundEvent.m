function [value,isterminal,direction] = groundEvent(t,x)
    h = norm(x(1:3))-6371e3;             
    value      = h; 
    isterminal = 1;       
    direction  = -1;      
end
