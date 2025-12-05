function [value,isterminal,direction] = touchDown(t,x)
    value      = x(3); 
    isterminal = 1;       
    direction  = 0;      
end