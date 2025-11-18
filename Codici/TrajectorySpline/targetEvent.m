
function [distance,isterminal,direction] = targetEvent(s,~,position,x_target)

    position = position.position ; 
    distance = norm(position(s) - x_target') - 0.1 ;
    isterminal = 1;       
    direction  = 0;    

end