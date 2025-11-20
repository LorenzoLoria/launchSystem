
function [distance,isterminal,direction] = targetEvent(theta,~,position,target)

    persistent countTargetEvent
    if isempty(countTargetEvent)
        countTargetEvent = 1 ;
    end


    position = position.position ; 
    distance = norm(position(:,countTargetEvent) - target') - 5e3 ;
    isterminal = 1;       
    direction  = 0;    


    countTargetEvent = countTargetEvent + 1 ;

end