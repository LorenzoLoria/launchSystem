
function [value,isterminal,direction] = arrestingEvents(s,m)

    
    % finalIdx = sVec(end)==s ;
    % idx = sum((sVec-s)<=0) - finalIdx ;
    % sInitial = sVec(idx) ;
    % sFinal = sVec(idx+1) ;  
    % 
    % pos = position.position ; 
    % pos = (pos(:,idx+1) - pos(:,idx)) / (sFinal - sInitial) * (s-sInitial) + pos(:,idx) ;


    value = [ m ];
    isterminal = [1];       
    direction  = [0];    


end