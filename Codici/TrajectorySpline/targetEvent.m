
function [distance,isterminal,direction] = targetEvent(theta,~,position,vSpline, target)

    
    thetaVec = vSpline.profile(:,1);
    finalIdx = thetaVec(end)==theta ;
    idx = sum((thetaVec-theta)<=0) - finalIdx;
    thetaInitial = thetaVec(idx) ; 
    thetaFinal = thetaVec(idx+1) ; 

    pos = position.position ; 
    pos = (pos(:,idx+1) - pos(:,idx)) / (thetaFinal - thetaInitial) * (theta-thetaInitial) + pos(:,idx) ;


    distance = norm(pos - target') - 5e3 ;
    isterminal = 1;       
    direction  = 0;    


end