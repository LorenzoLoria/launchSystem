function [value,isterminal,direction] = guidanceEvent(t, ~, guidanceTime, iterNumber)
  
 value = t - guidanceTime(iterNumber);
 % Stop the integration when the event occurs
 isterminal = 1; 
 % All directions are considered
 direction = 0; 

end