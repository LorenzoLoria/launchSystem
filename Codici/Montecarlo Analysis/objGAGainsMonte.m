function [objective] = objGAGainsMonte(gainGA,launcher,opt,mission,settings,windVelXFun,windVelYFun,stateCollocationRef,timeCollocationRef,maxThrust,maxGimball,thrustData)

[timeCollocation, stateCollocation] = totalTrajectoryControlled(launcher,opt,mission,settings,windVelXFun,windVelYFun,stateCollocationRef,timeCollocationRef,maxThrust,maxGimball,thrustData,gainGA);

objective = norm(stateCollocation(1:3,end,end) - mission.target.initialPointECI);

end