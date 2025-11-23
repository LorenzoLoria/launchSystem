function[objective] = objFunMultiStagesGA(x,mission,opt)

thrustDataVec = x;
%thrustDataVec is the variable containing informations about the
%optimisation paramets. As of right now, the trajectory optimisation
%variables are throttoling and angle. The best way to divide them is by
%using a optVarNum x 2 x nStages matrix, since it will help in keeping the
%code organised correctly.

[timeCollocation,stateCollocation] = totalTrajectory(mission,opt,thrustDataVec);

objective =    timeCollocation(end,end);

end