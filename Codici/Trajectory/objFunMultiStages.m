function[objective] = objFunMultiStages(x,mission,opt)

thrustDataVec = x;
%thrustDataVec is the variable containing informations about the
%optimisation paramets. As of right now, the trajectory optimisation
%variables are throttoling and angle. The best way to divide them is by
%using a optVarNum x 2 x nStages matrix, since it will help in keeping the
%code organised correctly.

thrustData = cell(opt.nStages);

for i = 1:nStages
    % Needs to be changed with the relative times
    tSpan = [0 3*24*3600];
    tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));    
    fThrust = griddedInterpolant(tVec, thrustDataVec(:,:,i), 'linear', 'none'); 
    thrustFunction= @(t) fThrust(t).';
    thrustData{i} = thrustFunction;
end

[timeCollocation,~] = totalTrajectory(mission,opt,thrustDataVec);

objective = timeCollocation(end,end) ;

end