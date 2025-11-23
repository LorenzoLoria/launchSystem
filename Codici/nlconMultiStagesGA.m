function[cin,ceq] = nlconMultiStagesGA(x,mission,opt)

thrustDataVec = x;

[timeCollocation,stateCollocation] = totalTrajectory(mission,opt,thrustDataVec);

A  = mission.capsule.Area;
cD = mission.capsule.supersonicCD;

for i = 1:opt.nStages
    time      = linspace(timeCollocation(1,i),timeCollocation(end,i),length(timeCollocation)-1);
    mdot      = diff(stateCollocation(end,:,i))./time;
    Thrust    = mdot .* opt.stage{i}.Isp .* misssion.environment.g0;
    Thrust    = [Thrust , Thrust(end)];
    r         = stateCollocation(1:3,:,i);
    rNorm     = sqrt( stateCollocation(1,:,i).^2 + stateCollocation(2,:,i).^2 +stateCollocation(3,:,i).^2 );
    h         = r - mission.environment.rEarth;
    v         = stateCollocation(4:6,:,i);
    vNorm     = sqrt( stateCollocation(4,:,i).^2 + stateCollocation(5,:,i).^2 +stateCollocation(6,:,i).^2 );
    m         = stateCollocation(7,:,i);
    rho       = mission.environment.rhoFun(h);
    D         = - 0.5 .* rho .*  v .* vNorm .* A .*cD;
    G         = - GM .* r ./rNorm.^3;
    acc       = (Thrust + D) ./ m + G;
    accMaxStage(i) = max(acc);
end
accMax = max(accMaxStage);

cin = [norm(stateCollocation(1:3,end,end) - mission.target)-1]; % [accMax-10*mission.environment.g0]
ceq = [];

end