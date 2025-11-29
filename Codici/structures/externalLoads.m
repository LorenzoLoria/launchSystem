function [q,an,at,T,D,angle,gamma] = externalLoads(timeCollocation,stateCollocation,mission,opt,thrustData)

vel = stateCollocation(4:6,:,1:end-1)-stateCollocation(4:6,1,1);
absVel = sqrt ( vel(1,:).^2+ vel(2,:).^2 + vel(3,:).^2 );

vel = vel(1:3,:);

pos = stateCollocation(1:3,:,1:end-1);

absH = sqrt ( pos(1,:,:).^2+ pos(2,:,:).^2 + pos(3,:,:).^2 )-mission.environment.rEarth;
absH = absH(:);

rhoVec = mission.environment.rhoFun(absH);

q = 0.5*rhoVec'.*(absVel).^2;

timeStage = timeCollocation(:,1:end-1);
timeStage = timeStage(:);
 
aCC = [diff(vel(1,:)) ; diff(vel(2,:));diff(vel(3,:))]./diff(timeStage)';

normRocket = vel./absVel ; 

an = dot(aCC,normRocket(:,1:end-1));

aTan = cross(aCC,normRocket(:,1:end-1));
at = sqrt ( aTan(1,:).^2+ aTan(2,:).^2 + aTan(3,:).^2 );


ThrustTotal = zeros(3,size(vel,2));
for i = 1:opt.nStages

    t0 = timeCollocation(1,i);
    tf = timeCollocation(end,i);
    tSpan = [t0 t0+tf]; %da rivedere nel caso i tempi non vadano bene
    
    tVec = linspace(tSpan(1),tSpan(end),size(thrustData,1));    
    fThrust = griddedInterpolant(tVec, thrustData(:,:,i), 'linear', 'none'); 
    thrustFunction= @(t) fThrust(t).';

    optVar = thrustFunction(timeCollocation(:,i));
    absThrust = optVar(1,:) .* opt.stage{i}.nEngines .*opt.stage{i}.engine.thrust;
    ThrustBRF = [absThrust;absThrust;absThrust].* [cos( deg2rad (optVar(3,:) )).*cos(deg2rad (optVar(2,:))); cos(deg2rad (optVar(3,:))).*sin(deg2rad (optVar(2,:))); sin(deg2rad (optVar(3,:)))];
    ThrustIRF = mission.Rfinal'*ThrustBRF;
    ThrustTotal(:,(i-1)*100+1: i*100) = ThrustIRF;
end    

absThrust2 = sqrt ( ThrustTotal(1,:).^2+ ThrustTotal(2,:).^2 + ThrustTotal(3,:).^2 );
angle = acos(dot(vel,ThrustTotal)./(absThrust2.*absVel) );
T = absThrust2;

[qmax,index] = max(q);

an = an(index);
at = at(index);
T = T(index);
angle = angle(index);
rho = rhoVec(index);
vQmax = absVel(index);

D = q * mission.capsule.supersonicCD* mission.capsule.Area;

D = D(index);

distmaxq = pos(:,:,1);

rVectorBoh = distmaxq(:,index);

gamma = acos ( dot(-rVectorBoh,vel(:,index))/(norm(rVectorBoh)*norm(vel(:,index))) )-pi/2 ;

end