function[cin, ceq] = nlconGA(x,mission)

thrustDataVec = x;
GM  = mission.environment.GM;
A   = mission.capsule.Area;
Cd  = mission.capsule.supersonicCD;
g0  = mission.environment.g0;
tSpan = [0 mission.launcher.engines{1}.isp];
tVec = linspace(tSpan(1),tSpan(end),size(thrustDataVec,1));
%thrustData = @(t) [interp1(tVec,thrustDataVec(:,1),t);interp1(tVec,thrustDataVec(:,2),t);interp1(tVec,thrustDataVec(:,3),t)];
F_thrust = griddedInterpolant(tVec, thrustDataVec, 'linear', 'none');  
thrustData = @(t) F_thrust(t).';  

x0 = [mission.initialPoint'; 0; 0; 0; mission.launcher.engines{1}.m0];

[ttL,xxL] = launcherTrajectory(x0,mission, thrustData,0,thrustDataVec);

Tmax = 4 * mission.launcher.engines{1}.thrust;

for i = 1:length(ttL)
    r(i) = norm(xxL(1:3,i));
    h = r(i) - mission.environment.rEarth;
    rho = mission.environment.gridInterp(h);
    v = xxL(4:6,i);
        optVar = thrustData(ttL(i)); % thrustdata dovrebbe essere una funzione vettoriale con Tx,Ty,Tz
    
        theta = atan2( x(3),sqrt(x(1)^2 + x(2)^2) );  
    phi   = atan2( x(2), x(1) );

    Rz = [ cos(phi)  -sin(phi)  0;
       sin(phi)   cos(phi)  0;
       0          0         1 ];

    Ry = [ cos(theta)  0  -sin(theta);
       0           1  0;
      sin(theta)  0  cos(theta) ];

    R = Rz * Ry;
    percVec = optVar(1);
    thetaGimball = optVar(2);
    gammaGimball = 0;

    ThrustBRF = percVec * 4 * mission.launcher.engines{1}.thrust*[cos(thetaGimball)*cos(gammaGimball); cos(thetaGimball)*sin(gammaGimball); sin(thetaGimball)];


    ThrustIRF = R*ThrustBRF;
    D = - 0.5 .* rho .* norm(v) .* A .* Cd .* v;
    G = - GM * r(i) /norm(r(i))^3;
    m(i) = xxL(7,i);
    acc(i) = norm((ThrustBRF + D)/m(i) + G);
    
end


x0C = xxL(1:6,end);

windDirection = -1+2*rand(3,1);
windDirection = windDirection/norm(windDirection);

[ttC,xxC] = ballisticTrajectory(x0C,mission,windDirection,ttL(end));

cin = [(norm(xxC(1:3,end) - mission.target))-0.6];

ceq =[];% [-min(acc);max(acc)-15*g0];
end