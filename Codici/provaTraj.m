clear
clc
close all

nOptPoints = 5 ;
x0 = [0 0 0] ; 
xf = [12 15 8] ;
target = xf ; 

initialGuess = [linspace(x0(2), xf(2), nOptPoints), linspace(x0(3), xf(3), nOptPoints), linspace(1, 30, nOptPoints)] ;

xCoord = linspace(x0(1), xf(1), nOptPoints+2) ;

limit = 20 ;
lowerBounds = [-limit*ones(nOptPoints,1), -limit*ones(nOptPoints,1), 1*ones(nOptPoints,1)] ;
upperBounds = [limit*ones(nOptPoints,1), limit*ones(nOptPoints,1), 100*ones(nOptPoints,1)] ;

fminconOptions = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter-detailed', ...
        'MaxFunctionEvaluations', 1e5, ...
        'MaxIterations', 1e5, ...
        'ConstraintTolerance', 1e-12, ...
        'StepTolerance', 1e-9, ...
        'OptimalityTolerance', 1e-10,...
        'FiniteDifferenceType','central');

gaOptions = optimoptions('ga', 'PlotFcn',{'gaplotbestf', 'gaplotbestindiv'}, 'display', 'iter','MaxStallGenerations', 10, ...
        'FunctionTolerance', 1e-6, 'EliteCount',  6,...
        'MaxGenerations', 1e4, 'PopulationSize', 40, ...
        'NonlinearConstraintAlgorithm', 'penalty',...
        'SelectionFcn', {@selectiontournament, 4},...
        'UseParallel',true,...
        'MutationFcn', {@mutationpower},...
        'CrossoverFcn', {@crossoverintermediate, 0.8});


optTraj = fmincon (@(x) fitnessFun(x, xCoord,target,x0,xf), initialGuess, ...
    [], [], [],[], lowerBounds, upperBounds, @(x) constraintFun(x,xCoord,target,x0,xf), fminconOptions) ;

% nVars = nOptPoints * 3 ;
% [optTraj, variabileInutile] = ga(  @(x) fitnessFun(x, xCoord,target,x0,xf) ...
%     ,nVars, [],[],[],[], lowerBounds, upperBounds ...
%     , @(x) constraintFun(x,xCoord,target,x0,xf), [] , gaOptions);
%% PLOT
yCoord = [x0(2) optTraj(1:nOptPoints) xf(2)];
yPoints = [xCoord ; yCoord]';
zCoord = [x0(3) optTraj(nOptPoints+1 : nOptPoints*2) xf(3)];
zPoints = [xCoord ; zCoord]';
velocity =[0.1 optTraj(2*nOptPoints+1 : end) optTraj(end)];
velocityPoints = [xCoord ; velocity]';

xQuery = linspace(x0(1), xf(1), 1000)';

[out] = splineGeneration(yPoints, xQuery) ;
yQuery = out.profile ;

[out] = splineGeneration(zPoints, xQuery) ;
zQuery = out.profile ;

[out] = splineGeneration(velocityPoints, xQuery) ;
velocityProfile = out.profile ;



figure(1)
plot(xCoord, yCoord)

figure(2)
plot(xCoord,zCoord)

figure(3)
plot(xCoord,velocity)

figure(4)
plot3(xCoord,yCoord,zCoord)

figure(5)
plot3(xQuery, yQuery(:,2), zQuery(:,2))
hold on
plot3(xCoord, [x0(2) initialGuess(1:length(initialGuess)/3) xf(2) ], [x0(3) initialGuess(length(initialGuess)/3 +1 :length(initialGuess)*2/3) xf(3) ])
[xSphere,ySphere,zSphere] = sphere ;
radius = 5 ; 
surf(radius * xSphere + 12 , radius * ySphere + 9.5 ,radius * zSphere + 3.4, 'EdgeAlpha',0.3 )
axis equal


%% 

pos = [xQuery yQuery(:,2) zQuery(:,2)]' ;
earthCenter = [12 ; 9.5 ; 3.4];
earthRadius = 5 ; 
distance = zeros(1,length(pos(1,:))) ;

for ii = 1:length(pos(1,:))
    distance(1,ii) = sqrt(sum((pos(:,ii) - earthCenter).^2)) ;

end







function [objective] = fitnessFun(x,xCoord,target,x0,xf)

fitnessFlag = 1 ;

[output] = trajectoryGeneration(x,xCoord,target,x0,xf,fitnessFlag);


vel = output.vel ; 
yCoord = output.yCoord ; 
zCoord = output.zCoord ; 
s = output.s ;

tFinal = trapz(s,abs(1./vel));


yOscillation = sqrt(sum(diff(yCoord).^2));
zOscillation = sqrt(sum(diff(zCoord).^2));


%objective = 1/m(end) + tFinal + yOscillation / 15 + zOscillation / 8 ; 
objective = tFinal + yOscillation / 5 + zOscillation / 3 ; 


end


function [mDot] = massIntegration(s,m,vModule, dvModule, dPosition)

rho = 1.225 ; 
surface = 1 ; 
CD = 0.2 ;
iSp = 300 ; 
g0 = 9.81 ;


direction = dPosition(s) ./ norm(dPosition(s)) ; 
thrust = m * direction * dvModule(s) - m * [0 ; 0 ; g0] - 0.5 * rho * (vModule(s)).^2 * surface * CD .* direction ;
mDot = - 1 / iSp / g0 * (norm(thrust));

if mDot>0
    keyboard
end



end


function [distance,isterminal,direction] = targetEvent(s,m,position,target)


    distance = norm(position(s) - target') - 1 ;


    isterminal = 1;       
    direction  = -1;      
end


function [C,Ceq] = constraintFun(x,xCoord,target,x0,xf)

Ceq = [] ;
fitnessFlag = 0 ; 

[output] = trajectoryGeneration(x,xCoord,target,x0,xf,fitnessFlag);


earthCenter = [12 ; 9.5 ; 3.4];
tol = 1.5 ;
earthRadius = 5 + tol; 


acceleration = output.acceleration ; 
thrustOld = output.thrust ; 
position = output.position ; 
s = output.s ;  
directionOld = output.direction ;

evaluationPoints = linspace(s(1), s(end), 5e2) ; 

thrust(1,:) = interp1(s,thrustOld(1,:),evaluationPoints);
thrust(2,:) = interp1(s,thrustOld(2,:),evaluationPoints);
thrust(3,:) = interp1(s,thrustOld(3,:),evaluationPoints);

direction(1,:) = interp1(s,directionOld(1,:),evaluationPoints);
direction(2,:) = interp1(s,directionOld(2,:),evaluationPoints);
direction(3,:) = interp1(s,directionOld(3,:),evaluationPoints);

s = interp1(s,s,evaluationPoints) ; 





pos = zeros(3,length(s)) ; 
distance = zeros(length(s),1) ; 
acc = zeros(length(s),1) ; 
steeringAngle = zeros(length(s),1) ; 


for ii = 1:length(s)
    
    pos(:,ii) = position(s(ii));
    distance(ii,1) = sqrt(sum((pos(:,ii) - earthCenter).^2)) ;
    acc(ii,1) = acceleration(s(ii));
    steeringAngle(ii,1) = acosd(dot(thrust(ii) /norm(thrust(ii)), direction(ii))) ;

    
end

g0 = 9.81 ;
maxAcc = 5 * g0 ; 

maxSteeringAngle = 45 ; 

maxThrust = 2500 ; 


%C = [earthRadius - min(distance) ; max(acc) - maxAcc ; max(steeringAngle) - maxSteeringAngle ; max(norm(thrust)) - maxThrust];
C = [earthRadius - min(distance) ; max(acc) - maxAcc ; max(norm(thrust)) - maxThrust ];

% maxAngle = max(steeringAngle) ;


end


function [output] = trajectoryGeneration(x,xCoord,target,x0,xf,fitnessFlag)



yCoord = [x0(2) x(1:length(x)/3) xf(2)];
zCoord = [x0(3) x(length(x)/3+1 : length(x)/3*2) xf(3) ];
vModule = [0.1 x(length(x)/3*2+1 : end) x(end)];

yPoints = [xCoord ; yCoord]' ;
zPoints = [xCoord ; zCoord]' ;
vModulePoints = [xCoord ; vModule]' ;


[out] = splineGeneration(yPoints) ;

y = out.polynomial ; 
dy = out.firstDer ;
ddy = out.secondDer ;

[out] = splineGeneration(zPoints) ;

z = out.polynomial ; 
dz = out.firstDer ;
ddz = out.secondDer ;

[out] = splineGeneration(vModulePoints) ;

vModule = out.polynomial ; 
dvModule = out.firstDer ;
ddvModule = out.secondDer ;


position = @(s) [ s ; y(s) ; z(s) ];
dPosition = @(s) [1 ; dy(s) ; dz(s) ];
dPositionModule = @(s) norm(dPosition(s)) ;


iSp = 300 ; 
g0 = 9.81 ; 

sSpan = [0, 100] ; 
m0 = 100 ;

opt = odeset('relTol', 1e-12, 'absTol', 1e-12, 'Events' , @(s,m)targetEvent(s,m,position,target), 'MaxStep', 0.1) ;

[s, m]=ode45(@(s,m) massIntegration(s, m, vModule, dvModule, dPosition), sSpan, m0, opt) ; 

vel = zeros(1,length(s));
direction = zeros(3,length(s));
thrust = zeros(3,length(s));


rho = 1.225 ;
surface = 1 ;
CD = 0.2 ; 

for ii = 1:length(s)

    vel(ii) = vModule(s(ii)) ;
    direction(:,ii) = dPosition(s(ii)) ./ norm(dPosition(s(ii))) ;
    thrust(:,ii) = m(ii) * direction(:,ii) * dvModule(s(ii)) - m(ii) * [0 ; 0 ; g0] - 0.5 * rho * (vModule(s(ii))).^2 * surface * CD .* direction(:,ii) ;


end





if fitnessFlag 

    output.vel = vel ;
    output.yCoord = yCoord ;
    output.zCoord = zCoord ;
    output.s = s ; 

else

    output.acceleration = dvModule ; 
    output.position = position ; 
    output.thrust = thrust ; 
    output.s = s ; 
    output.direction = direction ;

end














end













