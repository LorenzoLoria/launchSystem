clear all
close all
clc

mission = dataStruct ; 


% Initial condition for the simulation ( convention [r, theta, phi]
nOptPoints = 5 ;
rEarth = 6371 * 1e3 ;
x0 = [rEarth , 0, -pi/4] ; 
xf = [rEarth , pi, -pi/4] ;
target = xf ; 

% Initial Guess definition
initialGuessR = linspace(x0(1), x0(1) + 100 * 1e3, ceil(nOptPoints/2) + 1) ; 
initialGuessR = [initialGuessR, linspace( x0(1) + 100 * 1e3, xf(1), floor(nOptPoints/2) + 1) ]; 
initialGuessR = initialGuessR(2:end-1) ; 

initialGuessPhi = linspace(x0(3), pi/2, ceil(nOptPoints/2) + 1) ;
initialGuessPhi = [initialGuessPhi, linspace( pi/2, xf(3), floor(nOptPoints/2) + 1) ]; 
initialGuessPhi = initialGuessPhi(2:end-1) ; 

initialGuessV = linspace(1, 1e2, nOptPoints) ; 

initialGuess = [initialGuessR, initialGuessPhi , initialGuessV] ; 
initialGuess(6) = initialGuess(6) + 1e-6 ;
initialGuess(8) = initialGuess(8) - 1e-6 ;
initialGuess(9) = initialGuess(9) - 1e-6 ;
initialGuess(11) = initialGuess(11) + 1e-6 ;



thetaCoord = linspace(x0(2), xf(2), nOptPoints+2) ;

% Limits fro lower and upper boundary condition
heightLb = rEarth;
heightUb = rEarth + 200 * 1e3;
lowerBounds = [heightLb*ones(1,nOptPoints), -pi/2*ones(1,nOptPoints), 1*ones(1,nOptPoints)] ;
upperBounds = [heightUb*ones(1,nOptPoints), pi/2*ones(1,nOptPoints), 3e3*ones(1,nOptPoints)] ;

fminconOptions = optimoptions('fmincon', 'Display', 'iter-detailed', ...
        'MaxFunctionEvaluations', 1e5, ...
        'MaxIterations', 1000, ...
        'ConstraintTolerance', 1e-8, ...  
        'StepTolerance', 1e-6, ...
        'OptimalityTolerance', 1e-3,...   
        'FiniteDifferenceType','forward',...
        'DiffMinChange', 1e-5);  

a= 1 ;
optTraj = fmincon (@(x) fitnessFun(x, thetaCoord,x0, target, mission), initialGuess, ...
    [], [], [],[], lowerBounds, upperBounds, @(x) constraintFun(x,thetaCoord,x0,xf,mission), fminconOptions) ;



%% PLOT
rCoord = [x0(1) optTraj(1:nOptPoints) xf(1)];
phiCoord = [x0(3) optTraj(nOptPoints+1 : nOptPoints*2) xf(3)];
velocity =[10 optTraj(2*nOptPoints+1 : end) optTraj(end)];

rPoints       = [thetaCoord ; rCoord] ;
phiPoints     = [thetaCoord ; phiCoord] ;
vModulePoints = [thetaCoord ; vModule] ;

thetaQuery = linspace(x0(2), xf(2), 1000)';
nPointsPerInterval = 10 ;

[out] = splineGenerationEfficient(rPoints, nPointsPerInterval) ;
rQuery = out.profile ;

[out] = splineGenerationEfficient(phiPoints , nPointsPerInterval) ;
phiQuery = out.profile ;

[out] = splineGenerationEfficient(vModulePoints, nPointsPerInterval) ;
velocityProfile = out.profile ;



figure(1)
plot(thetaCoord, rCoord)

figure(2)
plot(xCoord,phiCoord)

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

