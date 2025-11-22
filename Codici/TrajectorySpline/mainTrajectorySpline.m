clear all
clc
close all

% Initial point and Target Point (da Modificare)
nOptPoints = 5 ;
x0 = [0 0 0] ; 
x0 = [0 0 120000] ;

x_target = [12 15 8] ;
x_target = [300000 1000 5000] ;

% Extrapolatin of Mission Data from Struct
mission = dataStruct;

% Random Initial Guess
initialGuessY = linspace(x0(2), x_target(2), nOptPoints+2) ; 
initialGuessY = initialGuessY(2:end-1) ; 
initialGuessZ = linspace(x0(3), x_target(3), nOptPoints+2) ; 
initialGuessZ = initialGuessZ(2:end-1) ; 
initialGuessV = linspace(1, 5, nOptPoints) ; 
initialGuessV = linspace(2500, 300, nOptPoints);

initialGuess = [initialGuessY, initialGuessZ , initialGuessV] ; 
% initialGuess = [linspace(x0(2), xf(2), nOptPoints), linspace(x0(3), xf(3), nOptPoints) , linspace(5, 10, nOptPoints)] ; 

xCoord = linspace(x0(1), x_target(1), nOptPoints+2) ;

% Lower and Upper Boundary Definition
limit = 20 ;
lowerBounds = [-50000*ones(nOptPoints,1), 0*ones(nOptPoints,1), 100*ones(nOptPoints,1)] ;
upperBounds = [50000*ones(nOptPoints,1), 150000*ones(nOptPoints,1), 4000*ones(nOptPoints,1)] ;

fminconOptions = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter-detailed', ...
        'MaxFunctionEvaluations', 1e5, ...
        'MaxIterations', 1000, ...
        'ConstraintTolerance', 1e-8, ...  
        'StepTolerance', 1e-6, ...
        'OptimalityTolerance', 1e-3,...   
        'FiniteDifferenceType','forward',...
        'DiffMinChange', 1e-5);    


optTraj = fmincon (@(x) fitnessFun(x, xCoord,x0,x_target,mission), initialGuess, ...
    [], [], [],[], lowerBounds, upperBounds, @(x) constraintFun(x,xCoord,x0,x_target, mission), fminconOptions) ;


%% PLOT
% Reconstruct the full coordinate vectors including start/end points
yCoord = [x0(2) optTraj(1:nOptPoints) x_target(2)];
zCoord = [x0(3) optTraj(nOptPoints+1 : nOptPoints*2) x_target(3)];
velocity =[2500 optTraj(2*nOptPoints+1 : end) optTraj(end)];

% Create point sets for spline generation
yPoints = [xCoord ; yCoord]';
zPoints = [xCoord ; zCoord]';
velocityPoints = [xCoord ; velocity]';

% Profile Generation
xQuery = linspace(x0(1), x_target(1), 1000)';

[out] = splineGeneration(yPoints, xQuery) ;
yQuery = out.profile ;

[out] = splineGeneration(zPoints, xQuery) ;
zQuery = out.profile ;

[out] = splineGeneration(velocityPoints, xQuery) ;
velocityProfile = out.profile ;

% Reconstruct the initial guess path for comparison
initialGuessYFull = [x0(2) initialGuess(1:length(initialGuess)/3) x_target(2) ];
initialGuessZFull = [x0(3) initialGuess(length(initialGuess)/3 +1 :length(initialGuess)*2/3) x_target(3) ];

figure(1)
plot(xCoord, yCoord)

figure(2)
plot(xCoord,zCoord)

figure(3)
plot(xCoord,velocity)

figure(4)
plot3(xCoord,yCoord,zCoord)



figure(5)
hold on
plot3(xQuery, yQuery(:,2), zQuery(:,2), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Optimized Path');
plot3(xCoord, initialGuessYFull,  initialGuessZFull, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Initial Guess');
[xSphere, ySphere, zSphere] = sphere;
radius = 5;
surf(radius * xSphere + 12, radius * ySphere + 9.5, radius * zSphere + 3.4, 'EdgeAlpha', 0.1);
title('3D Optimized Trajectory vs. Initial Guess');
xlabel('X Coordinate');
ylabel('Y Coordinate');
zlabel('Z Coordinate');
legend('show', 'Location', 'best');
axis equal;
grid on;



