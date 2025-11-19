clear
clc
close all


% Initial condition for the simulation ( convention [r, theta, phi]) 
% --> DA METTERE IN RADINTI
nOptPoints = 5 ;
rEarth = 6371 * 1e3 ;
x0 = [sqrt(2) * rEarth , 0, -45] ; 
xf = [sqrt(2) * rEarth , 180, -45] ;
target = xf ; 

% Initial Guess definition
initialGuessR = linspace(x0(1), x0(1) + 100 * 1e3, ceil(nOptPoints/2) + 1) ; 
initialGuessR = [initialGuessR, linspace( x0(1) + 100 * 1e3, x0(1), floor(nOptPoints/2) + 1) ]; 
initialGuessR = initialGuessY(2:end-1) ; 

initialGuessTheta = linspace(x0(2), xf(2), nOptPoints+2) ; 
initialGuessTheta = initialGuessTheta(2:end-1) ; 

initialGuessV = linspace(1, 1e2, nOptPoints) ; 
initialGuess = [initialGuessR, initialGuessTheta , initialGuessV] ; 


phiCoord = linspace(x0(3), 0, ceil(nOptPoints/2) + 1) ;
phiCoord = [phiCoord, fliplr(linspace( 0, x0(3),floor(nOptPoints/2) + 1) )];

% Limits fro lower and upper boundary condition
heightLb = rEarth;
heightUb = rEarth + 200 * 1e3;
lowerBounds = [heightLb*ones(nOptPoints,1), 0*ones(nOptPoints,1), 1*ones(nOptPoints,1)] ;
upperBounds = [heightUb*ones(nOptPoints,1), 360*limit*ones(nOptPoints,1), 100*ones(nOptPoints,1)] ;

fminconOptions = optimoptions('fmincon', 'Display', 'iter-detailed', ...
        'MaxFunctionEvaluations', 1e5, ...
        'MaxIterations', 1000, ...
        'ConstraintTolerance', 1e-8, ...  
        'StepTolerance', 1e-6, ...
        'OptimalityTolerance', 1e-3,...   
        'FiniteDifferenceType','forward',...
        'DiffMinChange', 1e-5);  

optTraj = fmincon (@(x) fitnessFun(x, xCoord,target,x0,xf), initialGuess, ...
    [], [], [],[], lowerBounds, upperBounds, @(phi) constraintFun(phi,phiCoord,target,x0,xf), fminconOptions) ;
