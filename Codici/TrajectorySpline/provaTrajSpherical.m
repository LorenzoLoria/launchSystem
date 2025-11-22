clear all
close all
clc

mission = dataStruct ; 


% Initial condition for the simulation ( convention [r, theta, phi]
nOptPoints = 5 ;
rEarth = 6371 * 1e3 ;
x0 = [rEarth , 0, pi/4] ; 
xf = [rEarth + 100e3 , pi, pi/4] ;
target = xf ; 

% Initial Guess definition
initialGuessR = linspace(x0(1), x0(1) + 100 * 1e3, ceil(nOptPoints/2) + 1) ; 
initialGuessR = [initialGuessR, linspace( x0(1) + 100 * 1e3, xf(1), floor(nOptPoints/2) + 1) ]; 
initialGuessR = initialGuessR(2:end-1) ; 

initialGuessPhi = linspace(x0(3), pi/2, ceil(nOptPoints/2) + 1) ;
initialGuessPhi = [initialGuessPhi, linspace( pi/2, xf(3), floor(nOptPoints/2) + 1) ]; 
initialGuessPhi = initialGuessPhi(2:end-1) ; 

initialGuessV = linspace(1e2, 3e2, nOptPoints) ; 

initialGuess = [initialGuessR, initialGuessPhi , initialGuessV] ; 


thetaCoord = linspace(x0(2), xf(2), nOptPoints+2) ;


figure(1)

x0Plot = [ x0(1) .* sin(x0(3)) .* cos(x0(2)) ;   %X
        x0(1) .* sin(x0(3)) .* sin(x0(2)) ;   %Y
        x0(1) .* cos(x0(3))];  
xfPlot = [ xf(1) .* sin(xf(3)) .* cos(xf(2)) ;   %X
        xf(1) .* sin(xf(3)) .* sin(xf(2)) ;   %Y
        xf(1) .* cos(xf(3))];  

plot3(x0Plot(1),x0Plot(2),x0Plot(3), 'xr', 'LineWidth',3)
hold on
axis equal
plot3(xfPlot(1),xfPlot(2),xfPlot(3), 'xr','LineWidth',3)
[xSphere,ySphere,zSphere] = sphere ;
radius = rEarth ; 
surf(radius * xSphere  , radius * ySphere ,radius * zSphere , 'EdgeAlpha',0.3 )

inGuess = [initialGuessR ; thetaCoord(2:end-1) ; initialGuessPhi] ;
for ii = 1:length(initialGuessR)
    pts(:,ii) = [ inGuess(1,ii) .* sin(inGuess(3,ii)) .* cos(inGuess(2,ii)) ;   %X
        inGuess(1,ii) .* sin(inGuess(3,ii)) .* sin(inGuess(2,ii)) ;   %Y
        inGuess(1,ii) .* cos(inGuess(3,ii))]; 
end

plot3(pts(1,:), pts(2,:), pts(3,:), 'ob', 'LineWidth',3)



%%




% Limits fro lower and upper boundary condition
heightLb = rEarth;
heightUb = rEarth + 200 * 1e3;
lowerBounds = [heightLb*ones(1,nOptPoints), 0*ones(1,nOptPoints), 300*ones(1,nOptPoints)] ;
upperBounds = [heightUb*ones(1,nOptPoints), pi*ones(1,nOptPoints), 3e3*ones(1,nOptPoints)] ;

fminconOptions = optimoptions('fmincon', 'Display', 'iter-detailed', ...
        'MaxFunctionEvaluations', 1e5, ...
        'MaxIterations', 1000, ...
        'ConstraintTolerance', 1e-12, ...  
        'StepTolerance', 1e-10, ...
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
vModulePoints = [thetaCoord ; velocity] ;

thetaQuery = linspace(x0(2), xf(2), 1000)';
nPointsPerInterval = 100 ;

[out] = splineGenerationEfficient(rPoints, nPointsPerInterval) ;
rQuery = out.profile' ;

[out] = splineGenerationEfficient(phiPoints , nPointsPerInterval) ;
phiQuery = out.profile' ;

[out] = splineGenerationEfficient(vModulePoints, nPointsPerInterval) ;
velocityProfile = out.profile' ;

%%

[output1] = trajectoryGenerationSpherical(optTraj, thetaCoord, x0, xf, mission, 1) ;


[output2] = trajectoryGenerationSpherical(optTraj, thetaCoord, x0, xf, mission, 0) ;


velocityPlot = figure(5) ;
plot(output1.theta, output1.vel)
title("Velocity spline")


accelerationPlot = figure(6) ;
plot(output1.theta, output2.acceleration)
title("Acceleration spline")


thrustPlot = figure(7) ;
for ii = 1:length(output2.thrust(1,:))

thrustNorm(ii) =  sqrt(sum(output2.thrust(:,ii).^2)) ;

end

plot(output2.theta, thrustNorm)
title("Thrust")



%%
figure(2)
plot(thetaCoord, rCoord)
title("Radius")


figure(3)
plot(thetaCoord,phiCoord)
title("Phi")


figure(4)
plot(thetaCoord,velocity)
title("Velocity")


% figure(5)
% plot3(thetaCoord,yCoord,zCoord)
% 
% figure(6)
% plot3(xQuery, yQuery(:,2), zQuery(:,2))
% hold on
% plot3(xCoord, [x0(2) initialGuess(1:length(initialGuess)/3) xf(2) ], [x0(3) initialGuess(length(initialGuess)/3 +1 :length(initialGuess)*2/3) xf(3) ])
% [xSphere,ySphere,zSphere] = sphere ;
% radius = 5 ; 
% surf(radius * xSphere + 12 , radius * ySphere + 9.5 ,radius * zSphere + 3.4, 'EdgeAlpha',0.3 )
% axis equal


figure(1)
for ii = 1 : length(rQuery)

pts(:,ii) = [ rQuery(2,ii) .* sin(phiQuery(2,ii)) .* cos(rQuery(1,ii)) ;   %X
        rQuery(2,ii) .* sin(phiQuery(2,ii)) .* sin(rQuery(1,ii)) ;   %Y
        rQuery(2,ii) .* cos(phiQuery(2,ii))]; 


end

%%

plot3(pts(1,:), pts(2,:),pts(3,:), 'g')






