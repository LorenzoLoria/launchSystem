clear all
close all 
clc

warning ('off','all')



lowerBounds = 1 * ones(8,1) ; 
upperBounds = 1000 * ones(8,1) ;
nVars = length(lowerBounds) ;
lowerBounds = [10 1 10 10 10 10 100 100]'; 
upperBounds = [1000 1 1000 1000 1000 1000 1000 1000]'; 


r0 = [0 5 50] ;
v0 = [0 0 0] ;
q0 = [1 0 0 0] ;
w0 = [0 0 0] ;
y0 = [r0 v0 q0 w0]' ;
tSpan = linspace(0,20,51);


targetPos = [10 0 0]' ;
targetVel = [0.1 0 -0.1]' ;



gains0 = [50 0 10 1 0 60 100 1000]' ;
intCon = 1:8 ;

% options_ga = optimoptions("ga", ...
%     "Display","iter", ...
%     "MaxGenerations",30, ...
%     "PopulationSize",20,...
%     "UseParallel",true,...
%     "FunctionTolerance", 1e-4,...
%     'EliteCount',  6,...
%     "CrossoverFraction", 0.9);
options_ga = optimoptions("ga", ...
    "Display","iter", ...
    "MaxGenerations",30, ...
    "PopulationSize",200,...
    "UseParallel",true,...
    "FunctionTolerance", 1e-4,...
    "MaxStallGenerations", 10,...
    'EliteCount',  6);




[gains] = ga(@(x) findGains(x,targetPos,targetVel,y0,tSpan,r0), nVars, [],[],[],[],lowerBounds,upperBounds, [] ,intCon, options_ga) ;







options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@touchDown);
[t, sol] = ode113(@(t,x) dynamicsAndControl(t,x,targetPos,targetVel,gains,r0), tSpan, y0,options);


%% PLOT 2D
close all

figure(1)
plot(sol(:,1), sol(:,3))
hold on 
axis equal
plot(r0(1),r0(3), 'xr', 'LineWidth', 3)
plot(targetPos(1),targetPos(3), 'xr', 'LineWidth', 3)


q = sol(:,7:10) ;
boosterAxis = [2.*(q(:,2).*q(:,4) + q(:,3).*q(:,1)) , 2.*(q(:,3).*q(:,4) - q(:,2).*q(:,1)) , q(:,1).^2 - q(:,2).^2 - q(:,3).^2 + q(:,4).^2] ;

arrowScale = 1 ;
for ii=1:length(sol(:,1))
quiver(sol(ii,1), sol(ii,3), arrowScale*boosterAxis(ii,1), arrowScale*boosterAxis(ii,3),'Color', 'g', 'AutoScale','off')

end

figure(2)
plot(t,sol(:,6))
title("vertical velocity")

figure(3)
plot(t,sol(:,4))
title("horizontal velocity")

