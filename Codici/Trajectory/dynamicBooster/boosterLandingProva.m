clear all
close all 
clc


r0 = [0 0 100] ;
v0 = [0 0 0] ;
q0 = [1 0 0 0] ;
w0 = [0 0 0] ;

y0 = [r0 v0 q0 w0]' ;
tSpan = linspace(0,200,501);


targetPos = [10 0 0]' ;
targetVel = [0 0 -0.3]' ;

kPos = [50 0 10]' ;
kVel = [1 0 60]' ;
KpAtt = 100 ;
KdAtt = 1000 ;


gains = [kPos ; kVel ; KpAtt ; KdAtt] ;
gains = [528 1 464 86 429 61 824 349]' ; 
%gains = [528 1 400 86 429 61 824 349]' ; 




options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@touchDown);
[t, sol] = ode113(@(t,x) dynamicsAndControl(t,x,targetPos,targetVel,gains,r0), tSpan, y0,options);


%% PLOT 3D

if r0(2)~=0

close all

figure(1)
plot3(sol(:,1), sol(:,2), sol(:,3))
hold on 
axis equal
plot3(r0(1),r0(2),r0(3), 'xr', 'LineWidth', 3)
plot3(targetPos(1),targetPos(2),targetPos(3), 'xr', 'LineWidth', 3)
xlabel("x axis")
ylabel("y axis")

q = sol(:,7:10) ;
boosterAxis = [2.*(q(:,2).*q(:,4) + q(:,3).*q(:,1)) , 2.*(q(:,3).*q(:,4) - q(:,2).*q(:,1)) , q(:,1).^2 - q(:,2).^2 - q(:,3).^2 + q(:,4).^2] ;

arrowScale = 1 ;
for ii=1:length(sol(:,1))
quiver3(sol(ii,1), sol(ii,2), sol(ii,3), arrowScale*boosterAxis(ii,1), arrowScale*boosterAxis(ii,2), arrowScale*boosterAxis(ii,3),'Color', 'g', 'AutoScale','off')

end
figure(2)
plot(t,sol(:,6))
title("vertical velocity")

figure(3)
plot(t,sol(:,4))
title("x velocity")

figure(4)
plot(t,sol(:,5))
title("y velocity")

else
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

end




