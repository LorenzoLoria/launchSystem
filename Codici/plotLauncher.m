%% Plot rocket
close all


hFuel  = sol1(2);
hOx  = sol1(3);
r1 = sol1(1);

thetaLow = linspace(pi,2*pi,1000);
thetaUpp = linspace(0,pi,1000);

hFuelVec = linspace(0,hFuel,2);
hOxVec = linspace(0,hOx,2);

h0 = 3;
figure

plot(r1*cos(thetaLow),h0 +r1*sin(thetaLow), 'r')
hold on
plot([r1,r1],[hFuelVec(1)+h0,hFuelVec(2)+h0],'r')
plot([-r1,-r1],[hFuelVec(1)+h0,hFuelVec(2)+h0],'r')
plot(r1*cos(thetaUpp),h0 +r1*sin(thetaUpp)+hFuel, 'r')

plot(r1*cos(thetaLow),h0 +r1*sin(thetaLow)+hFuel+r1+r1, 'b')
hold on
plot([r1,r1],[hOxVec(1)+h0+r1+hFuel+r1,hOxVec(2)+h0+hFuel+r1+r1],'b')
plot([-r1,-r1],[hOxVec(1)+h0+r1+hFuel+r1,hOxVec(2)+h0+hFuel+r1+r1],'b')
plot(r1*cos(thetaUpp),h0 +r1*sin(thetaUpp)+hFuel+r1+r1+hOx, 'b')
axis equal


hFuel2  = sol2(2);cd
hOx2  = sol2(3);
r2 = sol2(1);

hFuelVec2 = linspace(0,hFuel2,2);
hOxVec2 = linspace(0,hOx2,2);

h02 = hOxVec(2)+h0+hFuel+r1+r1+r1+r2+0.5;

thetaLow2 = linspace(pi,2*pi,1000);
thetaUpp2 = linspace(0,pi,1000);

plot(r2*cos(thetaLow),h02 +r2*sin(thetaLow), 'r')
hold on
plot([r2,r2],[hFuelVec2(1)+h02,hFuelVec2(2)+h02],'r')
plot([-r2,-r2],[hFuelVec2(1)+h02,hFuelVec2(2)+h02],'r')
plot(r2*cos(thetaUpp),h02 +r2*sin(thetaUpp)+hFuel2, 'r')

plot(r2*cos(thetaLow),h02 +r2*sin(thetaLow)+hFuel2+r2+r2, 'b')
hold on
plot([r2,r2],[hOxVec2(1)+h02+r2+hFuel2+r2,hOxVec2(2)+h02+hFuel2+r2+r2],'b')
plot([-r2,-r2],[hOxVec2(1)+h02+r2+hFuel2+r2,hOxVec2(2)+h02+hFuel2+r2+r2],'b')
plot(r2*cos(thetaUpp),h02 +r2*sin(thetaUpp)+hFuel2+r2+r2+hOx2, 'b')
axis equal

