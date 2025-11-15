
clearvars
clc
close all

[mission] = dataStruct;
tSpan = [0 3*24*3600];
x0=[0;0;mission.envirnoment.rEarth+100e3;8000;0;0];
[tt,xx] = ballisticTrajectory(x0,mission);
plot3(xx(1,:),xx(2,:),xx(3,:),"r", "LineWidth",1)
hold on

EarthPlot(mission.envirnoment.rEarth)

function [dxdt] = keplerian_rhs(~, x,mission)

%   KEPLERIAN_RHS  Evaluates the right-hand-side of a 2-body (keplerian) propagator
%   Evaluates the right-hand-side of a newtonian 2-body propagator.
%   x is the state


GM = mission.environment.GM;

rho = interp1(mission.envirnoment.altRange,mission.envirnoment.rho,norm(x(1:3))-mission.envirnoment.rEarth);

v = norm(x(4:6));

D = - 0.5 * rho * v^2 * mission.capsule.Area * mission.capsule.Cd .* x(4:6)./v;

% Initialize right-hand-side
dxdt = zeros(6,1);

% Extract positions
rr = x(1:3);

% Compute square distance and distance
dist2 = dot(rr, rr);
dist = sqrt(dist2);

% Position detivative is object's velocity
dxdt(1:3) = x(4:6);   

% Compute the gravitational acceleration using Newton's law + air drag
dxdt(4:6) = - GM * rr /(dist*dist2) + D./mission.capsule.weigth;

end


function [tt,xx] = ballisticTrajectory(x0,mission)

tSpan = [0 3*24*3600];

options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@groundEvent);

solution = ode45(@(t,x) keplerian_rhs(t,x,mission), tSpan, x0,options);

tt = linspace(solution.x(1), solution.x(end), 1000);

xx = deval(solution,tt);

end

function [value,isterminal,direction] = groundEvent(t,x)
    h = norm(x(1:3))-6371e3;             
    value      = h; 
    isterminal = 1;       
    direction  = -1;      
end


function EarthPlot(r)
    % EarthPlot - plot of the Earth sphere with specific radius and texture
    % Input:
    %   r - Earth radius [km]

    textureImage = imread('EarthTexture.jpg');

    % Genera la superficie sferica
    [theta, phi] = meshgrid(linspace(0, 2*pi, size(textureImage, 2)), linspace(0, pi, size(textureImage, 1)));
    X = r * sin(phi) .* cos(theta);
    Y = r * sin(phi) .* sin(theta);
    Z = r * cos(phi);
    
    surface(X, Y, Z, 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'CData', textureImage);
    axis equal;
    grid on
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
   
end

