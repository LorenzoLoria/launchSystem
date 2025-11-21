function [output] = trajectoryGenerationSpherical(x, thetaCoord, x0, xf, mission, fitnessFlag)
    % TRAJECTORYGENERATIONSPHERICAL
    % Generates a trajectory defined in Spherical Coordinates (r, theta, phi)
    % and converts it to Cartesian (x, y, z) for physics integration.
    %
    % INPUTS:
    %   x           : Optimization vector containing control points for [r, phi, v]
    %   sCoord      : The independent variable grid (Theta/Polar Angle )
    %   target      : Final target state (not strictly used here if xf is used)
    %   x0          : Initial state [r0, theta0, phi0] (Spherical)
    %   xf          : Final state [rf, thetaf, phif] (Spherical)
    %   fitnessFlag : Boolean to select output format 
    %   mission     : Struct containing mission/environment parameters (ANCORA DA AGGIUNGERE)
    %
    % OUTPUTS:
    %   output      : Struct containing time histories (pos, vel, acc, thrust, etc.)


persistent prevSimulation  prevVel prevRCoord prevPhiCoord prevTheta prevAccelerationModule prevPosition prevThrust prevVelocityUnitVector prevTrajectoryLength

if isempty(prevSimulation)
    prevSimulation = zeros(1,length(x)) ; 
end

if x == prevSimulation

   vel = prevVel ; 
   rCoord = prevRCoord ; 
   phiCoord = prevPhiCoord ;  
   theta = prevTheta ;
   accelerationModule = prevAccelerationModule ;
   position = prevPosition ;
   thrust = prevThrust ;
   velocityUnitVector = prevVelocityUnitVector ; 
   trajectoryLength = prevTrajectoryLength ;
   
else



    surface = mission.capsule.Area;
    CD      = mission.capsule.supersonicCD; 
    GM      = mission.environment.GM;
    m0      = mission.launcher.engines{1}.m0 ; 


    
    rCoord   = [x0(1) x(1:length(x)/3) xf(1)];
    phiCoord = [x0(3) x(length(x)/3+1 : length(x)/3*2) xf(3) ];
    vModule  = [10 x(length(x)/3*2+1 : end) x(end)];
    
    % Prepare data for spline generation
    rPoints       = [thetaCoord ; rCoord] ;
    phiPoints     = [thetaCoord ; phiCoord] ;
    vModulePoints = [thetaCoord ; vModule] ;




    % Variables passed as functions, transitioning to only vectors
    %
    %
    % discretizationLength = 0.1 ;
    %
    % [rSpline] = splineGenerationEfficient(rPoints, discretizationLength) ;
    %
    % r = griddedInterpolant(rSpline.profile(:,1),rSpline.profile(:,2),"linear","linear");
    % dr = griddedInterpolant(rSpline.dProfile(:,1),rSpline.dProfile(:,2),"linear","linear");
    % ddr = griddedInterpolant(rSpline.ddProfile(:,1),rSpline.ddProfile(:,2),"linear","linear");
    %
    % [phiSpline] = splineGenerationEfficient(phiPoints, discretizationLength) ;
    %
    % phi = griddedInterpolant(phiSpline.profile(:,1),phiSpline.profile(:,2),"linear","linear");
    % dPhi = griddedInterpolant(phiSpline.dProfile(:,1),phiSpline.dProfile(:,2),"linear","linear");
    % ddPhi = griddedInterpolant(phiSpline.ddProfile(:,1),phiSpline.ddProfile(:,2),"linear","linear");
    %
    % [vSpline] = splineGenerationEfficient(vModulePoints, discretizationLength) ;
    %
    % vModule = griddedInterpolant(vSpline.profile(:,1),vSpline.profile(:,2),"linear","linear");
    % dvModule = griddedInterpolant(vSpline.dProfile(:,1),vSpline.dProfile(:,2),"linear","linear");
    % Position Vector p(s)
    %  pos = @(theta) [ r(theta) .* sin(theta) .* cos(phi(theta)); ... % X
    %                   r(theta) .* sin(theta) .* sin(phi(theta)); ... % Y
    %                   r(theta) .* cos(theta) ];                      % Z
    %
    %  % First Derivative Position
    %  dPosition = @(theta) [dr(theta) * sin(theta) * cos(phi(theta)) + r(theta) * cos(theta) * cos(phi(theta)) - r(theta) * sin(theta) * sin(phi(theta)) * dPhi(theta); ... % dX
    %                    dr(theta) * sin(theta) * sin(phi(theta)) + r(theta) * cos(theta) * sin(phi(theta)) + r(theta) * sin(theta) * cos(phi(theta)) * dPhi(theta); ... % dY
    %                    dr(theta) * cos(theta) - r(theta) * sin(theta) ];  % dZ
    %
    % % Second Derivative Position
    %  ddPosition = @(theta) [  (ddr(theta)-r(theta)-r(theta)*dPhi(theta)^2)*sin(theta)*cos(phi(theta)) + 2*dr(theta)*cos(theta)*cos(phi(theta)) - 2*dr(theta)*sin(theta)*sin(phi(theta))*dPhi(theta) - 2*r(theta)*cos(theta)*sin(phi(theta))*dPhi(theta) - r(theta)*sin(theta)*sin(phi(theta))*ddPhi(theta); ... % d2X
    %                      (ddr(theta)-r(theta)-r(theta)*dPhi(theta)^2)*sin(theta)*sin(phi(theta)) + 2*dr(theta)*cos(theta)*sin(phi(theta)) + 2*dr(theta)*sin(theta)*cos(phi(theta))*dPhi(theta) + 2*r(theta)*cos(theta)*cos(phi(theta))*dPhi(theta) + r(theta)*sin(theta)*cos(phi(theta))*ddPhi(theta); ... % d2Y
    %                      ddr(theta)*cos(theta) - 2*dr(theta)*sin(theta) - r(theta)*cos(theta) ]; %d2Z
    %
    %  % Package structure for ODE
    %  position.position = pos ;
    %  position.dPosition = dPosition ;
    %  position.ddPosition = ddPosition ;
    %


    
% discretizationLength = pi/200 ;
nIntervals = length(thetaCoord) - 1 ;
nPointsPerInterval = 200 ;

theta = zeros(1,nPointsPerInterval + (nPointsPerInterval-1) * (nIntervals - 1)) ;
theta(1:nPointsPerInterval) = linspace(thetaCoord(1), thetaCoord(2), nPointsPerInterval) ;

for ii = 2 : nIntervals
    thetaInterval = linspace(thetaCoord(ii), thetaCoord(ii+1), nPointsPerInterval) ; 
    theta(1+nPointsPerInterval+(nPointsPerInterval - 1) * (ii-2) : nPointsPerInterval + (nPointsPerInterval-1)*(ii-1)) = thetaInterval(2:end) ; 
end

[rSpline] = splineGenerationEfficient(rPoints, nPointsPerInterval) ;
r = rSpline.profile(:,2)' ; 
dr = rSpline.dProfile(:,2)' ; 
ddr = rSpline.ddProfile(:,2)' ; 


[phiSpline] = splineGenerationEfficient(phiPoints, nPointsPerInterval) ;
phi = phiSpline.profile(:,2)' ; 
dPhi = phiSpline.dProfile(:,2)' ; 
ddPhi = phiSpline.ddProfile(:,2)' ; 


[vSpline] = splineGenerationEfficient(vModulePoints, nPointsPerInterval) ;
vModule = vSpline.profile(:,2)' ; 
dvModule = vSpline.dProfile(:,2)' ; 


pos = [ r .* sin(phi) .* cos(theta) ;   %X
        r .* sin(phi) .* sin(theta) ;   %Y
        r .* cos(phi)];               %Z

dPosition = [dr .* sin(phi) .* cos(theta) + r .* cos(theta) .* cos(phi) .*dPhi - r .* sin(theta) .* sin(phi) ; 
             dr .* sin(theta) .* sin(phi) + r .* cos(theta) .* sin(phi) + r .* sin(theta) .* cos(phi) .* dPhi ;
             dr .* cos(phi) - r .* sin(phi).* dPhi];

ddPosition = [(ddr - r  - r .* (dPhi).^2) .* sin(phi) .* cos(theta) + 2 .* dr .* cos(theta) .* cos(phi) .* dPhi - 2 .* dr .* sin(theta) .* sin(phi) - 2 .* r .* cos(phi) .* sin(theta) .* dPhi + r .* cos(theta) .* cos(phi) .*ddPhi ;
              (ddr - r - r .* (dPhi).^2) .* sin(theta) .* sin(phi) + 2 .* dr .* cos(theta) .* sin(phi) + 2 .* dr .* sin(theta) .* cos(phi) .* dPhi + 2 .* r .* cos(theta) .* cos(phi) .* dPhi + r .* sin(theta) .* cos(phi) .* ddPhi ; 
              ddr .* cos(phi) - 2 .* dr .* sin(phi).*dPhi - r .* cos(phi) .*dPhi.^2 - r.*sin(phi).*ddPhi];

     % Package structure for ODE
     position.position = pos ;
     position.dPosition = dPosition ;
     position.ddPosition = ddPosition ;

    
    trajectoryLength = zeros(1, length(pos(:,1))) ; 

    for ii = 1:length(pos(1,:))-1
    
        trajectoryLength(ii) = sqrt((pos(1,ii) - pos(1,ii+1)).^2 + (pos(2,ii) - pos(2,ii+1)).^2 + (pos(3,ii) - pos(3,ii+1)).^2) ; 
        
    end

    totalSLength = sum(trajectoryLength) ;
    sVec = cumsum(trajectoryLength) ;

    %sSpan = [0, 100] ;
    %thetaSpan = [thetaCoord(1), thetaCoord(end)] ;
    
    opt = odeset('relTol', 1e-6, 'absTol', 1e-2, 'Events' , @(s,m)arrestingEvents(s,m,position,sVec,xf)) ;
    

    % Run Integration
    % [theta, m]=ode45(@(theta,m) massIntegration(theta, m, vSpline, position, mission), theta, m0, opt) ; 
    [sVec, m]=ode45(@(s,m) massIntegration(s, m, vSpline, position, sVec, mission), sVec, m0, opt) ; 
    
    % Pre Allocation of Variables
    vel = zeros(1,length(theta));
    velocityUnitVector = zeros(3,length(theta));
    accelerationVector = zeros(3,length(theta));
    accelerationModule = zeros(1, length(theta));
    thrust = zeros(3,length(theta));
    
    
    for ii = 1:length(theta)

        % Atmospheric Density (Exponential approximation for better physics)
        h = norm(pos(:,ii)) - 6371000; 
        rho = 1.225 * exp(-h/7200); 

        % Velocity Module
        vel(ii) = vModule(ii) ;

        % Velocity Direction
        velocityUnitVector(:,ii) = dPosition(:,ii) ./ norm(dPosition(:,ii)) ; 

        % Acceleration Direction
        accelerationVector(:,ii) = dvModule(ii) * vModule(ii) / norm(dPosition(:,ii)) * velocityUnitVector(:,ii) + vModule(ii) * vModule(ii) / norm(dPosition(:,ii)) * ( (ddPosition(:,ii) * norm(dPosition(:,ii)) - dPosition(:,ii) * dot(dPosition(:,ii), ddPosition(:,ii)) / norm(dPosition(:,ii))) /norm(dPosition(:,ii))^2);

        % Acceleration Module
        accelerationModule(:,ii) = norm( accelerationVector(:,ii) ) ;

        % Gravity Vector
        gVector = - GM * pos(:,ii) /norm(pos(:,ii))^3;

        % Thrust Evaluation
        thrust(:,ii) = m(ii) * accelerationVector(:,ii) - m(ii) * gVector + 0.5 * rho * vModule(ii).^2 * surface * CD .* velocityUnitVector(:,ii) ;
    end   

    prevSimulation = x ;
    prevVel = vel ; 
    prevRCoord = rCoord ; 
    prevPhiCoord = phiCoord ; 
    prevTheta = theta ;
    prevAccelerationModule = accelerationModule ; 
    prevPosition = position ;
    prevThrust = thrust ; 
    prevVelocityUnitVector = velocityUnitVector ;
    prevTrajectoryLength = trajectoryLength ; 


    
end


    if fitnessFlag 

        output.vel      = vel ;
        output.rCoord   = rCoord ;
        output.phiCoord = phiCoord ;
        output.theta    = theta ; 
        output.position           = position.position ; 
        output.trajectoryLength = trajectoryLength ;
    
    else

        output.acceleration       = accelerationModule; 
        output.position           = position.position ; 
        output.thrust             = thrust ; 
        output.theta              = theta ; 
        output.velocityUnitVector = velocityUnitVector ;
        output.vel                = vel ;
    
    end

end

