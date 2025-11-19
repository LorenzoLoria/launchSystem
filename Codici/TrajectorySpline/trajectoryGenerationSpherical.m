function [output] = trajectoryGenerationSpherical(x, sCoord, target, x0, xf, mission, fitnessFlag)
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
    
    rCoord   = [x0(1) x(1:length(x)/3) xf(1)];
    phiCoord = [x0(3) x(length(x)/3+1 : length(x)/3*2) xf(3) ];
    vModule  = [0.1 x(length(x)/3*2+1 : end) x(end)];
    
    % Prepare data for spline generation
    rPoints       = [sCoord ; rCoord]' ;
    phiPoints     = [sCoord ; phiCoord]' ;
    vModulePoints = [sCoord ; vModule]' ;

    rQueryPoints = linspace(rCoord(1), rCoord(end), 100);
    [rSpline] = splineGeneration(rPoints, rQueryPoints) ;
    
    r = griddedInterpolant(rSpline.profile(:,1),rSpline.profile(:,2),"linear","linear");
    dr = griddedInterpolant(rSpline.dProfile(:,1),rSpline.dProfile(:,2),"linear","linear");
    ddr = griddedInterpolant(rSpline.ddProfile(:,1),rSpline.ddProfile(:,2),"linear","linear");
    
    phiQueryPoints = linspace(phiCoord(1), phiCoord(end), 100);
    [phiSpline] = splineGeneration(phiPoints, phiQueryPoints) ;
    
    phi = griddedInterpolant(phiSpline.profile(:,1),phiSpline.profile(:,2),"linear","linear");
    dphi = griddedInterpolant(phiSpline.dProfile(:,1),phiSpline.dProfile(:,2),"linear","linear");
    ddphi = griddedInterpolant(phiSpline.ddProfile(:,1),phiSpline.ddProfile(:,2),"linear","linear");
    
    vQueryPoints = linspace(vModule(1), vModule(end), 100);
    [vSpline] = splineGeneration(vModulePoints, vQueryPoints) ;
    
    vModule = griddedInterpolant(vSpline.profile(:,1),vSpline.profile(:,2),"linear","linear");
    dvModule = griddedInterpolant(vSpline.dProfile(:,1),vSpline.dProfile(:,2),"linear","linear");
    
    % Position Vector p(s)
    pos = @(s) [ r(s) .* sin(s) .* cos(phi(s)); ... % X
                 r(s) .* sin(s) .* sin(phi(s)); ... % Y
                 r(s) .* cos(s) ];                  % Z

    % First Derivative Position
    dPosition = @(s) [dr(s) * sin(s) * cos(phi(s)) + r(s) * cos(s) * cos(phi(s)) - r(s) * sin(s) * sin(phi(s)) * dphi(s); ... % dX
                      dr(s) * sin(s) * sin(phi(s)) + r(s) * cos(s) * sin(phi(s)) + r(s) * sin(s) * cos(phi(s)) * dphi(s); ... % dY
                      dr(s) * cos(s) - r(s) * sin(s) ];  % dZ 
   
   % Second Derivative Position
    ddPosition = @(s) [  (ddr(s)-r(s)-r(s)*dphi(s)^2)*sin(s)*cos(phi(s)) + 2*dr(s)*cos(s)*cos(phi(s)) - 2*dr(s)*sin(s)*sin(phi(s))*dphi(s) - 2*r(s)*cos(s)*sin(phi(s))*dphi(s) - r(s)*sin(s)*sin(phi(s))*ddphi(s); ... % d2X
                        (ddr(s)-r(s)-r(s)*dphi(s)^2)*sin(s)*sin(phi(s)) + 2*dr(s)*cos(s)*sin(phi(s)) + 2*dr(s)*sin(s)*cos(phi(s))*dphi(s) + 2*r(s)*cos(s)*cos(phi(s))*dphi(s) + r(s)*sin(s)*cos(phi(s))*ddphi(s); ... % d2Y
                        ddr(s)*cos(s) - 2*dr(s)*sin(s) - r(s)*cos(s) ]; %d2Z
    
    % Package structure for ODE
    position.position = pos ;
    position.dPosition = dPosition ; 
    position.ddPosition = ddPosition ; 
    
    iSp = 300 ;
    surface = 1 ;
    CD = 0.2 ; 
    
    %sSpan = [0, 100] ;
    sSpan = [sCoord(1), sCoord(end)] ;
    
    m0 = 90 ;
    
    opt = odeset('relTol', 1e-6, 'absTol', 1e-2, 'Events' , @(s,m)targetEvent(s,m,position,target), 'MaxStep', 0.5) ;
    
    % Run Integration
    [s, m]=ode45(@(s,m) massIntegration(s, m, vSpline, position, mission), sSpan, m0, opt) ; 
    
    % Pre Allocation of Variables
    vel = zeros(1,length(s));
    velocityUnitVector = zeros(3,length(s));
    accelerationVector = zeros(3,length(s));
    accelerationModule = zeros(1, length(s));
    thrust = zeros(3,length(s));
    
    
    for ii = 1:length(s)

        % Atmospheric Density (Exponential approximation for better physics)
        h = norm(pos(s(ii))) - 6371000; 
        rho = 1.225 * exp(-h/7200); 

        % Velocity Module
        vel(ii) = vModule(s(ii)) ;

        % Velocity Direction
        velocityUnitVector(:,ii) = dPosition(s(ii)) ./ norm(dPosition(s(ii))) ; 

        % Acceleration Direction
        accelerationVector(:,ii) = (dvModule(s(ii)) * vModule(s(ii)) / norm(dPosition(s(ii)))) * velocityUnitVector(:,ii) + (vModule(s(ii)) * vModule(s(ii)) / norm(dPosition(s(ii)))) * ( (ddPosition(s(ii)) * norm(dPosition(s(ii))) - dPosition(s(ii)) * (dot(dPosition(s(ii)), ddPosition(s(ii))) / norm(dPosition(s(ii))))) / (norm(dPosition(s(ii)))^2));

        % Acceleration Module
        accelerationModule(:,ii) = norm( accelerationVector(:,ii) ) ;

        % Gravity Vector
        gVector = - GM * pos(s(ii)) /norm(pos(s(ii)))^3;

        % Thrust Evaluation
        thrust(:,ii) = m(ii) * accelerationVector(:,ii) - m(ii) * gVector + 0.5 * rho * (vModule(s(ii))).^2 * surface * CD .* velocityUnitVector(:,ii) ;
    end   
    
    if fitnessFlag 
        output.vel      = vel ;
        output.rCoord   = rCoord ;
        output.phiCoord = phiCoord ;
        output.s        = s ; 
    
    else
        output.acceleration       = accelerationModule; 
        output.position           = position.position ; 
        output.thrust             = thrust ; 
        output.s                  = s ; 
        output.velocityUnitVector = velocityUnitVector ;
    
    end

end

