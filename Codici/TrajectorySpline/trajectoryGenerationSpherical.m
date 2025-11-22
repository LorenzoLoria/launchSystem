function [output] = trajectoryGenerationSpherical(x, thetaCoord, x0, target, mission, fitnessFlag)
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


persistent prevSimulation  prevVel prevRCoord prevPhiCoord prevTheta prevAccelerationModule
persistent prevPosition prevThrust prevVelocityUnitVector prevTrajectoryLength prevSVec unpheasibleFlag prevTBallisticFall

if isempty(prevSimulation)
    prevSimulation = zeros(1,length(x)) ;
end

if x == prevSimulation



    if unpheasibleFlag

        output.value = 1e9 ;
        output.unpheasibleFlag = unpheasibleFlag ;
        return

    else

        vel = prevVel ;
        rCoord = prevRCoord ;
        phiCoord = prevPhiCoord ;
        theta = prevTheta ;
        accelerationModule = prevAccelerationModule ;
        position = prevPosition ;
        thrust = prevThrust ;
        velocityUnitVector = prevVelocityUnitVector ;
        trajectoryLength = prevTrajectoryLength ;
        sVec = prevSVec;
        tBallisticFall = prevTBallisticFall ;

        output.unpheasibleFlag = unpheasibleFlag ;

        if fitnessFlag

            output.vel      = vel ;
            output.rCoord   = rCoord ;
            output.phiCoord = phiCoord ;
            output.theta    = theta ;
            output.sVec             = sVec ;
            output.position         = position.position ;
            output.trajectoryLength = trajectoryLength ;
            output.tBallisticFall   = tBallisticFall ;


        else

            output.acceleration       = accelerationModule;
            output.position           = position.position ;
            output.thrust             = thrust ;
            output.sVec               = sVec ;
            output.velocityUnitVector = velocityUnitVector ;
            output.vel                = vel ;

        end
    end



else

    surface = mission.capsule.Area;
    CD      = mission.capsule.supersonicCD;
    GM      = mission.environment.GM;
    m0      = mission.launcher.engines{1}.m0 ;


    rCoord   = [x0(1) x(1:(length(x)-1)/3)];
    phiCoord = [x0(3) x((length(x)-1)/3+1 : (length(x)-1)/3*2)];
    vModule  = [10 x((length(x)-1)/3*2+1 : end-1)];

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
        r .* cos(phi)];       %Z

    position.position = pos ;

    trajectoryLength = zeros(1, length(pos(:,1))) ;

    for ii = 1:length(pos(1,:))-1

        trajectoryLength(ii) = sqrt((pos(1,ii) - pos(1,ii+1)).^2 + (pos(2,ii) - pos(2,ii+1)).^2 + (pos(3,ii) - pos(3,ii+1)).^2) ;

    end

    sVec = cumsum(trajectoryLength) ;


    dPosition = zeros(3, length(sVec));
    ddPosition = zeros(3, length(sVec));
    dvModule = zeros(1, length(sVec)) ;

    for ii = 2 : length(sVec) - 1

        deltaS = trajectoryLength(ii);

        dPosition(:, ii) = (pos(:, ii + 1) - pos(:, ii - 1)) / (2 * deltaS) ;
        ddPosition(:,ii) = (pos(:, ii + 1)  - 2 * pos(:, ii) + pos(:, ii - 1)) / deltaS^2 ;
        dvModule(ii)   = (vModule(ii + 1) - vModule(ii - 1)) / (2 * deltaS) ;

    end

    dPosition(:, 1) = (-1.5 * pos(:,1) + 2 * pos(:,2) - 0.5 * pos(:,3)) / trajectoryLength(1) ;
    dPosition(:, end) = (1.5 * pos(:,end) - 2 * pos(:,end-1) + 0.5 * pos(:,end-2)) / trajectoryLength(end) ;
    ddPosition(:, 1) = ( 2 * pos(:,1) - 5 * pos(:, 2) + 4 * pos(:, 3) - pos(:,4)) / trajectoryLength(1)^2;
    ddPosition(:, end) = (-2 * pos(:,end) + 5 * pos(:, end-1) - 4 * pos(:, end-2) + pos(:,end-3)) / trajectoryLength(end)^2;
    dvModule(1) = (-1.5 * vModule(1) + 2 * vModule(2) - 0.5 * vModule(3)) / trajectoryLength(1) ;
    dvModule(end) = (+1.5 * vModule(end) - 2 * vModule(end-1) + 0.5 * vModule(end-2)) / trajectoryLength(end) ;

    position.dPosition = dPosition;
    position.ddPosition = ddPosition;

    terminalVelocityDirection = (dPosition(:,end)) / norm(dPosition(:,end)) ;
    initialCapsuleCondition = [pos(1,end) pos(2,end) pos(3,end) terminalVelocityDirection' * vModule(end)] ;


    tRelease = 0 ;
    windDirection = [1 ; 0 ; 0] ;
    [tBallisticFall,ballisticFall] = ballisticTrajectory(initialCapsuleCondition,mission,windDirection,tRelease) ;


    if norm(ballisticFall(1:3,end) - target) > 20*1e3

        unpheasibleFlag = 1 ;

        output.value = 1e9 ;
        output.unpheasibleFlag = unpheasibleFlag ;

        return

    else

        unpheasibleFlag = 0 ;
        output.unpheasibleFlag = unpheasibleFlag ;




        opt = odeset('relTol', 1e-6, 'absTol', 1e-6, 'Events' , @(s,m)arrestingEvents(s,m)) ;


        % Run Integration
        % [theta, m]=ode45(@(theta,m) massIntegration(theta, m, vSpline, position, mission), theta, m0, opt) ;
        [sVec, m]=ode45(@(s,m) massIntegration(s, m, vSpline, dvModule, position, sVec, mission), sVec, m0, opt) ;

        % Pre Allocation of Variables
        vel = zeros(1,length(sVec));
        velocityUnitVector = zeros(3,length(sVec));
        accelerationVector = zeros(3,length(sVec));
        accelerationModule = zeros(1, length(sVec));
        thrust = zeros(3,length(sVec));


        for ii = 1:length(sVec)

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
        prevSVec = sVec ;
        prevTBallisticFall = tBallisticFall ;





        if fitnessFlag

            output.vel      = vel ;
            output.rCoord   = rCoord ;
            output.phiCoord = phiCoord ;
            output.theta    = theta ;
            output.sVec             = sVec ;
            output.position         = position.position ;
            output.trajectoryLength = trajectoryLength ;
            output.tBallisticFall   = tBallisticFall ;

        else

            output.acceleration       = accelerationModule;
            output.position           = position.position ;
            output.thrust             = thrust ;
            output.sVec               = sVec ;
            output.velocityUnitVector = velocityUnitVector ;
            output.vel                = vel ;

        end

    end


end

