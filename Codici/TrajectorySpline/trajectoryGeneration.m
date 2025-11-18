function [output] = trajectoryGeneration(x,xCoord, x0,x_target,mission,fitnessFlag)
    
    yCoord = [x0(2) x(1:length(x)/3) x_target(2)];
    zCoord = [x0(3) x(length(x)/3+1 : length(x)/3*2) x_target(3) ];
    vModule = [2500 x(length(x)/3*2+1 : end) x(end)];
    
    % Create matrices for spline generation
    yPoints = [xCoord ; yCoord]' ;
    zPoints = [xCoord ; zCoord]' ;
    vModulePoints = [xCoord ; vModule]' ;
    
    % Spline Generation for y z and vModule
    yQueryPoints = linspace(yCoord(1), yCoord(end), 100);
    [ySpline] = splineGeneration(yPoints,yQueryPoints);

    y = griddedInterpolant(ySpline.profile(:,1),ySpline.profile(:,2),"linear","linear");
    dy = griddedInterpolant(ySpline.dProfile(:,1),ySpline.dProfile(:,2),"linear","linear");
    ddy = griddedInterpolant(ySpline.ddProfile(:,1),ySpline.ddProfile(:,2),"linear","linear");


    zQueryPoints = linspace(zCoord(1), zCoord(end), 100);
    [zSpline] = splineGeneration(zPoints, zQueryPoints) ;
    
    z = griddedInterpolant(zSpline.profile(:,1),zSpline.profile(:,2),"linear","linear");
    dz = griddedInterpolant(zSpline.dProfile(:,1),zSpline.dProfile(:,2),"linear","linear");
    ddz = griddedInterpolant(zSpline.ddProfile(:,1),zSpline.ddProfile(:,2),"linear","linear");
    

    vQueryPoints = linspace(vModule(1), vModule(end), 100);
    [vSpline] = splineGeneration(vModulePoints, vQueryPoints) ;
    
    vModule = griddedInterpolant(vSpline.profile(:,1),vSpline.profile(:,2),"linear","linear");
    dvModule = griddedInterpolant(vSpline.dProfile(:,1),vSpline.dProfile(:,2),"linear","linear");
    
    % Estrapolation of position, velocity and acceleration
    pos = @(s) [ s ; y(s) ; z(s) ];
    dPosition = @(s) [1 ; dy(s) ; dz(s) ];
    ddPosition = @(s) [0 ; ddy(s) ; ddz(s)];
    
    position.position = pos ;
    position.dPosition = dPosition ; 
    position.ddPosition = ddPosition ;
    
    % Unpack Mission Data from Struct
    rEarth  = mission.environment.rEarth;
    surface = mission.capsule.Area;
    CD  = mission.capsule.supersonicCD;
    GM  = mission.environment.GM;
    g0  = mission.environment.g0;

    %sSpan = [0, 100] ;
    sSpan = [xCoord(1), xCoord(end)] ;
    
    % Initial Value for integration
    m0 = mission.launcher.engines{1}.m0;
    
    optODE = odeset('relTol', 1e-6, 'absTol', 1e-4, 'Events' , @(s,m)targetEvent(s,m,position,x_target), 'MaxStep', 0.5) ;
    
    [s, m]=ode45(@(s,m) massIntegration(s, m, vSpline, position), sSpan, m0, optODE) ; 
    
    vel = zeros(1,length(s));
    velocityUnitVector = zeros(3,length(s));
    accelerationVector = zeros(3,length(s));
    accelerationModule = zeros(1, length(s));
    thrust = zeros(3,length(s));
    
    % Loop through the simulation steps
    for ii = 1:length(s)

       % Interpolate air density based on current altitude
       h   = norm(pos(s(ii)))-mission.environment.rEarth;  
       rho = mission.environment.gridInterp(h);
    
       vel(ii) = vModule(s(ii)) ;
       velocityUnitVector(:,ii) = dPosition(s(ii)) ./ norm(dPosition(s(ii))) ; 
       accelerationVector(:,ii) = (dvModule(s(ii)) * vModule(s(ii)) / norm(dPosition(s(ii)))) * velocityUnitVector(:,ii) + (vModule(s(ii)) * vModule(s(ii)) / norm(dPosition(s(ii)))) * ( (ddPosition(s(ii)) * norm(dPosition(s(ii))) - dPosition(s(ii)) * (dot(dPosition(s(ii)), ddPosition(s(ii))) / norm(dPosition(s(ii))))) / (norm(dPosition(s(ii)))^2) );
       accelerationModule(:,ii) = norm( accelerationVector(:,ii) ) ;

       % Gravity contribution
       gVector = - GM * pos(s(ii)) /norm(pos(s(ii)))^3;
       gVector = g0 / ( 1 + pos(s(ii))/rEarth).^2;

       % Calculate required thrust
       thrust(:,ii) = m(ii) * accelerationVector(:,ii) - m(ii) * gVector + 0.5 * rho * (vModule(s(ii))).^2 * surface * CD .* velocityUnitVector(:,ii) ;
     
    end
    
    % Output different data depending on whether the function is called
    if fitnessFlag 
    
        output.vel    = vel ;
        output.yCoord = yCoord ;
        output.zCoord = zCoord ;
        output.s      = s ; 
    else
        output.velocityUnitVector = velocityUnitVector ;
        output.acceleration = accelerationModule; 
        output.position = position.position ; 
        output.thrust = thrust ; 
        output.s = s ; 
    end

end
