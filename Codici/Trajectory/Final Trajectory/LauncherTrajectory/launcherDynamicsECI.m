
function dsdt = launcherDynamicsECI(t, x,thrustData, mission,stageNumber,opt,option2D,dimensions,engineVec)

% LAUNCHERDYNAMICS  3D launcher equations of motion.
%   This function computes the time derivative of the state vector for a 
%   multistage rocket in a 3D Cartesian coordinate system.
%   Assumption: Flat Earth approximation.
%
%   State Vector x (7x1):
%     x(1) = pos_x   [m]   Position North/East relative
%     x(2) = pos_y   [m]   Position Cross-range
%     x(3) = pos_z   [m]   Vertical Position 
%     x(4) = vx      [m/s] Velocity x-component
%     x(5) = vy      [m/s] Velocity y-component
%     x(6) = vz      [m/s] Velocity z-component
%     x(7) = m       [kg]  Instantaneous total mass
%
%   thrustData:
%
%     ThrustData(:,1) = Tx
%     ThrustData(:,2) = Ty
%     ThrustData(:,3) = Tz
%
%   mission.capsule:
%     .Area   = Reference cross-sectional area [m^2]
%     .Cd     = Drag coefficient [-] 
%
%   mission.environment:
%     .altRange = Vector of altitudes [m]
%     .rhoVal   = Vector of air densities [kg/m^3]
%     .g0       = Standard gravity [m/s^2]
%
%   mission.prop:
%     .Isp      = Specific Impulse [s]

    % Unpack State Vector and Mission Data
    r     = x(1:3);
    rMag   = sqrt(r'*r);
    % Add mean Wind profile
    h      = rMag-mission.environment.rEarth;
    if h<=100e3
        vxWind = mission.environment.windXFun(h/1000);
        vyWind = mission.environment.windYFun(h/1000);
    else
        vxWind = 0;
        vyWind = 0;
    end

    wind = mission.target.Rfinal'*[vxWind;vyWind;0];

    v     = x(4:6);
    m     = x(7);
    
    vRel = v + wind;
    vMag = sqrt(vRel'*vRel); 
    
    rMag = sqrt(r'*r);
    
    A   = mission.capsule.Area;
    g0  = mission.environment.g0; 
    GM  = mission.environment.GM;
 
    %rho = mission.environment.rhoFun(h);
    rho = 1.29*exp(-h/8433);

    dynamicPressure = 0.5 * rho * vMag^2;
    soundspeed = mission.aerodynamics.soundspeedFun(h);
    Mach = vMag/soundspeed;
    
    if Mach == 0
        Cd = 0.01;
    else
    [~,Cd,~,~] = CLCDcomputation(Mach,0,dynamicPressure,1,mission,stageNumber,dimensions,engineVec);
    end
    optVar = thrustData(t); 
    

if option2D == 1
    gammaGimball =optVar(2)*pi/180 ;
    thetaGimball =0;
else
    gammaGimball =optVar(2)*pi/180;
    thetaGimball =optVar(3)*pi/180;
end

percVec = optVar(1);

if stageNumber == 1
    staticContribution = (101325-mission.environment.pressure(h))*opt.stage{stageNumber}.engine.effAreaZero;
    isp = opt.stage{stageNumber}.engine.ispZero;
    thrustValue = opt.stage{stageNumber}.engine.thrustZero;
else
    staticContribution = 0;
    isp = opt.stage{stageNumber}.engine.ispVac;
    thrustValue = opt.stage{stageNumber}.engine.thrustVacum;
end

ThrustBRF = percVec * engineVec(1) *(thrustValue+staticContribution)* [cos(thetaGimball)*cos(gammaGimball); cos(thetaGimball)*sin(gammaGimball); sin(thetaGimball)];
ThrustIRF = mission.target.Rfinal'*ThrustBRF;


    % Drag contribution
    D = - 0.5 .* rho .* vMag .* A .* Cd .* vRel; 

    % Gravity contribution
    G = - GM * r /rMag^3;

    % mass flow rate
    mDot = -( percVec *thrustValue* engineVec(1)) / (g0 * isp); 

    % Equation of motion
    dsdt = zeros(7,1);
    
    % Velocity derivatives (Position rates)

    dsdt(1:3) = v;

    % Acceleration derivatives (Velocity rates)
    dsdt(4:6) = (ThrustIRF + D ) / m + G;  

    
    % Mass derivative
    dsdt(7) = mDot; 


end