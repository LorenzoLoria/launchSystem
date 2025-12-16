function dsdt = launcherDynamicsECIControlled(t, x,refStateFun, mission,stageNumber,opt,option2D,dimensions,engineVec,windVelXFun,windVelYFun,maxGimball,thrustData,gainGA,finsVec)
    

    % stateCollocation è quella del singolo stadio

    r      = x(1:3);

    rMag   = sqrt(r'*r);
    h      = rMag-mission.environment.rEarth; 
    if h<=100e3
        vxWind = windVelXFun(h/1000);
        vyWind = windVelYFun(h/1000);
    else
        vxWind = 0;
        vyWind = 0;
    end
    wind = mission.target.Rfinal'*[vxWind;vyWind;0];
    v      = x(4:6) ;

    m      = x(7);
    
    vRel = v + wind;
    vMag = sqrt(vRel'*vRel); 

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
    [~,Cd,~,~] = CLCDcomputation(Mach,0,dynamicPressure,1,mission,stageNumber,dimensions,engineVec,finsVec);
    end 

    if stageNumber == 1
        staticContribution = (101325-mission.environment.pressure(h))*opt.stage{stageNumber}.engine.effAreaZero;
        isp = opt.stage{stageNumber}.engine.ispZero;
        thrustValue = opt.stage{stageNumber}.engine.thrustZero;
        thrustMax = (thrustValue + staticContribution) * opt.stage{stageNumber}.nEngines;
    else
        staticContribution = 0;
        isp = opt.stage{stageNumber}.engine.ispVac;
        thrustValue = opt.stage{stageNumber}.engine.thrustVacum;
        thrustMax = (thrustValue + staticContribution) * opt.stage{stageNumber}.nEngines;
    end
    
    gain = gainGA([1:6]+6*(stageNumber-1))';

    if t < 15
        optVar = thrustData(t); 
        
    
        if option2D == 1
            gammaGimball = optVar(2)*pi/180 ;
            thetaGimball = 0;
        else
            gammaGimball =optVar(2)*pi/180;
            thetaGimball =optVar(3)*pi/180;
        end
        
        percVec = optVar(1);   
        ThrustBRF = percVec * engineVec(1) *(thrustValue+staticContribution) * [cos(thetaGimball)*cos(gammaGimball); cos(thetaGimball)*sin(gammaGimball); sin(thetaGimball)];
        ThrustIRF = mission.target.Rfinal'*ThrustBRF;
    else

        err = refStateFun(t+5)' - x(1:6);
    
        
        ThrustIRF = m * (gain(1:3) .* err(1:3) + gain(4:6) .* err(4:6)); %+ gain(7) .* err(7);
        dirThrust = ThrustIRF/sqrt(ThrustIRF' * ThrustIRF);

        if sqrt(ThrustIRF' * ThrustIRF) > thrustMax
            ThrustIRF = thrustMax .* dirThrust;
        end

        
        % compute angle wrt velocity
        %angleWrtVel = acosd(ThrustIRF'*v/(norm(ThrustIRF)*norm(v)));

        % if angleWrtVel > maxGimball
        %     angleWrtVel = maxGimball; %%%%% completare %%%%%
        % end

    end

    % Drag contribution
    D = - 0.5 .* rho .* vMag .* A .* Cd .* vRel; 

    % Gravity contribution
    G = - GM * r /rMag^3;

    % mass flow rate
    mDot = -( sqrt(ThrustIRF' * ThrustIRF)) / (g0 * isp); 

    % Equation of motion
    dsdt = zeros(7,1);
    
    % Velocity derivatives (Position rates)

    dsdt(1:3) = v;

    % Acceleration derivatives (Velocity rates)
    dsdt(4:6) = (ThrustIRF + D ) / m + G;  

    
    % Mass derivative
    dsdt(7) = mDot;

end