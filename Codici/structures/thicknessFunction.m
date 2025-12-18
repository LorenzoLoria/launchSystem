function [updatedStructuralMass, mStruct, tVec, idx] = thicknessFunction(mission, launcher, configuration, maxQData, internalActions)

% Function required to size the launcher thickess of all the different 
% components of the LV. Evaluation must be done in the most
% critical conditions: q, q-alpha, MECO

% ============================== DATA =====================================

r = [mission.capsule.radius];                  % radius of cylindrical section of each body [m] (VECTOR)
h = [mission.capsule.height];                         % length of each body [m] (VECTOR)
rhoProps = 0;                                      % density [kg/m^3]  

nStages = launcher(1) ;

for ii = 1:nStages
    enginesUsed{ii} =  mission.engines{launcher(1+ii)};
end

for ii = launcher(1):-1:1
    r = [r , configuration.geometry.stage{ii}.interstage.meanRadius, configuration.stage{ii}.fuelTankR,  ...
        configuration.stage{ii}.oxTankR, configuration.geometry.stage{ii}.radius] ;
    h = [ h, configuration.geometry.stage{ii}.interstage.length, configuration.stage{ii}.fuelTankH, configuration.stage{ii}.oxTankH,...
        configuration.geometry.stage{ii}.thrustFrame] ;
    
    rhoProps = [rhoProps, 0, enginesUsed{ii}.fuelDens, enginesUsed{ii}.oxDens, 0] ;

end


a              = maxQData.aMaxQ;
g0             = mission.environment.g0; 
nx             = a(1) / g0 + 1; 
SF             = mission.structure.safetyFactor; % safety factor

ultimateStress = mission.structure.ultimate; % ultimate stress for chosen material [Pa]
E              = mission.structure.E; % Young Modulus for chosen material [Pa]
shearAllowable = ultimateStress / 2; % ultimate shear stress for chosen material [Pa]
rhoMaterial    = mission.structure.rho; % density of the chosen material [kg/m^3] 
materialNumber = length(E);



N              = internalActions.N; % Axial Load [N] (from loadFinder) (VECTOR)
T              = internalActions.T; % Shear Load [N] (from loadFinder) (VECTOR)
M              = internalActions.M; % Bending Moment [Nm] (from loadFinder) (VECTOR)

% ============================ SOLUTION ===================================

% --- Associate the load to the correct component: the LV is divided into
% payload, interstage between payload and II stage, II stage (pressurized),
% interstage (II and I stage) and I stage (pressurized).

% --- Re-Definition of the loads

Nnew = N([4:2:10,13:2:19,22]);
Tnew = T([4:2:10,13:2:19,22]);
Mnew = M([4:2:10,13:2:19]);

Mtail = M([20:1:22]);
[~, idx] = max(abs(Mtail));
firstStageM = Mtail(idx);
Mnew(end+1) = firstStageM;

N = Nnew;
T = Tnew;
M = Mnew;

partsNumber = length(N);


pressurization = zeros(partsNumber, 1);

for ii = [3, 4, 7, 8]
    pressurization(ii) = mission.structure.tankPressure;
end


t = zeros(partsNumber, materialNumber);
tVec = zeros(partsNumber, 1);
mStruct = zeros(partsNumber, 1);
stressMatrix = zeros(partsNumber, 6);
idx = zeros(partsNumber, 1);

for i = 1 : partsNumber
    
    if pressurization(i) ~= 0
        for j = 1 : materialNumber
        % Minimum Allowable Thicknesses
        tAxial = abs((- N(i) / (2 * pi * r(i)) - M(i) / (pi * r(i)^2)) * SF + pressurization(i) * r(i) / 2 + rhoProps(i) * nx * g0 * r(i) * h(i) / 2) / ultimateStress(j);
        tShear = abs(T(i)) / (2 * pi * r(i) * shearAllowable(j)) * SF;

        t(i, j) = max(tAxial, tShear);

        % Consider buckling
        bucklingEq = @(t_var) ( ...
                    + abs(N(i) / (2 * pi * r(i) * t_var)) * SF ...
                    + abs(M(i) / (pi * r(i)^2 * t_var) ) * SF ...
                     )...
                    - ( ((9 * (t_var/r(i))^0.6 + 0.16 * (r(i)/h(i))^1.3 * (t_var/r(i))^0.3) ...
                    + min(0.191 * (pressurization(i)/E(j)) * (r(i)/t_var)^2, 0.229) ) * E(j) * t_var / r(i));
                    % - abs((rhoProps(i)) * nx * g0 * h(i) * r(i)) / (2 * t_var)) ...
                    % - abs((pressurization(i) * r(i)) / (2 * t_var))
        options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-39);
        tBuckling = fsolve(bucklingEq, 0.001,options);
        t(i, j) = max(t(i, j), tBuckling);
        end
        
        [tVec(i), idx(i)] = min( t(i, :) .* (t(i,:)>=mission.structure.minThickness) + mission.structure.minThickness .* (t(i,:)<mission.structure.minThickness) );
        % [tVec(i), idx(i)] = min( t(i, :) + inf * (t(i, :) < mission.structure.minThickness) );
        % if isinf(tVec(i))
        %     tVec(i) = mission.structure.minThickness;
        %     idx(i)= 1;
        % end
        % Hydrostatic pressure
        pHydro = rhoProps(i) .* nx .* g0 .* h(i);

        % Area of the cylinder
        A = pi * (r(i)^2 - (r(i)-tVec(i)).^2);
        
        % Inertia (valid for a cylinder)
        I = pi / 4 * (r(i)^4 - (r(i) - tVec(i)).^4);
        
        % Mass computation
        volume = pi .* h(i) .* (r(i)^2 - (r(i) - tVec(i)).^2) + 2*pi*(r(i)^3 - (r(i)-tVec(i))^3);
        mStruct(i) = volume .* rhoMaterial(idx(i));
        
        % Buckling for Pressurized Cylinders
        k0 = 9 * (tVec(i) / r(i)).^(0.6) + 0.16 .* (r(i) / h(i)).^(1.3) .* (tVec(i) / r(i)).^(0.3);
        kp = 0.191 * (pressurization(i) / E(idx(i))) .* (r(i) / tVec(i)).^2;
        
        sigmaBuckling = ((k0 + min(kp, 0.229)) * E(idx(i)) .* tVec(i)) / r(i);
        
        % if sigmaBuckling < ultimateStress
        %     error('Spessore troppo sottile')
        % end

        % Stresses (MAGNITUDE)
        longitudinalStress = N(i) / A; 
        bendingStress = M(i) ./ I .* r(i);
        shearStress = T(i) ./ A;
        hoopStress = (pressurization(i) .* r(i) + pHydro .* r(i)) ./ tVec(i);
        axialStress = hoopStress / 2;
        
        stressMatrix(i, :) = [longitudinalStress, bendingStress, shearStress, hoopStress, axialStress, sigmaBuckling];
        
    else 
        for j = 1 : materialNumber
        % Minimum Allowable Thicknesses
        tAxial = ( abs(N(i) / (2 * pi * r(i))) + abs(M(i) / (pi * r(i)^2)) ) / ultimateStress(j) * SF;
        tShear = T(i) / (2 * pi * r(i) * shearAllowable(j)) * SF;
    
        t(i, j) = max(tAxial, tShear);
    
        bucklingEq = @(t_var) ((abs(N(i)) / (pi*(r(i)^2-(r(i)-t_var)^2)) + abs(M(i)) / (pi/4*(r(i)^4-(r(i)-t_var)^4)) * r(i))) - E(j) * ( 9 * (t_var / r(i))^1.6 + 0.16 * (t_var / h(i))^1.3 ) ;
        %- 0.6 * (1 - 0.901 * (1 - exp(-1 / 16 * sqrt(r(i)/t_var)))) * E(j) * t_var / r(i);
            options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-6);
            tBuckling = fsolve(bucklingEq, t(i, j), options);
            
            t(i, j) = max(t(i, j), tBuckling);
        end
        

        [tVec(i), idx(i)] = min( t(i, :) .* (t(i,:)>=mission.structure.minThickness) + mission.structure.minThickness .* (t(i,:)<mission.structure.minThickness) );
        % [tVec(i), idx(i)] = min( t(i, :) + inf * (t(i, :) < mission.structure.minThickness) );
        % if isinf(tVec(i))
        %     tVec(i) = mission.structure.minThickness;
        %     idx(i)= 1;
        % end
        % Area 
        A = pi * (r(i)^2 - (r(i)-tVec(i)).^2);
        
        % Inertia
        I = pi / 4 * (r(i)^4 - (r(i) - tVec(i)).^4) ;
        
        % Stresses (MAGNITUDE)
        longitudinalStress = N(i) / A; 
        bendingStress = M(i) ./ I .* r(i);
        shearStress = T(i) ./ A;
        
        % Buckling for UnPressurized Cylinders
        
        sigmaBuckling = 0.6 * (1 - 0.901 * (1 - exp(-1 / 16 * sqrt(r(i)/tVec(i))))) * E(idx(i)) * tVec(i) / r(i); % NASA SP-8007
        
        stressMatrix(i, :) = [longitudinalStress, bendingStress, shearStress, 0, 0, sigmaBuckling];
    
        % Exctraction of the most critical result
        volume = pi .* h(i) * (r(i)^2 - (r(i) - tVec(i)).^2);
        mStruct(i) = volume .* rhoMaterial(idx(i));
    end
end

% =========================== ESTRAZIONE ==================================
updatedStructuralMass = sum(mStruct) ;

% Extract the results associated to interstage2, fuel2, ox2, interstage1,
% fuel1, ox1
mStruct = mStruct([2, 3, 4, 6, 7, 8]);
tVec = tVec([2, 3, 4, 6, 7, 8]);
idx = idx([2, 3, 4, 6, 7, 8]);
end