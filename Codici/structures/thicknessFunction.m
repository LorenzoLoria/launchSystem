function [t, tMax, stressMatrix, mStruct] = thicknessFunction(mission, nComponents)

% Function required to size the launcher thickess when tanks are
% pressurized during the flight. Evaluation must be done in the most
% critical conditions: q, q-alpha, MECO 

% --- INPUTS
% mission = struct containing LV data

% --- OUTPUT
% t = vector containing thickness of each component [m]
% tMax = maximum thickness of the launcher [m]
% stressMatrix = (6 x nComponents) matrix. Each column represents a
% component of the LV (e.g. Thrust structure, 1st stage etc.), each row
% represents a different kind of stress
% mStruct = vector containing mass of each structure [kg]

% ============================== DATA =====================================
% M              = mission.launcher.mass; % total mass of the element considered [kg]
r              = mission.launcher.diameter / 2; % radius of cylindrical section [m]
h              = mission.launcher.length; % length of the cylinder [m]
% hCM            = mission.launcher.hCM; % position of the center of mass [m]
nx             = launcher.nx; % longitudinal load factor
% nz             = launcher.nz; % lateral load factor
g0             = mission.environment.g0; 
SF             = mission.launcher.structures.SF; % safety factor
ultimateStress = mission.launcher.structures{ii}.ultimate; % ultimate stress for chosen material [Pa]
E              = mission.launcher.structures{ii}.E; % Young Modulus for chosen material [Pa]
shearAllowable = ultimateStress / 2; % ultimate shear stress for chosen material [Pa]
rhoOX          = mission.launcher.engines{ii}.oxDens; % density of the oxidizer [kg/m^3]
rhoFu          = mission.launcher.engines{ii}.fuDens; % density of the fuel [kg/m^3]
rhoMaterial    = mission.launcher.structures{ii}.rho; % density of the chosen material [kg/m^3] 
p              = mission.launcher.tankPressure; % pressure of tanks [kg/m^3]
longitudinalLoad = mission.launcher.structures.N; % Axial Load [N] (from loadFinder)
lateralLoad    = mission.launcher.structures.T; % Shear Load [N] (from loadFinder)
bendingMoment  = mission.launcher.structures.Mb; % Bending Moment [Nm] (from loadFinder)

% ============================ Solution ===================================

for i = 1 : nComponents
    if p(i) ~= 0
        % % Loads
        % longitudinalLoad = nx .* M * g0; % longitudinal load vector
        % lateralLoad = nz .* M .* g0; % lateral force vector
        % bendingMoment = lateralLoad .* hCM; % bending moment vector
        
        % Minimum Allowable Thicknesses
        tAxial = abs(- longitudinalLoad(i) / (2 * pi * r(i)) - bendingMoment(i) / (pi * r(i)^2) + p(i) * r(i) / 2 + rhoOX * nx * g0 * r(i) / 2) / ultimateStress(i) * SF;
        tShear = lateralLoad(i) / (2 * pi * r(i) * shearAllowable) * SF;
        
        
        if tAxial > tShear
            t(i) = tAxial;
        else
            t(i) = tShear;
        end

        % Area 
        A = pi * (r(i)^2 - (r(i)-t(i)).^2);
        
        % Inertia (valid for a cylinder)
        I = pi / 4 * (r(i)^4 - (r(i) - t(i)).^4);
        
        % Mass computation
        volume = pi .* h(i) .* (r(i)^2 - (r(i) - t(i)).^2);
        mStruct(i) = volume .* rhoMaterial(i);
        
        % Buckling for Pressurized Cylinders
        k0 = 9 * (t(i) / r(i)).^(0.6) + 0.16 .* (r(i) / h(i)).^(1.3) .* (t(i) / r(i)).^(0.3);
        kp = 0.191 * (p(i) / E(i)) .* (r(i) / t(i)).^2;
        
        
        if kp(ii) < 0.229
        else
            kp(ii) = 0.229;
        end
        
        
        sigmaBuckling = ((k0 + kp) * E .* t) / r;
        
        % Stresses (MAGNITUDE)
        longitudinalStress = longitudinalLoad / A * SF; 
        bendingStress = bendingMoment ./ I .* r * SF;
        shearStress = lateralLoad ./ A * SF;
        hoopStress = (p .* r + pHydro .* r) ./ t;
        axialStress = hoopStress / 2;
        
        stressMatrix(i, :) = [longitudinalStress, bendingStress, shearStress, hoopStress, axialStress, sigmaBuckling];

    else 

    % Loads
    % longitudinalLoad = nx .* M * g0; % longitudinal load vector
    % lateralLoad = nz .* M .* g0; % lateral force vector
    % bendingMoment = lateralLoad .* hCM; % bending moment vector
    
    % Area A and Inertia I are initially computed using the following
    % approximation:
    % A = 2 * pi * r * t;
    % I = pi * r^3 * t;
    
    % Minimum Allowable Thicknesses
    tAxial = abs(- longitudinalLoad(i) / (2 * pi * r(i)) - bendingMoment(i) / (pi * r(i)^2)) / ultimateStress(i) * SF;
    tShear = lateralLoad / (2 * pi * r * shearAllowable) * SF;
    
    
    if tAxial > tShear
        t(i) = tAxial;
    else
        t(i) = tShear;
    end
    
    % Area 
    A = pi * (r(i)^2 - (r(i)-t(i)).^2);
    
    % Inertia
    I = pi / 4 * (r(i)^4 - (r(i) - t(i)).^4);
    
    % Stresses (MAGNITUDE)
    longitudinalStress = longitudinalLoad / A * SF; 
    bendingStress = bendingMoment ./ I .* r * SF;
    shearStress = lateralLoad ./ A * SF;
    
    % Buckling for UnPressurized Cylinders
    
    sigmaBuckling = 9 * (t ./ r)^1.6 + 0.16 * (t ./ h)^1.3;
    
    stressMatrix(i, :) = [longitudinalStress, bendingStress, shearStress, 0, 0, sigmaBuckling];

    % Exctraction of the most critical result
    volume = pi .* h * (r(i)^2 - (r(i) - t(i)).^2);
    mStruct(i) = volume .* rhoMaterial;
    end
end
tMax = max(t);
end

% --- Data
% mass = 50000;
% height = 8;
% hCM = height / 2;
% radius = 2.5 / 2;
% thickness = 1e-3;
% loadFactorX = 2.5; % depends on the most critical load condition
% loadFactorZ = 0.2; % depends on the most critical load condition
% g0 = 9.81;
% safeyFactor = 1.25; % typically between 1.1 and 1.5
% sigmaAllowableSteel = 1034e6;
% sigmaAllowableAl = 448e6;
% shearAllowableSteel = sigmaAllowableSteel / 2;
% shearAllowableAl = sigmaAllowableAl / 2;
% stiffnessSteel = 207e9;
% stiffnessAl = 69e9;