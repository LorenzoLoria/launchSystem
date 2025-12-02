function [mission] = thicknessFunction(mission, ii)

% Function required to size the launcher thickess of all the different 
% components of the LV. Evaluation must be done in the most
% critical conditions: q, q-alpha, MECO 

% --- INPUTS
% mission       = struct containing LV data
% nComponents   = number of components in which the LV is divided

% --- OUTPUTS
% mStruct       = vector containing mass of each structure [kg]
% t             = vector containing thickness of each component [m]
% tMax          = maximum thickness of the launcher [m]
% stressMatrix  = (nComponents x 6) matrix. Each column represents a
%                component of the LV (e.g. Thrust structure, 1st stage etc.), 
%                each row represents a different kind of stress [MPa]


% ============================== DATA =====================================
nComponents    = mission.structure.nComponents;

r              = mission.structure.diameter / 2; % radius of cylindrical section of each body [m] (VECTOR)
h              = mission.structure.elementLength; % length of each body [m] (VECTOR)
a              = mission.structure.aMaxQ;
g0             = mission.environment.g0; 
nx             = a(1) / g0; 
SF             = mission.structure.SF; % safety factor

ultimateStress = mission.structure.ultimate; % ultimate stress for chosen material [Pa]
E              = mission.structure.E; % Young Modulus for chosen material [Pa]
shearAllowable = ultimateStress / 2; % ultimate shear stress for chosen material [Pa]
rhoMaterial    = mission.structure.rho; % density of the chosen material [kg/m^3] 

rhoOX          = mission.launcher.engines{ii}.oxDens; % density of the oxidizer [kg/m^3]
rhoFu          = mission.launcher.engines{ii}.fuelDens; % density of the fuel [kg/m^3]

p              = mission.structure.tankPressure; % pressure of component [Pa] [nComponents x 1]
N              = mission.structure.N; % Axial Load [N] (from loadFinder) (VECTOR)
T              = mission.structure.T; % Shear Load [N] (from loadFinder) (VECTOR)
M              = mission.structure.M; % Bending Moment [Nm] (from loadFinder) (VECTOR)

% ============================ SOLUTION ===================================

% --- Associate the load to the correct component: the LV is divided into
% payload, interstage between payload and II stage, II stage (pressurized),
% interstage (II and I stage) and I stage (pressurized).

% --- Re-Definition of the loads
payloadN = N(4);
payloadT = T(4);
payloadM = M(4);

interN = [];
interT = [];
interM = [];

for ii=6:2:size(N)-3
    interN = [interN, N(ii)];
    interT = [interT, T(ii)];
    interM = [interM, M(ii)];
end

firstStageN = N(end);
firstStageT = T(end);
firstStageM = M(end);

newN = [payloadN, interN, firstStageN];
newT = [payloadT, interT, firstStageT];
newM = [payloadM, interM, firstStageM];


t = zeros(nComponents, 1);
mStruct = zeros(nComponents, 1);
stressMatrix = zeros(nComponents, 6);

for i = 1 : nComponents
    if p(i) ~= 0
        
        % Minimum Allowable Thicknesses
        tAxial = abs(- N(i) / (2 * pi * r(i)) - M(i) / (pi * r(i)^2) + p(i) * r(i) / 2 + rhoOX * nx * g0 * r(i) * h(i) / 2) / ultimateStress * SF;
        tShear = T(i) / (2 * pi * r(i) * shearAllowable) * SF;

        t(i) = max(tAxial, tShear);

        % Consider buckling
        bucklingEq = @(t_var) ((N(i) / (pi*(r(i)^2-(r(i)-t_var)^2)) + M(i) / (pi/4*(r(i)^4-(r(i)-t_var)^4)) * r(i)) - ((p(i) * r(i) + rhoOX * nx * g0 * h(i) * r(i)) / (2 * t_var))) - (((9 * (t_var/r(i))^0.6 + 0.16 * (r(i)/h(i))^1.3 * (t_var/r(i))^0.3) + min(0.191 * (p(i)/E) * (r(i)/t_var)^2, 0.229)) * E * t_var / r(i));
        options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-6);
        tBuckling = fsolve(bucklingEq, t(i), options);

        t(i) = max(t(i), tBuckling);

        % Hydrostatic pressure
        pHydro = rhoOX .* nx .* g0 .* h(i);

        % Area of the cylinder
        A = pi * (r(i)^2 - (r(i)-t(i)).^2);
        
        % Inertia (valid for a cylinder)
        I = pi / 4 * (r(i)^4 - (r(i) - t(i)).^4);
        
        % Mass computation
        volume = pi .* h(i) .* (r(i)^2 - (r(i) - t(i)).^2);
        mStruct(i) = volume .* rhoMaterial;
        
        % Buckling for Pressurized Cylinders
        k0 = 9 * (t(i) / r(i)).^(0.6) + 0.16 .* (r(i) / h(i)).^(1.3) .* (t(i) / r(i)).^(0.3);
        kp = 0.191 * (p(i) / E) .* (r(i) / t(i)).^2;
        
        sigmaBuckling = ((k0 + min(kp, 0.229)) * E .* t(i)) / r(i);
        
        % if sigmaBuckling < ultimateStress
        %     error('Spessore troppo sottile')
        % end

        % Stresses (MAGNITUDE)
        longitudinalStress = N(i) / A * SF; 
        bendingStress = M(i) ./ I .* r(i) * SF;
        shearStress = T(i) ./ A * SF;
        hoopStress = (p(i) .* r(i) + pHydro .* r(i)) ./ t(i);
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
        tAxial = abs(- N(i) / (2 * pi * r(i)) - M(i) / (pi * r(i)^2)) / ultimateStress * SF;
        tShear = T(i) / (2 * pi * r(i) * shearAllowable) * SF;
    
        t(i) = max(tAxial, tShear);
    
        bucklingEq = @(t_var) ((N(i) / (pi*(r(i)^2-(r(i)-t_var)^2)) + M(i) / (pi/4*(r(i)^4-(r(i)-t_var)^4)) * r(i))) - (9 * (t_var ./ r(i))^1.6 + 0.16 * (t_var ./ h(i))^1.3) * E * t_var / r(i);
            options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-6);
            tBuckling = fsolve(bucklingEq, t(i), options);
            
            t(i) = max(t(i), tBuckling);

        % Area 
        A = pi * (r(i)^2 - (r(i)-t(i)).^2);
        
        % Inertia
        I = pi / 4 * (r(i)^4 - (r(i) - t(i)).^4);
        
        % Stresses (MAGNITUDE)
        longitudinalStress = N(i) / A * SF; 
        bendingStress = M(i) ./ I .* r(i) * SF;
        shearStress = T(i) ./ A * SF;
        
        % Buckling for UnPressurized Cylinders
        
        sigmaBuckling = (9 * (t(i) ./ r(i))^1.6 + 0.16 * (t(i) ./ h(i))^1.3) * E * t(i) / r(i);
        
        stressMatrix(i, :) = [longitudinalStress, bendingStress, shearStress, 0, 0, sigmaBuckling];
    
        % Exctraction of the most critical result
        volume = pi .* h(i) * (r(i)^2 - (r(i) - t(i)).^2);
        mStruct(i) = volume .* rhoMaterial;
    end
end

% =========================== ESTRAZIONE ==================================

mission.structure.mStruct = mStruct;
mission.structure.thickness = t;
mission.structure.stressMatrix = stressMatrix;

end

