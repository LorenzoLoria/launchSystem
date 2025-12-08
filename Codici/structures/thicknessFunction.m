function [mission] = thicknessFunction(mission, engineUsed)

% Function required to size the launcher thickess of all the different 
% components of the LV. Evaluation must be done in the most
% critical conditions: q, q-alpha, MECO 

% --- INPUTS
% mission                         = struct containing LV data
% engineUsed                      = parameter indicating which engine is
%                                   used

% --- OUTPUTS
% mission.structure.mStruct       = vector containing mass of each structure [kg]
% mission.structure.t             = vector containing thickness of each component [m]
% mission.structure.tMax          = maximum thickness of the launcher [m]
% mission.structure.stressMatrix  = (nComponents x 6) matrix. Each column represents a
%                                   component of the LV (e.g. Thrust structure, 1st stage etc.), 
%                                   each row represents a different kind of stress [MPa]


% ============================== DATA =====================================

r              = mission.structure.diameter / 2; % radius of cylindrical section of each body [m] (VECTOR)
h              = mission.structure.componentLength; % length of each body [m] (VECTOR)
a              = mission.structure.aMaxQ;
g0             = mission.environment.g0; 
nx             = a(1) / g0 + 1; 
SF             = mission.structure.SF; % safety factor

ultimateStress = mission.structure.ultimate; % ultimate stress for chosen material [Pa]
E              = mission.structure.E; % Young Modulus for chosen material [Pa]
shearAllowable = ultimateStress / 2; % ultimate shear stress for chosen material [Pa]
rhoMaterial    = mission.structure.rho; % density of the chosen material [kg/m^3] 
materialNumber = length(E);

rhoOX          = mission.launcher.engines{engineUsed}.oxDens; % density of the oxidizer [kg/m^3]
rhoFu          = mission.launcher.engines{engineUsed}.fuelDens; % density of the fuel [kg/m^3]

p              = mission.structure.pressurization; % pressure of component [Pa] [partsNumber x 1]
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

for ii=6:2:length(N)-3
    interN = [interN, N(ii)];
    interT = [interT, T(ii)];
    interM = [interM, M(ii)];
end

firstStageN = N(end);
firstStageT = T(end);
Mtail = M(end-2:end);
[~, idx] = max(abs(Mtail));
firstStageM = Mtail(idx);


N = [payloadN, interN, firstStageN];
T = [payloadT, interT, firstStageT];
M = [payloadM, interM, firstStageM];

partsNumber = length(N);

r = r * ones(partsNumber, 1);

t = zeros(partsNumber, materialNumber);
tVec = zeros(partsNumber, 1);
mStruct = zeros(partsNumber, 1);
stressMatrix = zeros(partsNumber, 6);
idx = zeros(partsNumber, 1);

for i = 1 : partsNumber
    
    if p(i) ~= 0
        for j = 1 : materialNumber
        % Minimum Allowable Thicknesses
        tAxial = abs((- N(i) / (2 * pi * r(i)) - M(i) / (pi * r(i)^2)) * SF + p(i) * r(i) / 2 + rhoFu * nx * g0 * r(i) * h(i) / 2) / ultimateStress(j);
        tShear = T(i) / (2 * pi * r(i) * shearAllowable(j)) * SF;

        t(i, j) = max(tAxial, tShear);

        % Consider buckling
        bucklingEq = @(t_var) abs( ...
                    - abs(N(i) / (pi*(r(i)^2-(r(i)-t_var)^2))) * SF ...
                    - abs(M(i) / (pi/4*(r(i)^4-(r(i)-t_var)^4)) * r(i)) * SF ...
                    + abs((p(i) * r(i)) / (2 * t_var)) ...
                    + abs((rhoFu * nx * g0 * h(i) * r(i)) / (2 * t_var)) ...
                    ) - (((9 * (t_var/r(i))^0.6 + 0.16 * (r(i)/h(i))^1.3 * (t_var/r(i))^0.3) + min(0.191 * (p(i)/E(j)) * (r(i)/t_var)^2, 0.229)) * E(j) * t_var / r(i));

        options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-6);
        tBuckling = fsolve(bucklingEq, t(i), options);
        t(i, j) = max(t(i, j), tBuckling);
        end

        [tVec(i), idx(i)] = min(t(i, :));

        % Hydrostatic pressure
        pHydro = rhoFu .* nx .* g0 .* h(i);

        % Area of the cylinder
        A = pi * (r(i)^2 - (r(i)-tVec(i)).^2);
        
        % Inertia (valid for a cylinder)
        I = pi / 4 * (r(i)^4 - (r(i) - tVec(i)).^4);
        
        % Mass computation
        volume = pi .* h(i) .* (r(i)^2 - (r(i) - tVec(i)).^2);
        mStruct(i) = volume .* rhoMaterial(idx(i));
        
        % Buckling for Pressurized Cylinders
        k0 = 9 * (tVec(i) / r(i)).^(0.6) + 0.16 .* (r(i) / h(i)).^(1.3) .* (tVec(i) / r(i)).^(0.3);
        kp = 0.191 * (p(i) / E(idx(i))) .* (r(i) / tVec(i)).^2;
        
        sigmaBuckling = ((k0 + min(kp, 0.229)) * E(idx(i)) .* tVec(i)) / r(i);
        
        % if sigmaBuckling < ultimateStress
        %     error('Spessore troppo sottile')
        % end

        % Stresses (MAGNITUDE)
        longitudinalStress = N(i) / A; 
        bendingStress = M(i) ./ I .* r(i);
        shearStress = T(i) ./ A;
        hoopStress = (p(i) .* r(i) + pHydro .* r(i)) ./ tVec(i);
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
        
        for j = 1 : materialNumber
        % Minimum Allowable Thicknesses
        tAxial = ( abs(- N(i) / (2 * pi * r(i))) - abs(M(i) / (pi * r(i)^2)) ) / ultimateStress(j) * SF;
        tShear = T(i) / (2 * pi * r(i) * shearAllowable(j)) * SF;
    
        t(i, j) = max(tAxial, tShear);
    
        bucklingEq = @(t_var) ((N(i) / (pi*(r(i)^2-(r(i)-t_var)^2)) + M(i) / (pi/4*(r(i)^4-(r(i)-t_var)^4)) * r(i))) - 0.6 * (1 - 0.901 * (1 - exp(-1 / 16 * sqrt(r(i)/t_var)))) * E(j) * t_var / r(i);
            options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-6);
            tBuckling = fsolve(bucklingEq, t(i, j), options);
            
            t(i, j) = max(t(i, j), tBuckling);
        end

        [tVec(i), idx(i)] = min(t(i, :));

        % Area 
        A = pi * (r(i)^2 - (r(i)-tVec(i)).^2);
        
        % Inertia
        I = pi / 4 * (r(i)^4 - (r(i) - tVec(i)).^4);
        
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

mission.structure.mStruct = mStruct;
mission.structure.thickness = tVec;
mission.structure.stressMatrix = stressMatrix;
mission.structure.materials = idx;

end