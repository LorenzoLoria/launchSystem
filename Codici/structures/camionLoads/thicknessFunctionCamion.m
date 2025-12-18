function [updatedStructuralMass, mStruct, tVec, idx] = thicknessFunctionCamion(mission, launcher, configuration, landLoads)

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

end

g0             = mission.environment.g0; 
SF             = mission.structure.safetyFactor; % safety factor

ultimateStress = mission.structure.ultimate; % ultimate stress for chosen material [Pa]
E              = mission.structure.E; % Young Modulus for chosen material [Pa]
shearAllowable = ultimateStress / 2; % ultimate shear stress for chosen material [Pa]
rhoMaterial    = mission.structure.rho; % density of the chosen material [kg/m^3] 
materialNumber = length(E);



N              = landLoads.N; % Axial Load [N] (from loadFinder) (VECTOR)
T              = landLoads.T; % Shear Load [N] (from loadFinder) (VECTOR)
M              = landLoads.M; % Bending Moment [Nm] (from loadFinder) (VECTOR)

% ============================ SOLUTION ===================================

% --- Associate the load to the correct component: the LV is divided into
% payload, interstage between payload and II stage, II stage (pressurized),
% interstage (II and I stage) and I stage (pressurized).

% --- Re-Definition of the loads

Nnew = N([2:2:18]);
Tnew = T([2:2:18]);
Mnew = M([2:2:16]);

Mtail = M([16:1:19]);
[~, idx] = max(abs(Mtail));
firstStageM = Mtail(idx);
Mnew(end+1) = firstStageM;

N = Nnew;
T = Tnew;
M = Mnew;

partsNumber = length(N);

t = zeros(partsNumber, materialNumber);
tVec = zeros(partsNumber, 1);
mStruct = zeros(partsNumber, 1);
stressMatrix = zeros(partsNumber, 6);
idx = zeros(partsNumber, 1);

for i = 1 : partsNumber
    for j = 1 : materialNumber
    % Minimum Allowable Thicknesses
    tAxial = ( abs(N(i) / (2 * pi * r(i))) + abs(M(i) / (pi * r(i)^2)) ) / ultimateStress(j) * SF;
    tShear = abs(T(i)) / (2 * pi * r(i) * shearAllowable(j)) * SF;

    t(i, j) = max(tAxial, tShear);

    bucklingEq = @(t_var) ((abs(N(i)) / (pi*(r(i)^2-(r(i)-t_var)^2)) + abs(M(i)) / (pi/4*(r(i)^4-(r(i)-t_var)^4)) * r(i))) - E(j) * ( 9 * (t_var / r(i))^1.6 + 0.16 * (t_var / h(i))^1.3 );
        options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-6);
        tBuckling = fsolve(bucklingEq, t(i, j), options);
        
        t(i, j) = max(t(i, j), tBuckling);
    end
    
    [tVec(i), idx(i)] = min(t(i,:));
    % [tVec(i), idx(i)] = min( t(i, :) .* (t(i,:)>=mission.structure.minThickness) + mission.structure.minThickness .* (t(i,:)<mission.structure.minThickness) );

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
    if i == 3 && i == 4 && i == 7 && i == 8
        volume = pi .* h(i) .* (r(i)^2 - (r(i) - tVec(i)).^2) + 2*pi*(r(i)^3 - (r(i)-tVec(i))^3);
    else
        volume = pi .* h(i) * (r(i)^2 - (r(i) - tVec(i)).^2);
    end
    mStruct(i) = volume .* rhoMaterial(idx(i));
end

% =========================== ESTRAZIONE ==================================
updatedStructuralMass = sum(mStruct) ;

% Extract the results associated to interstage2, fuel2, ox2, interstage1,
% fuel1, ox1
mStruct = mStruct([2, 3, 4, 6, 7, 8]);
tVec = tVec([2, 3, 4, 6, 7, 8]);
idx = idx([2, 3, 4, 6, 7, 8]);
end