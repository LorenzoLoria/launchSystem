function [internalActions] = loadsFinderMECO(mission, launcher, configuration, mecoData)

% Builds matrix A required to solve the linear system Ax = b where:
% A = matrix that encodes the equilibrium equations
% x = vector containing the internal actions unknown [N1 T1 M1 ... Nn Tn Mn]
% b = vector containing external forces
% 
% INPUTS:
% nComponents: number of components in which the LV is discretized;
% mission: struct containing the data needed
% ======================== DATA CONVERSION ================================

if launcher(1) == 1
    nElements = 17;
elseif launcher(1) == 2
    nElements = 21;
elseif launcher(1) == 3
    nElements = 25;
end

loadNodes = [2 3 5 7 9 11 12 14 16 18 20 21];

xcp = computeXcp(mission, configuration,launcher);
xcp_a = mission.aerodynamics.finsGeom.rootChord - computeFinXcp(mission, mecoData);

% Length of the element used for structural analysis

h = [mission.capsule.height/2, xcp - mission.capsule.height/2, mission.capsule.height - xcp];

for ii = launcher(1):-1:1
        h = [ h, configuration.geometry.stage{ii}.interstage.length/2, ...
            configuration.geometry.stage{ii}.interstage.length/2, ...
            configuration.stage{ii}.fuelTankH/2, configuration.stage{ii}.fuelTankH/2, ...
            configuration.stage{ii}.oxTankH/2, configuration.stage{ii}.oxTankH/2, ...
            configuration.geometry.stage{ii}.thrustFrame-xcp_a,xcp_a-configuration.geometry.stage{ii}.thrustFrame/2,...
            configuration.geometry.stage{ii}.thrustFrame/2];
end

m  = mecoData.massMECOVec;   

% Drag
drag    = mecoData.dMECO;
dragN = drag(1);
dragT = drag(2);

% Lift
lift    = mecoData.lMECO;   
liftN   = lift(1);
liftT   = lift(2);

% Gravity
g0      = mecoData.gMECO;
gN      = g0(1);
gT      = g0(2);

% Acceleration
a       = mecoData.aMECO;
aN      = a(1); 
aT      = a(2);    

% Fins' Lift
liftFins  = mecoData.liftFinsMECO;
liftFinsN = liftFins(1);
liftFinsT = liftFins(2);

% Fins' Drag
dragFins  = mecoData.dragFinsMECO;
dragFinsN = dragFins(1);
dragFinsT = dragFins(2);

% Number of nodes
nNodes    = nElements + 1;

% ===================== Creation of A Matrix ==============================
A = zeros(3 * nNodes, 3 * nNodes);

idxN = @(k) 3*(k-1) + 1;
idxT = @(k) 3*(k-1) + 2;
idxM = @(k) 3*(k-1) + 3;

A(idxN(1), idxN(1)) = 1; % 1*N1 = ... (0 nel vettore b)
A(idxT(1), idxT(1)) = 1; % 1*T1 = ... (0 nel vettore b)
A(idxM(1), idxM(1)) = 1; % 1*M1 = ... (0 nel vettore b)

for i = 1:nElements
    
    nodeUp = i;         % Nodo precedente (i)
    nodeDown = i + 1;   % Nodo attuale (i+1)
    
    % Indici delle RIGHE dove scriviamo l'equazione (Nodo a valle)
    rN = idxN(nodeDown);
    rT = idxT(nodeDown);
    rM = idxM(nodeDown);
    
    % Indici delle COLONNE (Variabili del nodo attuale i+1)
    cN_down = idxN(nodeDown); cT_down = idxT(nodeDown); cM_down = idxM(nodeDown);
    
    % Indici delle COLONNE (Variabili del nodo precedente i)
    cN_up = idxN(nodeUp); cT_up = idxT(nodeUp); cM_up = idxM(nodeUp);
    
    % --- EQUILIBRIO NORMALE (N) ---
    % N(i+1) - N(i) = CaricoAssiale
    A(rN, cN_down) =  1; 
    A(rN, cN_up)   = -1;
    
    % --- EQUILIBRIO TAGLIO (T) ---
    % T(i+1) - T(i) = CaricoTrasversale
    A(rT, cT_down) =  1;
    A(rT, cT_up)   = -1;
    
    % --- EQUILIBRIO MOMENTO (M) ---
    % M(i+1) = M(i) - T(i)*LunghezzaElemento + MomentoCarichi
    % M(i+1) - M(i) + T(i)*x = 0
    
    A(rM, cM_down) =  1;   % + M(i+1)
    A(rM, cM_up)   = -1;   % - M(i)
    A(rM, cT_up) = h(i);   % T(i+1) * lunghezza
    
end

% ==================== CREATION OF LOAD VECTOR b ==========================
b = zeros(3, nNodes);
b(:,3) = [dragN+liftN; dragT+liftT; 0];

k = 1;
for i = [2 5 7 9 12 14 16 18 21]
    b(:,i) = m(k) * [gN-aN; gT-aT; 0];
    k = k + 1;
end

b(:, 11) = [0.1*(dragFinsN+liftFinsN); 0.1*(dragFinsT+liftFinsT); 0];
b(:, 20) = [0.9*(dragFinsN+liftFinsN); 0.9*(dragFinsT+liftFinsT); 0];
b = b(:);

% =========================== SOLUTION ====================================
loads = A \ b;

internalActions.N = loads(1:3:end);
internalActions.T = loads(2:3:end);
internalActions.M = loads(3:3:end);
end