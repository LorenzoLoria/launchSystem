function [camionLoads] = camionLoads(mission, configuration, mer, launcher, staging)

% Function required to compute the loads to which the launcher is subject
% to when on the launchpad

% INPUT

% OUTPUTS
% landload : structure containing N, T, M

% =========================== SOLUTION ====================================

if launcher(1) == 1
    nElements = 14;
elseif launcher(1) == 2
    nElements = 18;
elseif launcher(1) == 3
    nElements = 22;
end

nNodes = nElements + 1;

% Nodes
loadNodes = [2:2:nNodes];

% Length of the element used for structural analysis
h = [2 / 3 * mission.capsule.height, 1 /3 * mission.capsule.height]; % centroid of cone

for ii = launcher(1):-1:1
    h = [ h, configuration.geometry.stage{ii}.interstage.length/2, ...
    configuration.geometry.stage{ii}.interstage.length/2, configuration.stage{ii}.fuelTankH/2, ...
    configuration.stage{ii}.fuelTankH/2, configuration.stage{ii}.oxTankH/2, ...
    configuration.stage{ii}.oxTankH/2, configuration.geometry.stage{ii}.thrustFrame/2,...
    configuration.geometry.stage{ii}.thrustFrame/2];
end

% Mass computation
m = [mission.capsule.weight];

for ii = launcher(1):-1:1
    m = [m, mer.stage{ii}.interStage, ...
        mer.stage{ii}.tankMassFuel+mer.stage{ii}.cryoInsuFuel, ...
        mer.stage{ii}.tankMassOx+mer.stage{ii}.cryoInsuOx, ...
        mer.stage{ii}.thrustFrame + mer.stage{ii}.engineWeight + ...
        mer.stage{ii}.avionics + mer.stage{ii}.wiring + mer.stage{ii}.tvc +  ...
        mer.stage{ii}.battery + mer.stage{ii}.pressurant];
end

% Gravity
gT      = 9.81;

% Number of nodes
nNodes    = nElements + 1;

% Posizione del centro di forza 

hCUMSUM = cumsum(h);
hCentroids = hCUMSUM(loadNodes - 1);

xcg = centerOfGravity(m, hCentroids);

% Posizioni bracci meccanici rispetto al naso
xc1 = hCUMSUM(5); % scelta a caso
xc1 = xcg - xc1;
xc2 = hCUMSUM(12); % scelta a caso
xc2 = xc2 - xcg;

% Bilancio forze
T1 = sum(m)*gT / (1 + xc1 / xc2);
T2 = T1 * xc1 / xc2;

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

k = 1;
for i = loadNodes
    b(:,i) = [0; m(k)*gT; 0];
    k = k + 1;
end

b(:,6) = b(:, 6) - [0;T1;0];
b(:,13) = b(:, 13) - [0;T2;0];

b = b(:);

% =========================== SOLUTION ====================================
loads = A \ b;

camionLoads.N = loads(1:3:end);
camionLoads.T = loads(2:3:end);
camionLoads.M = loads(3:3:end);

end