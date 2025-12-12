function [landLoads] = landLoads(mission)

% Function required to compute the loads to which the launcher is subject
% to when on the launchpad

% INPUT

% OUTPUTS
% landload : structure containing N, T, M

% ======================== DATA CONVERSION ================================

% =========================== SOLUTION ====================================

if launcher(1) == 1
    nElements = 8;
elseif launcher(1) == 2
    nElements = 12;
elseif launcher(1) == 3
    nElements = 16;
end

nNodes = nElements + 1;

% Nodes
loadNodes = [2,3:2:nNodes-1,nNodes-1];

xcp = computeXcp(mission, configuration,launcher);
xcp_a = mission.aerodynamics.finsGeom.rootChord - computeFinXcp(mission, maxQData);

% Computation of centroids
centroidCapsule = mission.capsule.height / 3;

% Length of the element used for structural analysis
h = [centroidCapsule];

for ii = launcher(1):-1:1
    if ii == 1
        h = [ h, configuration.geometry.stage{ii}.interstage.length/2, configuration.geometry.stage{ii}.interstage.length/2, (configuration.geometry.stage{ii}.tanksLength-configuration.stage{1}.engine.length)/2,(configuration.geometry.stage{ii}.tanksLength-configuration.stage{1}.engine.length)/2];
    else
        h = [ h, configuration.geometry.stage{ii}.interstage.length/2, configuration.geometry.stage{ii}.interstage.length/2, configuration.geometry.stage{ii}.tanksLength/2,configuration.geometry.stage{ii}.tanksLength/2];
    end
end

h(end) = (configuration.geometry.stage{1}.tanksLength-configuration.stage{1}.engine.length)/2;

% Computation of side area
sideArea = [mission.capsule.radius*mission.capsule.height];

for ii = launcher(1):-1:1
    if ii == 1
        sideArea = [sideArea, 2*configuration.geometry.stage{ii}.interstage.length*configuration.geometry.stage{ii}.interstage.length, 2*(configuration.geometry.stage{ii}.tanksLength-configuration.geometry.stage{ii}.engine.lenght)*configuration.geometry.stage{ii}.radius];
    else
        sideArea = [sideArea, 2*configuration.geometry.stage{ii}.interstage.length*configuration.geometry.stage{ii}.interstage.length, 2*configuration.geometry.stage{ii}.tanksLength*configuration.geometry.stage{ii}.radius];
    end
end

m = [mission.capsule.weight];

for ii = launcher(1):-1:1
    m = [m, mer.stage{ii}.interstage, configuration.geometry.stage{ii}.mStage-mer.stage{ii}.interstage];
end


% Gravity
gN      = 9.81;

% Ground Wind Model
vss = 9.5 * h.^0.2; % [m/s]
effectiveWindSpeed = sqrt((1.25 * vss)^2 + (2.56 * vss)^2);
effWindDynPressure = 0.5 * 1.29 * (2.85 * vss)^2;

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

k = 1;
for i = loadNodes(2:end-1)
    b(:,i) = m(k) * [gN; 0; 0];
    k = k + 1;
end

b(:, [2 3]) = b(:, [3 2]); % inversione richiesta siccome il CG del payload è prima del CP del corpo

b = b(:);

% =========================== SOLUTION ====================================
loads = A \ b;

landLoads.N = loads(1:3:end);
landLoads.T = loads(2:3:end);
landLoads.M = loads(3:3:end);

end