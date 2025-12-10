function [mission, b] = loadsFinder(mission)

% Builds matrix A required to solve the linear system Ax = b where:
% A = matrix that encodes the equilibrium equations
% x = vector containing the internal actions unknown [N1 T1 M1 ... Nn Tn Mn]
% b = vector containing external forces
% 
% INPUTS:
% nComponents: number of components in which the LV is discretized;
% mission: struct containing the data needed
% ======================== DATA CONVERSION ================================

nElements = mission.structure.nElements;
m         = mission.structure.massMaxQVec;   

% Drag
drag    = mission.structure.dMaxQ;
dragN = drag(1);
dragT = drag(2);

% Lift
lift    = mission.structure.lMaxQ;   
liftN   = lift(1);
liftT   = lift(2);

% Gravity
g0      = mission.structure.gMaxQ;
gN      = g0(1);
gT      = g0(2);

% Acceleration
a       = mission.structure.aMaxQ;
aN      = a(1); 
aT      = a(2);    

% Length of the elements
h       = mission.structure.elementLength;

% Fins' Lift
liftFins  = mission.structure.liftFinsMaxQ;
liftFinsN = liftFins(1);
liftFinsT = liftFins(2);

% Fins' Drag
dragFins  = mission.structure.dragFinsMaxQ;
dragFinsN = dragFins(1);
dragFinsT = dragFins(2);

% Number of nodes
nNodes    = nElements + 1;
loadNodes = mission.structure.loadNodes; % indicates the nodes associated to loads

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
b(:,2) = [dragN+liftN; dragT+liftT; 0];

k = 1;
for i = loadNodes(2:end-1)
    b(:,i) = m(k) * [gN-aN; gT-aT; 0];
    k = k + 1;
end

b(:, [2 3]) = b(:, [3 2]); % inversione richiesta siccome il CG del payload è prima del CP del corpo

b(:,end-1) = [dragFinsN+liftFinsN; dragFinsT+liftFinsT; 0];
%b(:,end-1) = [dragFinsN+liftFinsN; -1.119698091469362e+06; 0];
b = b(:);

% =========================== SOLUTION ====================================
loads = A \ b;

mission.structure.N = loads(1:3:end);
mission.structure.T = loads(2:3:end);
mission.structure.M = loads(3:3:end);
end