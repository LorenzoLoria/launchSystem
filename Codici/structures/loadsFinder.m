function [N, T, M, A] = loadsFinder(nComponents, launcher)

% Builds matrix A required to solve the linear system Ax = b where:
% A = matrix that encodes the equilibrium equations
% x = vector containing the internal actions unknown [N1 T1 M1 ... Nn Tn Mn]
% b = vector containing external forces
% 
% INPUTS:
% nComponents: number of components in which the LV is discretized;
% mission: struct containing the data needed
% ======================== DATA CONVERSION ================================
m       = mission.launcher.mass;    
mPay    = m(1); % payload mass
m2      = m(2); % second stage
m1      = m(3); % first  stage
drag    = launcher.drag;         
aN      = launcher.accelerationAxial; 
aT      = launcher.accelerationNormal;
gamma   = launcher.gamma;        
alpha   = launcher.alpha;      
lift    = launcher.lift;         
g0      = launcher.g0;
h       = launcher.elementLength; % length of the element
nNodes  = nComponents + 1;

% ===================== Creation of A Matrix ==============================
A = zeros(3 * nNodes, 3 * nNodes);

idxN = @(k) 3*(k-1) + 1;
idxT = @(k) 3*(k-1) + 2;
idxM = @(k) 3*(k-1) + 3;

A(idxN(1), idxN(1)) = 1; % 1*N1 = ... (0 nel vettore b)
A(idxT(1), idxT(1)) = 1; % 1*T1 = ... (0 nel vettore b)
A(idxM(1), idxM(1)) = 1; % 1*M1 = ... (0 nel vettore b)

for i = 1:nComponents
    
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
    A(rM, cT_down)   = h(i);   % T(i+1) * lunghezza
    
end

% ==================== CREATION OF LOAD VECTOR b ==========================
b = zeros(3 * nNodes, 1);

% Node 1
b(1) = 0;
b(2) = 0;
b(3) = 0;

% Node 2
b(4) = - drag * cos(alpha) - lift * sin(alpha);
b(5) = - drag * sin(alpha) + lift * cos(alpha);
b(6) = 0;

% Node 3
b(7) = - mPay * (g0 * sin(gamma) + aN);
b(8) = - mPay * (g0 * cos(gamma) + aT);
b(9) = 0;

% Node 4
b(10) = 0;
b(11) = 0;
b(12) = 0;

% Node5 
b(13) = -m2 * (aN + g0 * sin(gamma));
b(14) = -m2 * (aT + g0 * cos(gamma));
b(15) = 0;

% Node 6
b(16) = 0;
b(17) = 0;
b(18) = 0;

% Node 7
b(19) = - m1 * (aN + g0 * sin(gamma));
b(20) = - m1 * (aT + g0 * cos(gamma));
b(21) = 0;

% Node 8
b(22) = - drag * cos(alpha) - lift * sin(alpha);
b(23) = lift * cos(alpha) - drag * sin(alpha);
b(24) = 0; 

% Node 9
b(25) = 0;
b(26) = 0;
b(27) = 0;

% =========================== SOLUTION ====================================
loads = A \ b;

N = loads(1:3:end);
T = loads(2:3:end);
M = loads(3:3:end);

end