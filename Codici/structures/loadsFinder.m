function [N, T, M] = loadsFinder(nComponents, launcher)
% Builds matrix A required to solve the linear system Ax = b where:
% A = matrix that encodes the equilibrium equations
% x = vector containing the internal actions unknown [N1 T1 M1 ... Nn Tn Mn]
% b = vector containing external forces
%
% CONVENTION CHANGE:
% Node 1 = NOSE (TOP)
% Node nNodes = TAIL (BOTTOM)

% ======================== DATA CONVERSION ================================
m       = launcher.mass;         
drag    = launcher.drag;         
aA      = launcher.accelerationAxial; 
aN      = launcher.accelerationNormal;
gamma   = launcher.gamma;        
theta   = launcher.theta;        
thrust  = launcher.thrust;       
lift    = launcher.lift;         
g0      = launcher.g0;
h       = launcher.stagesDimensions; 
hCP     = launcher.pressureCenterPosition;
nNodes  = nComponents + 1;

% ===================== Creation of A Matrix ==============================
A = zeros(3 * nNodes, 3 * nNodes);

colN = @(i) 3*(i-1) + 1;  
colT = @(i) 3*(i-1) + 2;  
colM = @(i) 3*(i-1) + 3;  

% --- Fill matrix A with equilibrium conditions for each component
% Element i connects Node i (Top) to Node i+1 (Bottom)
for i = 1:nComponents
    
    % Row indices for component i
    rN = 3*(i-1) + 1;   
    rT = 3*(i-1) + 2;   
    rM = 3*(i-1) + 3;   
    
    % ============================
    % AXIAL EQUILIBRIUM (N)
    % -N_top + N_bottom = Loads
    % ============================
    A(rN, colN(i))     = A(rN, colN(i))     + 1;   % +Ni (Top)
    A(rN, colN(i+1))   = A(rN, colN(i+1))   - 1;   % -N(i+1) (Bottom)
    
    % ==============================
    % TRANSVERSE EQUILIBRIUM (T)
    % +T_top - T_bottom = Loads
    % ==============================
    A(rT, colT(i))     = A(rT, colT(i))     - 1;   % -Ti (Top)
    A(rT, colT(i+1))   = A(rT, colT(i+1))   + 1;   % +T(i+1) (Bottom)
    
    % ============================
    % MOMENT EQUILIBRIUM (M)
    % --- MODIFIED ---
    % Calcolando il momento rispetto al nodo inferiore (i+1):
    % M_top - M_bottom - T_top * h = ...
    % (Nota: la struttura dei segni dipende dalla convenzione, qui mantengo
    %  la coerenza con la riga originale: M_i - M_i+1 - TermineTaglio)
    %
    % La differenza cruciale: il braccio h agisce sul taglio al nodo i (TOP),
    % non i+1.
    % ============================
    A(rM, colM(i))     = A(rM, colM(i))     - 1;   % -Mi
    A(rM, colM(i+1))   = A(rM, colM(i+1))   + 1;   % +M(i+1)
    
    % CHANGE: Use colT(i) instead of colT(i+1) because shear at the top (Node i)
    % creates a moment arm relative to the bottom (Node i+1).
    A(rM, colT(i))     = A(rM, colT(i))     + h(i); % -T(i)*h(i)
end

% ==================== CREATION OF LOAD VECTOR b ==========================
b  = zeros(3 * nNodes, 1); 

for i = 1:nComponents
    rN = 3*(i-1) + 1;   
    rT = 3*(i-1) + 2;   
    rM = 3*(i-1) + 3;   
    
    % I segni di b rimangono generalmente validi (rappresentano la variazione
    % di carico lungo l'elemento), assumendo che i vettori drag/lift siano
    % definiti coerentemente.
    
    b(rN) = - drag / nComponents - m(i)*g0*sin(gamma) - m(i) * aA;
    
    b(rT) = - lift / nComponents + m(i)*g0*cos(gamma) + m(i) * aN;
    
    b(rM) = ( m(i)*g0*cos(gamma) + m(i) * aN) * h(i) / 2 - lift * hCP(i);
end

% ============================ THRUST =====================================

% rN_last = 3*(nComponents-1) + 1;
% rT_last = 3*(nComponents-1) + 2; 

% b(rN_last) = b(rN_last) + thrust * cos(theta); 
% b(rT_last) = b(rT_last) + thrust * sin(theta); 

% =================== BOUNDARY CONDITIONS (NOSE = NODE 1) =================
% --- MODIFIED ---
% Il nodo 1 (Naso) è ora l'estremità libera.
% Le condizioni al contorno (spostate nelle righe finali della matrice per comodità)
% devono imporre N1=0, T1=0, M1=0.

rBC_N = 3*nComponents + 1;
rBC_T = 3*nComponents + 2;
rBC_M = 3*nComponents + 3;

% Imponiamo condizioni sul nodo 1 (colonne 1, 2, 3)
A(rBC_N, colN(1)) = 1;  % N_1 = 0
b(rBC_N)          = 0;

A(rBC_T, colT(1)) = 1;  % T_1 = 0
b(rBC_T)          = 0;

A(rBC_M, colM(1)) = 1;  % M_1 = 0
b(rBC_M)          = 0;

% =========================== SOLUTION ====================================
loads = A \ b;

N = loads(1:3:end);
T = loads(2:3:end);
M = loads(3:3:end);

end