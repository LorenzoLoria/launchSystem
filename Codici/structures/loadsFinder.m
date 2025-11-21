function loads = loadsFinder(nComponents, launcher)
% Builds matrix A required to solve the linear system Ax = b where:
% A = matrix that encodes the equilibrium equations for
%     axial (N), transverse (T) and bending (M) internal loads
% x = vector containing the internal actions unknown
% b = vector containing all the forces acting on a particular element
%     (size nComponents x 1)

% --- INPUTS:
% nComponents = number of components in which the launcher is discretized
% x = vector of component lengths or lever arms (size nComponents x 1)
% launcher = struct containing required data

% --- OUTPUT:
% loads = [N1 T1 M1 N2 T2 M2 ... Nn Tn Mn]'  (3*nComponents x 1). Solution
% of the linear system

% ======================== DATA CONVERSION ================================
m       = launcher.mass;         % [nComponents x 1]
drag    = launcher.drag;         % [nComponents x 1]
aA      = launcher.accelerationAxial; % scalar (global acceleration)
aN      = launcher.accelerationNormal;
alpha   = launcher.alpha;        % scalar angle
theta   = launcher.theta;        % scalare (spinta globale)
thrust  = launcher.thrust;       % scalare
lift    = launcher.lift;         % [nComponents x 1]
g0      = launcher.g0;
h       = launcher.stagesDimensions; % [nComponents x 1]

% Numero di nodi (dal basso all'alto: 1 ... nComponents+1)
nNodes = nComponents + 1;

% ===================== Creation of A Matrix ==============================

% --- Create matrix A (initialized to zero)
% Now we have N,T,M for all nodes 1..nNodes
A = zeros(3 * nNodes, 3 * nNodes);

% Helper functions to map (node i, type N/T/M) to column index in x
colN = @(i) 3*(i-1) + 1;  % column of Ni
colT = @(i) 3*(i-1) + 2;  % column of Ti
colM = @(i) 3*(i-1) + 3;  % column of Mi

% --- Fill matrix A with equilibrium conditions for each component
for i = 1:nComponents
    
    % Row indices for component i
    rN = 3*(i-1) + 1;   % axial equilibrium
    rT = 3*(i-1) + 2;   % transverse equilibrium
    rM = 3*(i-1) + 3;   % moment equilibrium
    
    % ============================
    % AXIAL EQUILIBRIUM (N)
    %   -Ni + N(i+1) = ...
    % ============================
    A(rN, colN(i))     = A(rN, colN(i))     - 1;   % -Ni
    A(rN, colN(i+1))   = A(rN, colN(i+1))   + 1;   % +N(i+1)
    
    % ==============================
    % TRANSVERSE EQUILIBRIUM (T)
    %   +Ti - T(i+1) = ...
    % ==============================
    A(rT, colT(i))     = A(rT, colT(i))     + 1;   % +Ti
    A(rT, colT(i+1))   = A(rT, colT(i+1))   - 1;   % -T(i+1)
    
    % ============================
    % MOMENT EQUILIBRIUM (M)
    %   Mi - M(i+1) - T(i+1)*x(i) = ...
    % dove x(i) è la lunghezza / leva del componente i
    % ============================
    A(rM, colM(i))     = A(rM, colM(i))     + 1;   % +Mi
    A(rM, colM(i+1))   = A(rM, colM(i+1))   - 1;   % -M(i+1)
    A(rM, colT(i+1))   = A(rM, colT(i+1))   - h(i);% -T(i+1)*x(i)
end

% ==================== CREATION OF LOAD VECTOR b ==========================

% Now b has 3*nNodes rows (3 per component + 3 for boundary at top node)
b  = zeros(3 * nNodes, 1);  % initialize RHS vector

for i = 1:nComponents
    
    % Row indices for component i (coerenti con A)
    rN = 3*(i-1) + 1;   % axial equilibrium
    rT = 3*(i-1) + 2;   % transverse equilibrium
    rM = 3*(i-1) + 3;   % moment equilibrium
    
    % --- Axial equilibrium RHS ---
    % b_N(i) = -drag_i - m_i*g0*sin(alpha_i) - m_i*a_i
    b(rN) = - drag(i) ...
            - m(i)*g0*sin(alpha) ...
            - m(i)*aA;
    
    % --- Transverse equilibrium RHS ---
    % b_T(i) = -lift_i - m_i*g0*cos(alpha_i)
    b(rT) = - lift(i) ...
            - m(i)*g0*cos(alpha) + m(i) * aN;
    
    % --- Moment equilibrium RHS ---
    % b_M(i) = (m_i*g0*cos(alpha_i) - lift_i)*x(i)/2
    b(rM) = ( m(i)*g0*cos(alpha) - lift(i) - m(i) * aN) * h(i) / 2;
end

% Thrust contribution of the first element (nodo 1 in basso)
b(1) = b(1) + thrust * cos(theta); % contributes to axial equilibrium of element 1
b(2) = b(2) + thrust * sin(theta); % contributes to transverse equilibrium of element 1

% =================== BOUNDARY CONDITIONS AT TOP NODE =====================
% Nodo finale = nNodes = nComponents+1, free end:
% N(nNodes) = 0, T(nNodes) = 0, M(nNodes) = 0

rBC_N = 3*nComponents + 1;
rBC_T = 3*nComponents + 2;
rBC_M = 3*nComponents + 3;

A(rBC_N, colN(nNodes)) = 1;  % N_{n+1} = 0
b(rBC_N)               = 0;

A(rBC_T, colT(nNodes)) = 1;  % T_{n+1} = 0
b(rBC_T)               = 0;

A(rBC_M, colM(nNodes)) = 1;  % M_{n+1} = 0
b(rBC_M)               = 0;

% =========================== SOLUTION ====================================
loads = A \ b;

end
