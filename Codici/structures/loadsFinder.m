function loads = loadsFinder(nComponents, x, launcher)
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
a       = launcher.acceleration; % [nComponents x 1]
alpha   = launcher.alpha;        % [nComponents x 1]
theta   = launcher.theta;        % scalare (spinta globale)
thrust  = launcher.thrust;       % scalare
lift    = launcher.lift;         % [nComponents x 1]
g0      = launcher.g0;

% ===================== Creation of A Matrix ==============================

% --- Create matrix A (initialized to zero)
A = zeros(3 * nComponents, 3 * nComponents);

% Helper functions to map (component i, type N/T/M) to column index in x
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
    A(rN, colN(i)) = A(rN, colN(i)) - 1;   % -Ni
    
    if i < nComponents
        A(rN, colN(i+1)) = A(rN, colN(i+1)) + 1;   % +N(i+1)
    end
    
    % ==============================
    % TRANSVERSE EQUILIBRIUM (T)
    %   +Ti - T(i+1) = ...
    % ==============================
    A(rT, colT(i)) = A(rT, colT(i)) + 1;   % +Ti
    
    if i < nComponents
        A(rT, colT(i+1)) = A(rT, colT(i+1)) - 1;   % -T(i+1)
    end
    
    % ============================
    % MOMENT EQUILIBRIUM (M)
    %   Mi - M(i+1) - T(i+1)*x(i) = ...
    % dove x(i) è la lunghezza / leva del componente i
    % ============================
    A(rM, colM(i)) = A(rM, colM(i)) + 1;   % +Mi
    
    if i < nComponents
        % -M(i+1)
        A(rM, colM(i+1)) = A(rM, colM(i+1)) - 1;
        
        % -T(i+1)*x(i)  (moment of shear about the section)
        A(rM, colT(i+1)) = A(rM, colT(i+1)) - x(i);
    end
end

% ==================== CREATION OF LOAD VECTOR b ==========================

b  = zeros(3 * nComponents, 1);  % initialize RHS vector

for i = 1:nComponents
    
    % Row indices for component i (coerenti con A)
    rN = 3*(i-1) + 1;   % axial equilibrium
    rT = 3*(i-1) + 2;   % transverse equilibrium
    rM = 3*(i-1) + 3;   % moment equilibrium
    
    % --- Axial equilibrium RHS ---
    % b_N(i) = -drag_i - m_i*g0*sin(alpha_i) - m_i*a_i
    b(rN) = - drag(i) ...
            - m(i)*g0*sin(alpha) ...
            - m(i)*a;
    
    % --- Transverse equilibrium RHS ---
    % b_T(i) = -lift_i - m_i*g0*cos(alpha_i)
    b(rT) = - lift(i) ...
            - m(i)*g0*cos(alpha);
    
    % --- Moment equilibrium RHS ---
    % b_M(i) = (m_i*g0*cos(alpha_i) - lift_i)*x(i)/2
    b(rM) = ( m(i)*g0*cos(alpha) - lift(i) ) * x(i) / 2;
end

% Thrust contribution of the first element
b(1) = b(1) + thrust * cos(theta);
b(2) = b(2) + thrust * sin(theta);

% =========================== SOLUTION ====================================
loads = A \ b;

end