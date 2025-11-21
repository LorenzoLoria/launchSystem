function [b,loads] = loadsFinder_freefree(nComponents, launcher)
% Free-free beam launcher internal loads solver
% Each node has [N, T, M]; total 3*(nComponents+1) unknowns.
% We set N1=T1=M1=0 and Nn+1=Tn+1=Mn+1=0 (free ends).

m       = launcher.mass(:);
drag    = launcher.drag(:);
a       = launcher.acceleration;
alpha   = launcher.alpha;
theta   = launcher.theta;
thrust  = launcher.thrust;
lift    = launcher.lift(:);
g0      = launcher.g0;
x       = launcher.stagesDimensions;

% -------------------- Matrix sizes --------------------
nNodes = nComponents + 1;
nUnknowns = 3 * nNodes;        % all internal actions
A = zeros(3*nComponents, nUnknowns);
b = zeros(3*nComponents, 1);

% Helper anonymous maps
colN = @(i) 3*(i-1) + 1;
colT = @(i) 3*(i-1) + 2;
colM = @(i) 3*(i-1) + 3;

% -------------------- Fill equations --------------------
for i = 1:nComponents
    rN = 3*(i-1) + 1;
    rT = 3*(i-1) + 2;
    rM = 3*(i-1) + 3;

    % --- Axial equilibrium ---
    A(rN, colN(i))   = -1;
    A(rN, colN(i+1)) =  1;

    % --- Shear equilibrium ---
    A(rT, colT(i))   =  1;
    A(rT, colT(i+1)) = -1;

    % --- Moment equilibrium ---
    A(rM, colM(i))   =  1;
    A(rM, colM(i+1)) = -1;
    A(rM, colT(i+1)) = -x(i);

    % --- RHS ---
    b(rN) = -drag(i) - m(i)*g0*sin(alpha) - m(i)*a;
    b(rT) = -lift(i) - m(i)*g0*cos(alpha);
    b(rM) = (m(i)*g0*cos(alpha) - lift(i)) * x(i)/2;
end

% -------------------- Apply free-free BCs --------------------
% N1=T1=M1 = 0, Nn=Tn=Mn = 0
knownIdx = [1,2,3, 3*nNodes-2, 3*nNodes-1, 3*nNodes]; % indices of known DOFs
unknownIdx = setdiff(1:nUnknowns, knownIdx);

% Eliminate known DOFs
A_reduced = A(:, unknownIdx);
x_unknown = A_reduced \ b;

% Build full vector (with zeros at boundaries)
loads_full = zeros(nUnknowns,1);
loads_full(unknownIdx) = x_unknown;

loads = loads_full;
end
