function A = loadsFinder(launcher, nComponents, F_ext)

% Builds matrix A required to solve the linear system Ax = b where:
% A = matrix that allows the definition of the equilibrium equations for
% longitudinal, lateral and bending loads
% x = vector of internal loads required to counter-balance the effects of
% the external loads
% b = vector of external loads
    
% --- INPUTS
% launcher = struct containing launcher data
% nComponents = number of components in which the LV has been divided

% --- OUTPUTS
% x = vector (3*nComponents, 1) containing the internal loads of the
% structure

% --- SOLUTION ------------------------------------------------------------ 
    
% --- Create matrix A (initialized to zero)
A = zeros(3 * nComponents, 3 * nComponents);

% --- Fill matrix A with equilibrium conditions for each component
for i = 1:nComponents
    % Axial equilibrium (forces along x direction)
    A(i, i) = 1; % Axial force of component i
    if i < nComponents
        A(i, i + 1) = -1; % Transfer of axial force to the next component
    end
    
    % Transversal equilibrium (forces along z direction)
    A(nComponents + i, nComponents + i) = 1; % Transversal force of component i
    if i < nComponents
        A(nComponents + i, nComponents + i + 1) = -1; % Transfer of transversal force
    end
    
    % Moment equilibrium (torques)
    A(2 * nComponents + i, 2 * nComponents + i) = 1; % Moment of component i
    if i < nComponents
        A(2 * nComponents + i, 2 * nComponents + i + 1) = -1; % Transfer of moment
    end
end

% --- Apply continuity conditions for forces and moments
for i = 1:nComponents - 1
    % Continuity in axial direction
    A(i, i + 1) = -1;
    A(i + 1, i) = 1;
    
    % Continuity in transversal direction
    A(nComponents + i, nComponents + i + 1) = -1;
    A(nComponents + i + 1, nComponents + i) = 1;
    
    % Continuity in moments
    A(2 * nComponents + i, 2 * nComponents + i + 1) = -1;
    A(2 * nComponents + i + 1, 2 * nComponents + i) = 1;
end

% --- Solve the system Ax = b
% F_ext is the vector of external forces applied (should be of length 3*nComponents)
% Solve the linear system to find the internal reactions x
x = A \ F_ext;  % Solve the system (A*x = b)

end