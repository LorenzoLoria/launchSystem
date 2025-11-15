function A = loadsFinder(launcher, nComponents)

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
A = zeros(3 * nComponents, 3 * nComponents);

for i = 1:nComponents
    % Longitudinal Equilibrium
    A(i, i) = 1; % Longitudinal load on component i
    if i < nComponents
        A(i, i + 1) = -1; % Longitudinal load transmission to the next component
    end
    
    % Lateral Equilibrium 
    A(nComponents + i, nComponents + i) = 1; % Lateral load on component i
    if i < nComponents
        A(nComponents + i, nComponents + i + 1) = -1; % Longitudinal load transmission to the next component
    end
    
    % Bending Moment equilibrium
    A(2 * nComponents + i, 2 * nComponents + i) = 1; % bending moment on component i
    if i < nComponents
        A(2 * nComponents + i, 2 * nComponents + i + 1) = -1; % Bending moment transmission to the next component
    end
end

% --- Continuity conditions
for i = 1:nComponents - 1
    % Continuity in longitudinal direction
    A(i, i + 1) = -1;
    A(i + 1, i) = 1;
    
    % Continuity in lateral direction
    A(nComponents + i, nComponents + i + 1) = -1;
    A(nComponents + i + 1, nComponents + i) = 1;
    
    % Continuity in bending moments
    A(2 * nComponents + i, 2 * nComponents + i + 1) = -1;
    A(2 * nComponents + i + 1, 2 * nComponents + i) = 1;
end

end