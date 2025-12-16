close all
clear
clc
% --- Full script start ---
rodLength = 10; 

% --- 1. Definizione Forze (Input Generalizzato) ---
forcePositions = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 10]; 
axialForces =    [+10; +10; +10; +10; +10; -10; -10; -10; -10; -10; 0]; 
transversalForces = [+1; -1; +1; -1; +1; -1; +1; -1; -4; +4; 0];
momentCouples = zeros(size(forcePositions)); % Aggiunta di una coppia, se necessario

% Assicurati che tutti i vettori delle forze abbiano la stessa lunghezza
if length(forcePositions) ~= length(axialForces) || length(forcePositions) ~= length(transversalForces)
    error('Input vectors for forces must have the same length.');
end

% --- 2. Global Equilibrium Check ---
sumAxialForces = sum(axialForces);
sumTransversalForces = sum(transversalForces);
sumMoment = sum(momentCouples) + sum(transversalForces .* forcePositions); 

fprintf('--- Equilibrium Check ---\n');
fprintf('Sum of Axial Forces (ΣFx): %.2f\n', sumAxialForces);
fprintf('Sum of Transversal Forces (ΣFy): %.2f\n', sumTransversalForces);
fprintf('Sum of Moments with respect to x=0 (ΣM): %.2f\n', sumMoment);
tolerance = 1e-9;
if (abs(sumAxialForces) > tolerance) || (abs(sumTransversalForces) > tolerance) || (abs(sumMoment) > tolerance)
    warning('ATTENTION: The rod is NOT in equilibrium. The calculated diagrams are correct for the input forces, but the system is NOT self-stable as required.');
end
fprintf('------------------------------\n');

% --- 3. Internal Actions Calculation ---
epsilon = 1e-6; 
sectionPositions = [0; rodLength]; 
for i = 1:length(forcePositions)
    if forcePositions(i) > 0 
        sectionPositions = [sectionPositions; forcePositions(i) - epsilon];
    end
    if forcePositions(i) < rodLength 
        sectionPositions = [sectionPositions; forcePositions(i) + epsilon];
    end
end
sectionPositions = unique(sort(sectionPositions));

normalForces = zeros(size(sectionPositions));
shearForces = zeros(size(sectionPositions));
bendingMoments = zeros(size(sectionPositions));

for j = 1:length(sectionPositions)
    x = sectionPositions(j);
    
    leftIndices = forcePositions < x;
    
    normalForces(j) = sum(axialForces(leftIndices)); 
    shearForces(j) = sum(transversalForces(leftIndices)); 
    
    momentsFromForces = sum(transversalForces(leftIndices) .* (x - forcePositions(leftIndices)));
    momentsFromCouples = sum(momentCouples(leftIndices));
    bendingMoments(j) = momentsFromForces + momentsFromCouples;
end

% --- 4. Plotting the Diagrams ---
figure('Name', 'Internal Actions Diagrams');
subplot(3, 1, 1);
plot(sectionPositions, normalForces, 'b-o', 'LineWidth', 2, 'MarkerSize', 4);
grid on;
title('Normal Force Diagram N(x) (Compression < 0)');
xlabel('Position x');
ylabel('N');
hold on;
plot([0 rodLength], [0 0], 'k--');
ylim([min(normalForces)-1 max(normalForces)+1]);

subplot(3, 1, 2);
plot(sectionPositions, shearForces, 'r-o', 'LineWidth', 2, 'MarkerSize', 4);
grid on;
title('Shear Force Diagram T(x) (Positive if to the left ↑)');
xlabel('Position x');
ylabel('T');
hold on;
plot([0 rodLength], [0 0], 'k--');
ylim([min(shearForces)-1 max(shearForces)+1]);

subplot(3, 1, 3);
plot(sectionPositions, bendingMoments, 'g-o', 'LineWidth', 2, 'MarkerSize', 4);
grid on;
title('Bending Moment Diagram M(x) (Bottom fibers tension > 0)');
xlabel('Position x');
ylabel('M');
hold on;
plot([0 rodLength], [0 0], 'k--');
ylim([min(bendingMoments)-1 max(bendingMoments)+1]);
sgtitle(['Internal Actions for Rod of Length L = ', num2str(rodLength)], 'FontWeight', 'bold');
% --- Full script end ---