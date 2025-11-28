% -------------------------------------------------------------------------
% ---------------------- MONTECARLO ANALYSIS ------------------------------
% -------------------------------------------------------------------------
clear all
clc
close all

[mission,optimisation] = dataStruct;

% Definition of the uncertanties

sizeMC = 30;

structuralMassUncertainty = 0.15 * randn(sizeMC,1);
propellantMassUncertainty = 0.05 * randn(sizeMC,1);
specificImpulseUncertainty = 0.01 * randn(sizeMC,1);

% assembly of all possible combination
Matr = [];
for i = 1 :  optimisation.nStages
    structuralMass(:,i) = optimisation.stage{i}.structuralMass + structuralMassUncertainty .* optimisation.stage{i}.structuralMass;
    propellantMass(:,i) = optimisation.stage{i}.mProp + propellantMassUncertainty .* optimisation.stage{i}.mProp;
    specificImpulse(:,i) = optimisation.stage{i}.Isp + specificImpulseUncertainty .* optimisation.stage{i}.Isp;
    [T, A, N] = ndgrid(structuralMass(:,i), propellantMass(:,i), specificImpulse(:,i));
    Matr =[Matr ,T(:),A(:),N(:)];
end

% Shuffle of the Matrix elements
for i=1:size(Matr,2)
    shuffledMatrix(:,i) = Matr(randperm(size(Matr,1)),i);
end

% Start of Montecarlo analysis
parfor i = 1:size(shuffledMatrix,1)
    


%%
Matr1 = Matr(:,1);
Matr2 = Matr(:,2);
Matr3 = Matr(:,3);
% shuffle
Matr1=Matr1(randperm(length(Matr1)));
Matr2=Matr2(randperm(length(Matr2)));
Matr3=Matr3(randperm(length(Matr3)));
Matr=[Matr1 Matr2 Matr3];
