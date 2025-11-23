function [epsS] = StructuralMassIndex(nStage,PropellantMass)

% This function estimates the structural mass index

    if nStage == 1
        PropPoints = [20000 39000 40000 95000 100000 120000 140000 150000 160000 165000 170000 180000 190000 220000 250000 290000 400000 700000 700000 700000 2000000];
        epsPoints  = [0.068 0.075 0.09 0.055 0.115 0.07 0.057 0.065 0.05 0.065 0.07 0.044 0.12 0.062 0.05 0.07 0.07 0.05 0.055 0.06 0.06];
        [p] = polyfit(PropPoints,epsPoints,1);
        epsS = polyval(p,PropellantMass);
    elseif nStage == 2
        PropPoints = [17000 18000 19000 20000 20500 28000 35000 36000 44000 85000 90000 160000 450000];
        epsPoints  = [0.172 0.158 0.096 0.097 0.116 0.076 0.088 0.115 0.05 0.061 0.068 0.072 0.068];
        [p] = polyfit(PropPoints,1./epsPoints,1);
        epsS = 1/polyval(p,PropellantMass);
    elseif nStage == 3
        PropPoints = [10000 25000 45000 110000];
        epsPoints  = [0.101 0.096 0.083 0.1];
        [p] = polyfit(PropPoints,epsPoints,1);
        epsS = polyval(p,PropellantMass);
    end


end