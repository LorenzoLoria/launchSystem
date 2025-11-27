function [epsS] = StructuralMassIndex(nStage,propellantMass)

% Function to compute the structural mass index (SFORZA book for values, pp 255)

% INPUTS:
% nStage : stage being considered
% propellantMass : mass of the propellant
% OUTPUT:
% epsS : structural mass index

if nStage == 1
    PropPoints = [37800 39200 96120 101000 112700 117800 155000 156260 158000 174000 177800 226000 284089 419400 71800 2145700];
    epsPoints  = [0.0794 0.0965 0.0591 0.1277 0.0720 0.0569 0.0516 0.0658 0.0772 0.0718 0.1361 0.0646 0.0741 0.0742 0.0493 0.0621];
    [p] = polyfit(PropPoints,epsPoints,1);
    epsS = polyval(p,propellantMass);
elseif nStage == 2
    PropPoints = [6000 9700 16600 16900 16930 20830 28900 34000 35000 86000 95400 156000 456100 21100];
    epsPoints  = [0.1583 0.1227 0.2048 0.1834 0.1087 0.1077 0.0830 0.0966 0.1286 0.0640 0.0721 0.0751 0.0712 0.1315];
    [p] = polyfit(PropPoints,1./epsPoints,1);
    epsS = 1/polyval(p,propellantMass);
elseif nStage == 3
    PropPoints = [10700 22000 46600 106300];
    epsPoints  = [0.1121 0.1070 0.0898 0.1097];
    [p] = polyfit(PropPoints,epsPoints,1);
    epsS = polyval(p,propellantMass);
end


end