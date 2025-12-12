function Xcg = computeXCG(mission, configuration, launcher, mer, maxQData)
% Computes the center of gravity Xcg of a three-stage launcher with a nose cone
% Assumption: Uniform Weight Distribution 
%
% Inputs:
%   N        : number of stages 
%   lc1, lc2 : length of stage 1 and 2
%   lco      : length of conical nose
%   m        : mass of the launcher at time t 
%   mc2      : wet mass of stage 2 
%   mco      : wet mass of payload 
%
% Output:
%   Xcg      : center of gravity measured from the top of the launcher

N = launcher(1);

% Payload Data
% m = massMaxQ;
lco = mission.capsule.height;
mco = mission.capsule.weight;
Xco = lco / 2; 


if N == 3

    lc3 = configuration.geometry.stage{3}.tanksLength;
    lc2 = configuration.geometry.stage{2}.tanksLength;
    lc1 = configuration.geometry.stage{1}.tanksLength-configuration.stage{1}.engine.length;
    li1 = configuration.geometry.stage{1}.interstage.length;
    li2 = configuration.geometry.stage{2}.interstage.length;
    li3 = configuration.geometry.stage{3}.interstage.length;

    m3 = configuration.stage{3}.mStage-mer.stage{3}.interStage; 
    m2 = configuration.stage{2}.mStage-mer.stage{2}.interStage;
    m1 = configuration.stage{1}.mStage-mer.stage{1}.interStage + maxQData.massMaxQ - configuration.totalMass;
    mi1 = mer.stage{1}.interStage;
    mi2 = mer.stage{2}.interStage;
    mi3 = mer.stage{3}.interStage;

    X6 = lco + li3/2;
    X5 = lco + li3 + lc3/2;               
    X4 = lco + li3 + lc3 + li2/2;
    X3 = lco + li3 + lc3 + li2 + lc2 / 2;
    X2 = lco + li3 + lc3 + li2 + lc2 + li1 / 2;
    X1 = lco + li3 + lc3 + li2 + lc2 + li1 + lc1 / 2;  

    momentSum = m1*X1 + mi1*X2 + m2*X3 + mi2*X4 + m3*X5 + mi3*X6;

    Xcg = (momentSum + mco*Xco) / (m1 + m2 + m3 + mi1 + mi2 + mi3 + mco);

elseif N == 2

    lc2 = configuration.geometry.stage{2}.tanksLength;
    lc1 = configuration.geometry.stage{1}.tanksLength;
    li1 = configuration.geometry.stage{1}.interstage.length;
    li2 = configuration.geometry.stage{2}.interstage.length;

    m2 = configuration.stage{2}.mStage-mer.stage{2}.interStage;
    m1 = configuration.stage{1}.mStage-mer.stage{1}.interStage + maxQData.massMaxQ - configuration.totalMass;
    mi1 = mer.stage{1}.interStage;
    mi2 = mer.stage{2}.interStage;

    X4 = lco + li2/2;
    X3 = lco + li2 + lc2 / 2;
    X2 = lco + li2 + lc2 + li1 / 2;
    X1 = lco + li2 + lc2 + li1 + lc1 / 2;

    momentSum = m1*X1 + mi1*X2 + m2*X3 + mi2*X4;

    Xcg = (momentSum + mco*Xco) / (m1 + m2 + mi1 + mi2 + mco);

elseif N == 1

    lc1 = configuration.geometry.stage{1}.tanksLength-configuration.stage{1}.engine.length;
    li1 = configuration.geometry.stage{1}.interstage.length;

    m1 = configuration.stage{1}.mStage-mer.stage{1}.interStage + maxQData.massMaxQ - configuration.totalMass;
    mi1 = mer.stage{1}.interStage;

    X2 = lco + li1 / 2;
    X1 = lco + li1 + lc1 / 2;

    momentSum = m1*X1 + mi1*X2;

    Xcg = (momentSum + mco*Xco) / (m1 + mi1 + mco);
end

% if N == 2
%     lc2 = configuration.geometry.stage{2}.tankLength;
%     lc1 = configuration.geometry.stage{1}.tankLength;
%     X2 = lco + lc2 / 2;
%     X1 = lco + lc2 + lc1 / 2;
%     m1 = configuration.stage{2}.mStage;
%     m2 = configuration.stage{1}.mStage + mission.structure.massMaxQ - configuration.totalMass;
% 
%     Xcg = (m1*X1 + m2 * X2 + mco * Xco) / mission.structure.massMaxQ;
% end


end