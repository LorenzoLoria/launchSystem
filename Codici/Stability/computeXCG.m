function Xcg = computeXCG(mission, opt)
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

N = opt.nStages;

% Payload Data
% m = massMaxQ;
lco = mission.capsule.height;
mco = mission.capsule.weigth;
Xco = lco / 2; 


if N == 3
     
    lc3 = opt.stage{3}.length;
    lc2 = opt.stage{2}.length;
    lc1 = opt.stage{1}.length;
    li1 = mission.structures{1}.lengthInterstage;
    li2 = mission.structures{2}.lengthInterstage;
    li3 = mission.structures{3}.lengthInterstage;
  
    m3 = opt.stage{3}.mStage; 
    m2 = opt.stage{2}.mStage;
    m1 = opt.stage{1}.mStage + mission.structure.massMaxQ - opt.m0Tot;
    mi1 = mission.structures{1}.mInterstage;
    mi2 = mission.structures{2}.mInterstage;
    mi3 = mission.structures{3}.mInterstage;
    
    X6 = lco + li3/2;
    X5 = lco + li3 + lc3/3;               
    X4 = lco + li3 + lc3 + li2/2;
    X3 = lco + li3 + lc3 + li2 + lc2 / 2;
    X2 = lco + li3 + lc3 + li2 + lc2 + li1 / 2;
    X1 = lco + li3 + lc3 + li2 + lc2 + li1 + lc1 / 2;  
    
    momentSum = m1*X1 + mi1*X2 + m2*X3 + mi2*X4 + m3*X5 + mi3*X6;

    Xcg = (momentSum + mco*Xco) / (m1 + m2 + m3 + mi1 + mi2 + mi3 + mco);

elseif N == 2

    lc2 = opt.stage{2}.length;
    lc1 = opt.stage{1}.length;
    li1 = mission.structures{1}.lengthInterstage;
    li2 = mission.structures{2}.lengthInterstage;
  
    m2 = opt.stage{2}.mStage;
    m1 = opt.stage{1}.mStage + mission.structure.massMaxQ - opt.m0Tot;
    mi1 = mission.structures{1}.mInterstage;
    mi2 = mission.structures{2}.mInterstage;
                 
    X4 = lco + li2/2;
    X3 = lco + li2 + lc2 / 2;
    X2 = lco + li2 + lc2 + li1 / 2;
    X1 = lco + li2 + lc2 + li1 + lc1 / 2;
    
    momentSum = m1*X1 + mi1*X2 + m2*X3 + mi2*X4;
    
    Xcg = (momentSum + mco*Xco) / (m1 + m2 + mi1 + mi2 + mco);

elseif N == 1

    lc1 = opt.stage{1}.length;
    li1 = mission.structures{1}.lengthInterstage;
  
    m1 = opt.stage{1}.mStage + mission.structure.massMaxQ - opt.m0Tot;
    mi1 = mission.structures{1}.mInterstage;
                 
    X2 = lco + li1 / 2;
    X1 = lco + li1 + lc1 / 2;
    
    momentSum = m1*X1 + mi1*X2;
    
    Xcg = (momentSum + mco*Xco) / (m1 + mi1 + mco)
end



end