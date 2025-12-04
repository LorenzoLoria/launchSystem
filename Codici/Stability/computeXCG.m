function Xcg = computeXCG(mission, opt, m)
% Computes the center of gravity Xcg of a two-stage launcher with a nose cone
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
    lco = mission.capsule.height;
    mco = mission.capsule.weigth;
    Xco = lco / 2; 

    
    if N == 3
        mc3 = opt.stage{3}.engine; 
        mc2 = opt.stage{2}.engine; 
        lc3 = opt.stage{3}.length;
        lc2 = opt.stage{2}.length;
        lc1 = opt.stage{1}.length;
      
        m3 = mc3; 
        m2 = mc2;
        m1 = m - mc3 - mc2 - mco; 

        X3 = lco + lc3/2;               
        X2 = lco + lc3 + lc2/2;         
        X1 = lco + lc3 + lc2 + lc1/2;   
        
        momentSum = (m1*X1) + (m2*X2) + (m3*X3);

    elseif N == 2
        mc2 = opt.stage{2}.engine; 
        lc2 = opt.stage{2}.length;
        lc1 = opt.stage{1}.length;
        
      
        m2 = mc2;
        m1 = m - m2 - mco;
        
        X2 = lco + lc2/2;
        X1 = lco + lc2 + lc1/2;
        
        momentSum = (m1*X1) + (m2*X2);

    elseif N == 1
    
        lc1 = opt.stage{1}.length;
        
        m1 = m - mco;
        
        X1 = lco + lc1/2;
        
        momentSum = (m1*X1);
        
    end

    Xcg = (momentSum + mco*Xco) / m;

end