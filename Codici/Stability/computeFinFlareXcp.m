function Xcp = computeFinFlareXcp(mission,opt, hf, db)
% computeFinFlareXcp   Computes the new center of pressure (xCP) for a 
%                 fin-stabilized launch vehicle with a flare derived from
%                 Sforza
%
% Inputs:
%   N   : number of stages [m/s]
%   lc2, lc1   : length of stage 2 and stage 1 [m/s]
%   lco   : length of payload [m]
%   alpha   : angle of attack [deg]
%   d   : reference diameter [m]
%   hf   : length of flare [m]
%   db   : diameter of base of flare [m]
%
% Output:
%   xCP : center of pressure of launch vehicle from nose [m]

N = opt.nStages;
alpha = mission.alpha;
lco = mission.capsule.height;
lc1 = opt.stage{1}.length;
lc2 = opt.stage{2}.length;
d = mission.structure.diameter;
hf = 0;
db = 0;

S  = pi*d^2/4;   % reference surface area [m^2]
Sb = pi*db^2/4;  % base surface area for volume normalization [m^2]
dm = (db + d)/2; 

Kf = 1/2;

if N == 2

    l = lco + lc1 + lc2;
    xCP_over_d = ((2/3)*(lc2/d)*(S/Sb)+ (Kf + 1 - S/Sb)*(l/d) - (hf/d)*( (dm^2)/(d^2) - 1 )*(S/Sb)) / (1 + Kf);
    Xcp = xCP_over_d * d;

elseif N == 1 

    Xcp_over_l = 0.63*(1-(sin(alpha))^2)+ 0.5*(lco+lc2)/lco*(sin(alpha))^2;
    Xcp = Xcp_over_l*(lco+lc2);
    
else 
    Xcp = lco * 2/3; 
end


   
   
end
