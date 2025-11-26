function Xcg = computeXCG(lc1, lc2, lco, mc1, mc2, mco, m_dot1, m_dot2, t, t1, t2)
% Computes the center of gravity Xcg of a two-stage launcher with a nose cone
%
% Inputs:
%   lc1, lc2 : length to CG of stage 1 and 2
%   lco      : length to CG of conical nose
%   mc1, mc2, mco : masses of stage 1, stage 2, and nose cone
%   m_dot1, m_dot2 : propellant mass flow rates of stage 1 and 2
%   t        : current time
%   t1, t2   : staging times of stage 1 and stage 2
%
% Output:
%   Xcg      : center of gravity measured from the top of the launcher

% Stage 1
if t <= t1
    m1 = mc1 - m_dot1*t;
    X1 = lc1;
    m2 = mc2;
    X2 = lc2;
elseif t > t1 && t <= t2
    m1 = 0;
    X1 = 0;
    m2 = mc2 - m_dot2*(t - t1); % remaining mass of stage 2
    X2 = lc2;
else
    m1 = 0;
    X1 = 0;
    m2 = 0;
    X2 = 0;
end

% CG calculation
Xcg = (m1*X1 + m2*X2 + mco*lco) / (m1 + m2 + mco);

end
