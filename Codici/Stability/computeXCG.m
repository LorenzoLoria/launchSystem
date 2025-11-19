function Xcg = computeXCG(lc1, lc2, lf, lc3, lco, R, rf, mc1, mc2, mf, mc3, mco, m)

% Computes the center of gravity Xcg of a three stage launcher
%
% Inputs:
%   lc1 : length of section 1
%   lc2 : length of section 2
%   lf  : length of fin/flared section
%   lc3 : length of cylindrical section 3
%   lco : length of conical/ogive nose
%   R   : base radius of flare
%   rf  : tip radius of flare
%
%   mc1, mc2, mf, mc3, mco : masses of each section
%   m   : total mass (sum of all components)
%
% Output:
%   Xcg : center of gravity measured from the top of the launcher

    % --- Compute section CG locations (hc1, hc2, hf, hc3, hco) ---
    hc1 = lc1 / 2;

    hc2 = lc1 + lc2/2;

    % Fin/flared section CG (given formula)
    hf = lc1 + lc2 + lf/4 * ( (R^2 + 2*R*rf + 3*rf^2) / (R^2 + R*rf + rf^2) );

    hc3 = lc1 + lc2 + lf + lc3/2;

    hco = lc1 + lc2 + lf + lc3 + lco/4;

    % --- Solve CG equation:  m*(l - Xcg) = sum(m_i * h_i)  ---
    numerator = mco*hco + mc3*hc3 + mf*hf + mc2*hc2 + mc1*hc1;

    % Total length l of vehicle
    l = lc1 + lc2 + lf + lc3 + lco;

    % Compute Xcg
    Xcg = l - numerator / m;

end
