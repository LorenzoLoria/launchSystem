function Cp_max = cpmax_modified_newtonian(M1, gamma)
% cpmax_modified_newtonian
%   Computes C_{p,max} for the modified Newtonian model using
%   the stagnation pressure behind a normal shock (Rayleigh-Pitot formula),
%   followed by isentropic compression to stagnation.
%
%   C_{p,max} = (2 / (gamma M_inf^2)) * (p0/p_inf - 1)

    % --- Normal shock relations (state 1 -> 2) ---
    % Static pressure ratio p2/p1
    p2_p1 = 1 + 2*gamma/(gamma+1) * (M1^2 - 1);

    % Downstream Mach number M2
    num = 1 + 0.5*(gamma-1)*M1^2;
    den = gamma*M1^2 - 0.5*(gamma-1);
    M2_sq = num / den;
    M2 = sqrt(M2_sq);

    % --- Isentropic compression from state 2 to stagnation (2 -> 0_2) ---
    p02_p2 = (1 + 0.5*(gamma-1)*M2^2)^(gamma/(gamma-1));

    % Total stagnation-to-freestream pressure ratio p0/p1
    p0_p1 = p2_p1 * p02_p2;

    % Modified Newtonian stagnation Cp_max
    Cp_max = (2 / (gamma * M1^2)) * (p0_p1 - 1);
end
