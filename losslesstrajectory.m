function[v, h, h_apogee]= losslesstrajectory(m_s, m, m0, t_b, T, g0)

% --- INPUTS:
% m_s: vector of the structure mass of the different stage;
% m: vector of mass of each stage;
% m0: total mass of the launch vehicle;
% t_b: vector of burning times for each stage;
% T: vector of thrusts of each stage;
% g0: gravitational acceleration

% --- OUTPUTS:
% v: vector containing velocities after complete burn; [m/s]
% h: vector containing the heights reached after complete burn; [km]
% h_apogee: maximum reachable height [km]

% --- ASSUMPTIONS:
% No drag losses;
% Constant mass flow rate;
% Stage starts after the one before finished firing

% --- Solution
n = length(m_s); % number of stages
mp = []; % propellant mass vector
mpdot = []; % mass flow rate vector
c = []; % characteristic velocity vector
v0 = 0; % m/s
h0 = 0; % m
v = [v0];
h = [h0];
m_remaining = m0;
m_stage = 0;
for i = 1:n
    mp = [mp, m(i)-m_s(i)];
    m_before = m_remaining;
    m_after = m_remaining - mp(i);
    mpdot = [mpdot, mp(i)/t_b(i)];
    c = [c, T(i)/mpdot(i)];
    MR = m_after / m_before; 
    v(i+1) = v(i) - c(i) * log(MR) - g0 * t_b(i);
    h(i+1) = h(i) + v(i) * t_b(i) + c(i) * t_b(i) / mp(i) * m_before * (MR * log(MR) - MR + 1) - 0.5 * g0 * t_b(i)^2;
    m_remaining = m_after - m_s(i);
end
h_apogee = h(end) + 0.5 * v(end)^2 / g0;

h = 1e-3 * h;
h_apogee = h_apogee * 1e-3;