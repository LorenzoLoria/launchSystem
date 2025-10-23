clear;
clc;

%% PARETO FRONT

% Sometimes we have to find the optimal payload and the optimal altitude
% --> the problem is that sometimes we could have a condition in which our
% output has an opposite behaviour (for ex. for an optimal altitude we have
% a non optimal payload) --> in this case we obtain the so called PARETO
% FRONT so we collect not only one optimal solution but a set of optimal
% solutions that combined will give us a certain plot --> in this case, to 
% find a proper solution, the designer should add other requirements of the
% mission or to understand which of the solutions could be reliable or not.
%
% In this case we will have no more a scalar objectove function, but we
% will have to perform a multi-objective optimisation ('gamultiobj').
% 
% Here we have two objectives to oprimise: payload and altitute
%
% In this case we have an objective function in a vectorial form --> we can 
% transform it in a scalar form introducing a weight factor and so we
% return to the previous case and we can perform a single objactive
% optimisation --> we will collect different solutions depending on the 
% weight that we set at the beginning --> then we can perform a fitting to
% obtain the Pareto front.
% 
% Alternatively we can use a proper function (called "multi-objective 
% generic algorithm") that is a function that is optimised to perform the
% Pareto front without using this transformation of the objective function.




% Definition of the problem --> OBJ: optimise altitude and payload
rocket.Isp = 280;
rocket.m0 = 12300;
rocket.T = 450000;
% rocket.mpay = 800;

% Initial condition
rocket.v0 = 0;
rocket.h0 = 0;

% Definition of options

% Note that ('PlotFcn','gaplotpareto') is an external function that Matlab 
% uses to generate the Pareto front at each step --> so you can see in live 
% how the Pareto front change during the optimisation
optionsGA_mult = optimoptions('gamultiobj', 'Display', 'iter', 'ConstraintTolerance', 1e-6, 'FunctionTolerance', 1e-6, 'PlotFcn','gaplotpareto');


[x_pareto, out_w, exitflag, output] = gamultiobj(@(x) fun_pareto(x, rocket), 1, [], [], [], [], 0.7, 0.9, @(x) non_con(x, rocket), optionsGA_mult)
% where:
%  - x_pareto: optimal propellant mass fraction
%  - out_w: set of solutions that the Pareto front gives us (evaluation of
%           the objective function for the set of variables defined by x_pareto)





% PLOTS 

% 1st figure --> Payload vs Apogee 
figure('Color', 'w')
grid on
hold on
plot(-out_w(:, 2), -out_w(:, 1)/1000, 'k.')
xlabel('Payload [Kg]')
ylabel('Apogee [Km]')

% 2nd figure 
% left --> comparison between the propellant mass fraction and the
%          payload, to see how the payload changes changing the  
%          propelland mass fraction
% right --> comparison between the propellant mass fraction and the
%           altitude, to see how the altitude changes changing the  
%           propelland mass fraction
figure('Color', 'w')
grid on
hold on
xlabel('Propellant mass fraction')

yyaxis left
plot(x_pareto, -out_w(:, 2), 'k-')
ylabel('Payload [Kg]')

yyaxis right
plot(x_pareto, -out_w(:, 1)/1000, '.', 'MarkerEdgeColor', "#D95319")
ylabel('Altitude [Km]')





% ---------------------------- FUNCTIONS ----------------------------------

% Definition of objective function
function [obj] = fun_pareto(var, rocket)

eps_p = var;

rocket.mp = eps_p * rocket.m0; % propellant mass
rocket.tb = rocket.mp * rocket.Isp * 9.81 / rocket.T; % burning time
rocket.c = rocket.Isp * 9.81; 

MR = (rocket.m0 - rocket.mp) / rocket.m0; % mass ratio

v_ct = rocket.v0 - rocket.c * log(MR) - 9.81 * rocket.tb; % speed at cut-off time
h_ct = rocket.h0 + rocket.v0 * rocket.tb + rocket.c * rocket.m0 / rocket.mp * rocket.tb * (MR * log(MR) - MR + 1) - 9.81 * rocket.tb^2 / 2; % altitude at cut-off time

h_apogee = h_ct + 1/2 * v_ct^2 / 9.81;
m_pay = rocket.m0 - rocket.mp - MER_solid(rocket.mp);

% The difference is here --> we want the maximisation of h_apogee and the
% maximisation of the payload
obj = [- h_apogee, -m_pay]; 

end




% Definition of the structural mass
function [ms] = MER_solid(mp)

pc = 7; % assumption

params_Me = [800.340063456389   0.0970643753730160  ...
             -162.431518885189  -0.0180992880437423  -8.48759341267327  -6441.47167561289];

Minert = @(params_Me, mp, pc) params_Me(1) + params_Me(2).* mp + params_Me(3) .* pc ...
    + params_Me(4) .* (pc + params_Me(5)).*(mp + params_Me(6));

if mp > 5000
    ms = Minert(params_Me, mp, pc);
else
    ms = 140.94.*log(mp) - 823.29;
end

end




% NON-linear constraint function
function [c, ceq] = non_con(var, rocket)

eps_p = var;

rocket.mp = eps_p * rocket.m0; % propellant mass
rocket.tb = rocket.mp * rocket.Isp * 9.81 / rocket.T; % burning time
rocket.c = rocket.Isp * 9.81; 

MR = (rocket.m0 - rocket.mp) / rocket.m0; % mass ratio

v_ct = rocket.v0 - rocket.c * log(MR) - 9.81 * rocket.tb; % speed at cut-off time
h_ct = rocket.h0 + rocket.v0 * rocket.tb + rocket.c * rocket.m0 / rocket.mp * rocket.tb * (MR * log(MR) - MR + 1) - 9.81 * rocket.tb^2 / 2; % altitude at cut-off time

h_apogee = h_ct + 1/2 * v_ct^2 / 9.81;
m_pay = rocket.m0 - rocket.mp - MER_solid(rocket.mp);


% The only difference will be the inequality non linear constraint --> h_apogee and m_pay should be >0 (so we put a minus) 
c = [- h_apogee, - m_pay];
ceq = [];

end