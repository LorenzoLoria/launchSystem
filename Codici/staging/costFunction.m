function [objective] = costFunction(z,x0,  target, input, param, reference)

theta_vec = z(1:param.nGrid);
tf = z(end);

% === Setup controllo piecewise-constant ===
tkn  = linspace(0, tf, param.nGrid); % nodi temporali dei tratti

% Funzione: dato t, restituisce il theta del tratto corrente

theta_of_t = @(t) interp1(tkn,theta_vec,t,"linear");

% === Integrazione ODE con theta(t) variabile ===
userFun = @(t, x) dynamics_function(t, x, theta_of_t, input, param, reference);
tSpan = [0, tf];         
solution = ode45(userFun, tSpan, x0); % risolve il sistema di equazioni 

% contenuto in userFun integrando tra 0 e tf0. 

% Extract the solution on uniform grid:
t = linspace(solution.x(1), solution.x(end), param.nGrid); % crea vettore del tempo uniformemente equispaziato

objective = t(end);

end