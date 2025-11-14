function [cneq, ceq] = nonlinear_const(z, x0, target, input, param, reference)

theta_vec = z(1:param.thetaGrid);
tf = z(end);

% === Setup controllo piecewise-constant ===
tkn  = linspace(0, 1, param.thetaGrid); % nodi temporali dei tratti

% Funzione: dato t, restituisce il theta del tratto corrente
warpFactor = tkn.^input.warpValuePower; 
tWarp = 0 + (tf-0) * warpFactor;
theta_of_t = @(t) interp1(tWarp, theta_vec, t, "linear", "extrap");

%theta_of_t = @(t) interp1(tkn,theta_vec,t,"linear","extrap");

% === Integrazione ODE con theta(t) variabile ===
userFun = @(t, x) dynamics_function(t, x, theta_of_t, input, param, reference);
tSpan = [0, tf];         
solution = ode45(userFun, tSpan, x0); % risolve il sistema di equazioni 

% contenuto in userFun integrando tra 0 e tf0. 

% Extract the solution on uniform grid:
t = linspace(solution.x(1), solution.x(end), param.nGrid); % crea vettore del tempo uniformemente equispaziato

x = deval(solution, t); % valuta la soluzione continua dell’ODE nei tempi scelti
yfinal = x(2,end); 
dxfinal = x(3,end); 
dyfinal = x(4,end); 

% A=-eye(param.thetaGrid);
% for i=1:param.thetaGrid-1
%     A(i,i+1)=1;
% end
% tStep=(solution.x(end)-solution.x(1))/(param.nGrid-1);
% dtheta=(A*z(1:param.thetaGrid)')./tStep;
% B=A;
% dtehta_diff=B*dtheta;
% for j=1:param.thetaGrid-1
%     dtheta_diff(j)=abs(dtehta_diff(j));
% end
% 
 cneq = [] ;%dtheta_diff-0.2*ones(param.thetaGrid-1,1);  %No inequality constraints
ceq = [ (yfinal - target.y)/target.y ...
        (dxfinal - target.vx)/target.vx ...
        (dyfinal - target.vy)/target.vx];  %Boundary Condition
        

end