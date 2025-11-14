%% Flipped Classroom: Optimization of the Trajectory
clear; % Clear workspace
clc; % Clear command window
close all; % Close all figures

%%  --- DATA ----------------------------------------------------------------
inputVal.F = 3e6; % Thrust force in Newtons
inputVal.M0 = 1.4e5; % Initial mass in kilograms
inputVal.mDot = 801.6; % Mass flow rate in kg/s
inputVal.A = 7.069; % Cross-sectional area in m^2
inputVal.Cd = 0.5; % Drag coefficient

% Target Values
target.y = 200000; % Target altitude in meters
target.vx = 7.796e3; % Target horizontal velocity in m/s
target.vy = 0; % Target vertical velocity in m/s

% Reference Values
reference.rho = 1.225; % Air density in kg/m^3
reference.h = 8.44e3; % Reference height in meters

% Parameters
param.g = 9.80665; % Acceleration due to gravity in m/s^2
param.nGrid = 200; % Number of elements in the grid
param.thetaGrid = 8; % Number of theta values
tf_max = inputVal.M0 / inputVal.mDot; % Maximum burn time

inputVal.warpValuePower = 1;
%% Initial Condition

x0 = [0; 0; 0; 0.01];

theta0 = linspace(deg2rad(89.998), deg2rad(-20), param.thetaGrid);

tf0 = 0.95 * tf_max;  % Initial guess for final time

tkn = linspace(0, 1, param.thetaGrid); % Time vector for interpolation

warpFactor = tkn.^inputVal.warpValuePower; 
tWarp = 0 + (tf0-0) * warpFactor;
theta_of_t = @(t) interp1(tWarp, theta0, t, "linear", "extrap"); % Interpolated theta function

% === ODE Integration with variable theta(t) ===
userFun = @(t, x) dynFunctRocket(t, x, theta_of_t, inputVal, param, reference); % Dynamics function
tSpan = [0, tf0];         
solution = ode45(userFun, tSpan, x0); % Solve the system of equations 

% Extract the solution on a uniform grid:
sol.t = linspace(solution.x(1), solution.x(end), param.nGrid); % Create a uniformly spaced time vector

xx = deval(solution, sol.t); % Evaluate the continuous solution of the ODE at chosen times
sol.x = xx(1,:); % Extract x positions
sol.y = xx(2,:); % Extract y positions
sol.dx = xx(3,:); % Extract horizontal velocities
sol.dy = xx(4,:); % Extract vertical velocities
sol.v = sqrt(sol.dx.^2 + sol.dy.^2); % Calculate the resultant velocity
guess.tf = sol.t(end); % Store the final time guess

%% --------------------Initial Guess Figure--------------------------------

figure
plot(sol.x / 1000, sol.y / 1000) % Plot trajectory in km
hold on
yline(target.y / 1000) % Plot target altitude
hold off
ylim([0, 250]); % Set y-axis limits
xlabel("Horizontal Distance [Km]") % Label x-axis
ylabel("Vertical Distance [Km]") % Label y-axis
grid on % Enable grid

figure
plot(sol.t, sol.dx) % Plot horizontal velocity over time
hold on
yline(target.vx) % Plot target horizontal velocity
xlabel('Time [s]') % Label x-axis
ylabel('Velocity [m/s]') % Label y-axis
title('Velocity Profile Over Time') % Title of the plot
grid on % Enable grid

figure
plot(sol.t, sol.dy) % Plot vertical velocity over time
hold on
yline(target.vy) % Plot target vertical velocity
xlabel('Time [s]') % Label x-axis
ylabel('Velocity [m/s]') % Label y-axis
title('Velocity Profile Over Time') % Title of the plot
grid on % Enable grid

% Check if the initial guess is satisfactory
satisfactory = input('Is the initial guess satisfactory? (yes/no):', 's'); % User input for guess satisfaction

if strcmpi(satisfactory, 'yes') % If the guess is satisfactory
    disp('Continuing with the optimization...'); % Display message
    close all % Close all figures
else
    disp('Stopping the execution. Please adjust the initial guess.'); % Display message
    return; % Stop the execution if the guess is not satisfactory
end

%% ---------------------Optimisation Section------------------------------

% The for cycle is intended to increase the number of elements in the theta
% grid by updating the initial guess at every cycle in order to facilitate 
% the convergence of fmincon  
theta0new=theta0;
tWarpnew=tWarp;



for n=1:4
    theta0old=theta0new;
    tWarpold=tWarpnew;
    z0 = [theta0old guess.tf]; % Initial guess for optimization
    
    problem.x0 = z0; % Initial guess for optimization variables
    
    % Parameters for fmincon
    problem.Aineq = []; 
    problem.Bineq = [];  % No linear inequality constraints
    problem.Aeq = []; 
    problem.Beq = [];  % No equality constraints
    lb = [-pi/2 * ones(1, param.thetaGrid) 1]; % Lower bounds for theta and tf
    ub = [+pi/2 * ones(1, param.thetaGrid) tf_max]; % Upper bounds for theta and tf
    problem.lb = lb; % Set lower bounds
    problem.ub = ub; % Set upper bounds
    problem.solver = 'fmincon'; % Optimization solver
    problem.options = optimoptions('fmincon', ...
       'Algorithm', 'interior-point', 'EnableFeasibilityMode', true, ...
       'Display', 'iter', 'MaxFunctionEvaluations', 2e6, ...
       'StepTolerance', 1e-19, 'OptimalityTolerance', 1e-8, 'ConstraintTolerance', 1e-10); % Solver options
    
    problem.objective = @(z) z(end); % Objective function to minimize (final time)
    problem.nonlcon = @(z) nonlinear_const(z, x0, target, inputVal, param, reference); % Nonlinear constraints
    
    [xSoln, fVal, exitFlag, output] = fmincon(problem); % Run optimization
    
    tfOpt = fVal; % Optimal final time
    thetaOpt_vec = xSoln(1:param.thetaGrid); % Optimal theta values
    tknOpt = linspace(0, tfOpt, param.thetaGrid); % Time vector for optimal solution
    
    param.thetaGrid=ceil(param.thetaGrid*1.6);
    
    tkn = linspace(0, 1, param.thetaGrid); % Time vector for interpolation
    
    warpFactor = tkn.^inputVal.warpValuePower; 
    tWarpnew = 0 + (tf0-0) * warpFactor;
    
    theta0new=interp1(tknOpt,thetaOpt_vec,tWarpnew,'spline');
    theta_of_t = @(t) interp1(tWarpnew, theta0new, t, "spline"); % Interpolated theta function
    guess.tf=tfOpt;
    if n==1
        theta0fCd=xSoln;
    end
end
%% -------------------------Plot Results----------------------------------

figure(4)
% Plot optimal theta profile vs time
plot(tknOpt, thetaOpt_vec * 180 / pi, 'b-o', 'LineWidth', 1.5); % Optimal theta in degrees
hold on
% Plot initial theta profile vs time
plot(tWarpold, theta0old * 180 / pi, 'r--'); % Initial theta in degrees
xlabel('Time[s]'); % Label x-axis
ylabel('Theta [°]'); % Label y-axis
legend(["Optimal Solution", "Initial Guess"]); % Legend
grid on; % Enable grid
title('Optimal Theta vs Initial Guess'); % Title of the plot

%% ------------ Final Run to Retrieve Results------------------------------

thetatOpt = @(t) interp1(tknOpt, thetaOpt_vec, t, "linear", "extrap"); 
userFun = @(t, x) dynamics_function(t, x, thetatOpt, inputVal, param, reference); 
tSpanOpt = [0, fVal]; 
x0opt = x0; 
solutionOpt = ode45(userFun, tSpanOpt, x0opt); 
solOpt.t = linspace(solutionOpt.x(1), solutionOpt.x(end), param.nGrid); 
xxOpt = deval(solutionOpt, solOpt.t);

solOpt.x = xxOpt(1,:); 
solOpt.y = xxOpt(2,:);
solOpt.dx = xxOpt(3,:); 
solOpt.dy = xxOpt(4,:); 

%% ----------------------- Optimisation Plots ------------------------------

figure(5)
plot(solOpt.x / 1000, solOpt.y / 1000, 'b-', 'LineWidth', 1.5, ...
     'DisplayName', 'Optimal Trajectory'); % Plot optimal trajectory in km
hold on
yline(target.y / 1000, 'r--', 'LineWidth', 1.5, ...
      'DisplayName', 'Target Altitude'); % Plot target altitude
hold off
ylim([0, 250]); 
xlabel("Horizontal Distance [Km]");
ylabel("Altitude [Km]"); 
title("Trajectory"); 
legend('Location', 'northeastoutside'); 
grid on; 

figure(6)
plot(solOpt.t, solOpt.dx, 'b-', 'LineWidth', 1.5);
hold on
yline(target.vx, 'r--', 'LineWidth', 1.5); 
hold off
xlabel('Time [s]'); 
ylabel('Horizontal Velocity [m/s]'); 
title('Horizontal Velocity $v_x$','Interpreter','latex'); 
legend(["Optimal Velocity $v_x$","Target $v_x$"], 'Location', 'northeastoutside', 'Interpreter', 'latex'); 
grid on; 

figure(7)
plot(solOpt.t, solOpt.dy, 'b-', 'LineWidth', 1.5);
hold on
yline(target.vy, 'r--', 'LineWidth', 1.5); 
hold off
xlabel('Time [s]'); 
ylabel('Vertical Velocity [m/s]');
title('Vertical Velocity $v_y$','Interpreter','latex'); 
legend(["Optimal Velocity $v_y$","Target $v_y$"], 'Location', 'northeastoutside', 'Interpreter', 'latex');
grid on; 

%% Thrust-over-Weight Ratio Evolution

massEvolution = inputVal.M0 - inputVal.m_dot .* solOpt.t;
weight = massEvolution * param.g;
thrustOverWeight = inputVal.F ./ weight;

figure()
plot(solOpt.t, thrustOverWeight, 'b', 'LineWidth', 1.5)
hold on 
yline(1, '--k', 'LineWidth', 1.5)
grid on
title('T/W ratio evolution')
xlabel('Time [s]')
ylabel('T/W')
legend('T/W', 'T/W=1', 'Location','northeastoutside')

%% evolution of theta with cd
param.thetaGrid = 8; % Number of theta values
Ncurves = 9;
Cd  = linspace(0.2,1,Ncurves);

figure; hold on


% Crea una colormap con Ncurves colori distinti
cmap = turbo(Ncurves);         % o parula, jet, viridis, etc.

for i=1:length(Cd)
    inputVal.Cd=Cd(i);
    z0 = [theta0fCd]; % Initial guess for optimization
    
    problem.x0 = z0; % Initial guess for optimization variables
    
    % Parameters for fmincon
    problem.Aineq = []; 
    problem.Bineq = [];  % No linear inequality constraints
    problem.Aeq = []; 
    problem.Beq = [];  % No equality constraints
    lb = [-pi/2 * ones(1, param.thetaGrid) 1]; % Lower bounds for theta and tf
    ub = [+pi/2 * ones(1, param.thetaGrid) tf_max]; % Upper bounds for theta and tf
    problem.lb = lb; % Set lower bounds
    problem.ub = ub; % Set upper bounds
    problem.solver = 'fmincon'; % Optimization solver
    problem.options = optimoptions('fmincon', ...
       'Algorithm', 'interior-point', 'EnableFeasibilityMode', true, ...
       'Display', 'iter', 'MaxFunctionEvaluations', 2e6, ...
       'StepTolerance', 1e-19, 'OptimalityTolerance', 1e-8, 'ConstraintTolerance', 1e-8); % Solver options
    
    problem.objective = @(z) z(end); % Objective function to minimize (final time)
    problem.nonlcon = @(z) nonlinear_const(z, x0, target, inputVal, param, reference); % Nonlinear constraints
    
    [xSoln, fVal, exitFlag, output] = fmincon(problem); % Run optimization
    tfOpt = fVal; % Optimal final time
    thetaOpt_vec = xSoln(1:param.thetaGrid); % Optimal theta values
    tknOpt = linspace(0, tfOpt, param.thetaGrid); % Time vector for optimal solution
    plot(tknOpt,xSoln(1:param.thetaGrid), 'Color', cmap(i,:))
    hold on

end


colormap(turbo)
cb = colorbar;
clim([min(Cd) max(Cd)])           % imposta i limiti
cb.Ticks = Cd;                       % etichette esatte
cb.TickLabels = string(Cd);
cb.Label.String = 'Cd';
grid on;
xlabel('x'); ylabel('y'); zlabel('z')


