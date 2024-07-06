% Define the ODE function


% Set parameters
g = 9.81; % gravitational acceleration
D = 2;    % some constant
d = 0.2;  % some constant

% Set initial condition
h0 = 8;

% Time span for simulation
tspan = [0 128];
ode = @(t, h) -sqrt(2 * g) * (d / D)^2 * sqrt(h);
% Call ode45 with the odeprint option
options = odeset('OutputFcn',@odeprint);
[t, h] = ode45(ode, tspan, h0, options);

% Plot the solution
figure;
plot(t, h, '-o');
xlabel('Time');
ylabel('h(t)');
title('Solution of the ODE');

%