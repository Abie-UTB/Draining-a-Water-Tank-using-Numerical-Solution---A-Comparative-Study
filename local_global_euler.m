% Given code for the water drainage system

% Constants
g = 9.81;  % Gravitational acceleration
D = 2;     % Diameter of the cylinder
d = 0.2;   % Diameter of the small hole at the bottom of the cylinder
h = 5;     % Step size
t = 0:h:127;  % Time interval to be evaluated

% Initialization
y_euler = zeros(1, length(t));
y_exact = (sqrt(8) - sqrt(g/2) * (d / D)^2 * t).^2;

y_euler(1) = 8;   % Initial level of water

n = length(t) - 1;

% Function to be solved
y_dot = @(t, y) -sqrt(2 * g) * (d / D)^2 * sqrt(y);

% Euler method
for i = 1:n
    y_euler(i + 1) = y_euler(i) + h * y_dot(t(i), y_euler(i));
end

% Calculate local and global errors
local_errors = abs(y_exact - y_euler(1:length(t)));

global_error = max(local_errors);

% Plotting
figure;

% Plot the numerical solution
subplot(2, 1, 1);
plot(t, y_euler, 'b-', 'LineWidth', 2);
hold on;
plot(t, y_exact, 'r--', 'LineWidth', 2);
legend('Numerical Solution (Euler)', 'Exact Solution');
xlabel('Time');
ylabel('Water Level');
title('Numerical vs. Exact Solutions');

% Plot the local error
subplot(2, 1, 2);
plot(t, local_errors, 'g-', 'LineWidth', 2);
xlabel('Time');
ylabel('Local Error');
title('Local Error over Time');

% Display the global error
fprintf('Global Error: %f\n', global_error);
