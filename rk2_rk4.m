% Constants
g = 9.81;  % Gravitational acceleration
D = 2;     % Diameter of the cylinder
d = 0.2;   % Diameter of the small hole at the bottom of the cylinder
h = 5;     % Step size
tend = 127; % End time

% Initialization
t = 0:h:tend;  % Time interval to be evaluated
n = length(t) - 1;

% Function to be solved
y_dot = @(t, y) -sqrt(2 * g) * (d / D)^2 * sqrt(y);

% True solution (for comparison)
y_true = zeros(1, length(t));
y_true(1) = 8;   % Initial level of water

for i = 1:n
    y_true(i + 1) = (sqrt(8) - sqrt(g/2) * (d / D)^2 * t(i+1)).^2;
end

% RK2 method
y_rk2 = zeros(1, length(t));
y_rk2(1) = 8;   % Initial level of water

for i = 1:n
    % Compute slope k1
    k1 = y_dot(t(i), y_rk2(i));
    
    % Compute slope k2 using k1
    k2 = y_dot(t(i) + h, y_rk2(i) + h * k1);
    
    % Update the solution using weighted average of slopes
    y_rk2(i + 1) = y_rk2(i) + h * (k1 + k2) / 2;
end

% RK4 method
y_rk4 = zeros(1, length(t));
y_rk4(1) = 8;   % Initial level of water

for i = 1:n
    % Compute slope k1
    k1 = y_dot(t(i), y_rk4(i));
    
    % Compute slope k2 using k1
    k2 = y_dot(t(i) + 0.5 * h, y_rk4(i) + 0.5 * h * k1);
    
    % Compute slope k3 using k2
    k3 = y_dot(t(i) + 0.5 * h, y_rk4(i) + 0.5 * h * k2);
    
    % Compute slope k4 using k3
    k4 = y_dot(t(i) + h, y_rk4(i) + h * k3);
    
    % Update the solution using weighted average of slopes
    y_rk4(i + 1) = y_rk4(i) + h * (k1 + 2*k2 + 2*k3 + k4) / 6;
end

% Implicit Midpoint method using given formula
y_implicit_midpoint = zeros(1, length(t));
y_implicit_midpoint(1) = 8;   % Initial level of water

for i = 1:n
    % Compute k1 using the given formula
    k1 = y_dot((t(i) + t(i + 1))/2, y_implicit_midpoint(i) + (h/2) * y_dot(t(i), y_implicit_midpoint(i)));

    % Update the solution using the given formula
    y_implicit_midpoint(i + 1) = y_implicit_midpoint(i) + h * k1;
end

% Calculate absolute errors
error_rk2 = abs(y_true - y_rk2);
error_rk4 = abs(y_true - y_rk4);
error_implicit_midpoint = abs(y_true - y_implicit_midpoint);

% Plotting the results
figure;

% Plotting the water level
subplot(2,1,1);
plot(t, y_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True Solution');
hold on;
plot(t, y_rk2, 'bo', 'LineWidth', 1.5, 'DisplayName', 'RK2');
plot(t, y_rk4, 'r*', 'LineWidth', 1.5, 'DisplayName', 'RK4');
plot(t, y_implicit_midpoint, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Implicit Midpoint');
xlabel('Time (s)');
ylabel('Water Level');
grid on;
legend('show');

% Plotting the error
subplot(2,1,2);
plot(t, error_rk2, 'b-', 'LineWidth', 1.5, 'DisplayName', 'RK2 Error');
hold on;
plot(t, error_rk4, 'r-', 'LineWidth', 1.5, 'DisplayName', 'RK4 Error');
plot(t, error_implicit_midpoint, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Implicit Midpoint Error');
xlabel('Time (s)');
ylabel('Absolute Error [m]');
grid on;
legend('show');
