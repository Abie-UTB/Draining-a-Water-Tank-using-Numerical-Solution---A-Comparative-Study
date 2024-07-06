clear all
close all
clc

% Constants
g = 9.81;  % Gravitational acceleration
D = 2;     % Diameter of the cylinder
d = 0.2;   % Diameter of the small hole at the bottom of the cylinder
h = 5;     % Step size
t = 0:h:127;  % Time interval to be evaluated
y_true = zeros(1, length(t));
y_true(1) = 8;   % True solution for the initial level of water
n = length(t) - 1;

y_dot = @(x, y) (-sqrt(2 * g) * (d / D)^2 * sqrt(y)); % Function to be solved

% True solution (for comparison)
for i = 1:n
    y_true(i + 1) = (sqrt(8) - sqrt(g/2) * (d / D)^2 * t(i+1)).^2;
end

% 4th Order Runge-Kutta method
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


% Calculate errors
error_rk4 = abs(y_true - y_rk4);
error_implicit_midpoint = abs(y_true - y_implicit_midpoint);

% Plotting
figure;
subplot(2, 1, 1);
plot(t, y_true, '-k', 'LineWidth', 2, 'DisplayName', 'True Solution');
hold on;
plot(t, y_rk4, '-^', 'Color', 'Green', 'LineWidth', 1, 'DisplayName', 'RK4');
plot(t, y_implicit_midpoint, '-d', 'Color', 'Magenta', 'LineWidth', 1, 'DisplayName', 'Implicit Midpoint');
xlabel('Time [sec]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
ylabel('Level of water [m]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
legend('Location', 'Best');
%title('Comparison of Solutions: True vs Numerical Methods', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 12);
grid on;

subplot(2, 1, 2);
plot(t, error_rk4, '-^', 'Color', 'Green', 'LineWidth', 1, 'DisplayName', 'RK4 Error');
hold on;
plot(t, error_implicit_midpoint, '-d', 'Color', 'Magenta', 'LineWidth', 1, 'DisplayName', 'Implicit Midpoint Error');
xlabel('Time [sec]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
ylabel('Local Error', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
legend('Location', 'Best');
%title('Local Error: Numerical Methods vs True Solution', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 12);
grid on;

% Calculate the global errors (e.g., using the L2 norm)
global_error_rk4 = norm(error_rk4, 2) / sqrt(length(y_true));
global_error_implicit_midpoint = norm(error_implicit_midpoint, 2) / sqrt(length(y_true));

% Display or use the global errors as needed
disp(['Global Error (RK4): ' num2str(global_error_rk4)]);
disp(['Global Error (Implicit Midpoint): ' num2str(global_error_implicit_midpoint)]);
% Find the maximum errors
max_error_rk4 = max(error_rk4);
max_error_implicit_midpoint = max(error_implicit_midpoint);

% Display the results
fprintf('Maximum Error (RK4): %f\n', max_error_rk4);
fprintf('Maximum Error (Implicit Midpoint): %f\n', max_error_implicit_midpoint);
