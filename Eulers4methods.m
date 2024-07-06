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

% Explicit Euler method
y_explicit_euler = zeros(1, length(t));
y_explicit_euler(1) = 8;   % Initial level of water

for i = 1:n
    y_explicit_euler(i + 1) = y_explicit_euler(i) + h * y_dot(t(i), y_explicit_euler(i));
end

% Implicit Euler method
y_impl_eu = zeros(1, length(t));
y_impl_eu(1) = 8;   % Initial level of water

for i = 1:n
    % Calculate slope at tn, yn
    k1 = y_dot(t(i), y_impl_eu(i));
    % Update yn+1 using the implicit method
    y_impl_eu(i+1) = fzero(@(y) y - y_impl_eu(i) - h * y_dot(t(i+1), y), y_impl_eu(i));
end



% Implicit Midpoint method using given formula
y_implicit_mp = zeros(1, length(t));
y_implicit_mp(1) = 8;   % Initial level of water

for i = 1:n
    % Compute k1 using the given formula
    k1 = y_dot((t(i) + t(i + 1))/2, y_implicit_mp(i) + (h/2) * y_dot(t(i), y_implicit_mp(i)));

    % Update the solution using the given formula
    y_implicit_mp(i + 1) = y_implicit_mp(i) + h * k1;
end

% Modified Euler method (Heun's method)
y_modified_euler = zeros(1, length(t));
y_modified_euler(1) = 8;   % Initial level of water

for i = 1:n
    k1 = y_dot(t(i), y_modified_euler(i));
    k2 = y_dot(t(i) + h, y_modified_euler(i) + h * k1);
    y_modified_euler(i + 1) = y_modified_euler(i) + 0.5 *h* (k1 + k2);
end

% Calculate errors
error_explicit_euler = abs(y_true - y_explicit_euler);
error_implicit_euler = abs(y_true - y_impl_eu);
error_implicit_midpoint = abs(y_true - y_implicit_mp);
error_modified_euler = abs(y_true - y_modified_euler);

% Plotting
figure('Name', 'Comparison of Solutions and Errors', 'NumberTitle', 'off');

% Subplot 1: Comparison of Solutions
subplot(2, 1, 1);
plot(t, y_true, '-k', 'LineWidth', 2, 'DisplayName', 'True Solution');
hold on;
plot(t, y_explicit_euler, '-o', 'Color', 'Blue', 'LineWidth', 1, 'DisplayName', 'Explicit Euler');
plot(t, y_impl_eu, '-s', 'Color', 'Red', 'LineWidth', 1, 'DisplayName', 'Implicit Euler');
plot(t, y_implicit_mp, '-^', 'Color', 'Green', 'LineWidth', 1, 'DisplayName', 'Implicit Midpoint');
plot(t, y_modified_euler, '-d', 'Color', 'Magenta', 'LineWidth', 1, 'DisplayName', 'Modified Euler');
xlabel('Time [sec]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
ylabel('Level of water [m]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
legend('Location', 'Best');
%title('Comparison of Solutions: True vs Numerical Methods', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 12);
grid on;

% Subplot 2: Comparison of Errors
subplot(2, 1, 2);
plot(t, error_explicit_euler, '-o', 'Color', 'Blue', 'LineWidth', 1, 'DisplayName', 'Explicit Euler Error');
hold on;
plot(t, error_implicit_euler, '-s', 'Color', 'Red', 'LineWidth', 1, 'DisplayName', 'Implicit Euler Error');
plot(t, error_implicit_midpoint, '-^', 'Color', 'Green', 'LineWidth', 1, 'DisplayName', 'Implicit Midpoint Error');
plot(t, error_modified_euler, '-d', 'Color', 'Magenta', 'LineWidth', 1, 'DisplayName', 'Modified Euler Error');
xlabel('Time [sec]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
ylabel('Local Error [m]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
legend('Location', 'Best');
%title('Local Error: Numerical Methods vs True Solution', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 12);
grid on;

% Adjust figure properties for better visualization
set(gcf, 'Position', [100, 100, 800, 800]); % Set figure size
set(gca, 'FontName', 'Arial', 'FontSize', 10); % Set font properties

% Calculate the global errors (e.g., using the L2 norm)
global_error_explicit_euler = norm(error_explicit_euler, 2) / sqrt(length(y_true));
global_error_implicit_euler = norm(error_implicit_euler, 2) / sqrt(length(y_true));
global_error_implicit_midpoint = norm(error_implicit_midpoint, 2) / sqrt(length(y_true));
global_error_modified_euler = norm(error_modified_euler, 2) / sqrt(length(y_true));

% Display or use the global errors as needed
disp(['Global Error (Explicit Euler): ' num2str(global_error_explicit_euler)]);
disp(['Global Error (Implicit Euler): ' num2str(global_error_implicit_euler)]);
disp(['Global Error (Implicit Midpoint): ' num2str(global_error_implicit_midpoint)]);
disp(['Global Error (Modified Euler): ' num2str(global_error_modified_euler)]);
% Find the maximum errors
max_error_explicit_euler = max(error_explicit_euler);
max_error_implicit_euler = max(error_implicit_euler);
max_error_implicit_midpoint = max(error_implicit_midpoint);
max_error_modified_euler = max(error_modified_euler);

% Display the results
fprintf('Maximum Error (Explicit Euler): %f\n', max_error_explicit_euler);
fprintf('Maximum Error (Implicit Euler): %f\n', max_error_implicit_euler);
fprintf('Maximum Error (Implicit Midpoint): %f\n', max_error_implicit_midpoint);
fprintf('Maximum Error (Modified Euler): %f\n', max_error_modified_euler);