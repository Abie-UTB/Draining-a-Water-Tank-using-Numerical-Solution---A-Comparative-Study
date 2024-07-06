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

% Implicit Midpoint method
y_implicit_midpoint = zeros(1, length(t));
y_implicit_midpoint(1) = 8;   % Initial level of water
tic
for i = 1:n
    % Compute k1 using the given formula
    k1 = y_dot((t(i) + t(i + 1))/2, y_implicit_midpoint(i) + (h/2) * y_dot(t(i), y_implicit_midpoint(i)));

    % Update the solution using the given formula
    y_implicit_midpoint(i + 1) = y_implicit_midpoint(i) + h * k1;
end
toc
% 4th order range kutta method
y_rk4 = zeros(1, length(t));
y_rk4(1) = 8;
tic
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
toc
% 3-step Implicit Adams-Moulton method
y_adams_moulton_3step = zeros(1, length(t));
y_adams_moulton_3step(1:3) = y_true(1:3);  % Use true solution for initial steps
tic
for i = 3:n
    % Use fsolve to solve the implicit equation iteratively
    implicit_eq = @(x) x - y_adams_moulton_3step(i) - ...
        (h / 12) * (5 * y_dot(t(i + 1), x) + 8 * y_dot(t(i), y_adams_moulton_3step(i)) - ...
        y_dot(t(i - 1), y_adams_moulton_3step(i - 1)));
    y_adams_moulton_3step(i + 1) = fsolve(implicit_eq, y_adams_moulton_3step(i));
end
toc
% 3-step Explicit Adams-Moulton method
y_adams_moulton_explicit_3step = zeros(1, length(t));
y_adams_moulton_explicit_3step(1:3) = y_true(1:3);  % Use true solution for initial steps
tic
for i = 3:n
    % Calculate the next value using the explicit formula
    y_adams_moulton_explicit_3step(i + 1) = y_adams_moulton_explicit_3step(i) + ...
        h * (23 * y_dot(t(i), y_adams_moulton_explicit_3step(i)) - 16 * y_dot(t(i - 1), y_adams_moulton_explicit_3step(i - 1)) ...
        + 5 * y_dot(t(i - 2), y_adams_moulton_explicit_3step(i - 2))) / 12;
end
toc
% Predictor-Corrector Method (Adams-Bashforth-Moulton)
y_pred_corrected = zeros(1, length(t));
y_pred_corrected(1) = 8;   % Initial level of water
tic
for i = 1:n
    % Predictor Step (Adams-Bashforth)
    predictor = y_pred_corrected(i) + h * (3/2 * y_dot(t(i), y_pred_corrected(i)) - 1/2 * y_dot(t(i), y_pred_corrected(i)));
    
    % Corrector Step (Adams-Moulton)
    corrector = y_pred_corrected(i) + h * (1/2 * (y_dot(t(i + 1), predictor) + y_dot(t(i), y_pred_corrected(i))));
    
    y_pred_corrected(i + 1) = corrector;
end
toc
% Plotting the results
figure;
subplot(2,1,1);
plot(t, y_implicit_midpoint, 'b-', 'LineWidth', 2, 'DisplayName', 'Implicit Midpoint');
hold on;
% Plotting the results (continued)
plot(t, y_adams_moulton_3step, 'g-', 'LineWidth', 2, 'DisplayName', '3-step Implicit Adams-Moulton');
plot(t, y_adams_moulton_explicit_3step, 'r-', 'LineWidth', 2, 'DisplayName', '3-step Explicit Adams-Bashforth');
plot(t, y_rk4, '-ob', 'LineWidth', 1, 'DisplayName', '4th Order Runge-Kutta');
plot(t, y_pred_corrected, 'm-', 'LineWidth', 2, 'DisplayName', 'Predictor-Corrector (Adams-Bashforth-Moulton)');

% Add labels, legend, grid, etc.
xlabel('Time (s)');
ylabel('Water Level [m]');
grid on;
legend('show');
%title('Water Drainage Comparison');

% Calculate absolute errors
error_modified_euler = abs(y_true - y_implicit_midpoint);
error_adams_moulton_implicit = abs(y_true - y_adams_moulton_3step);
error_adams_moulton_explicit = abs(y_true- y_adams_moulton_explicit_3step);
error_rk4 = abs(y_true - y_rk4);
error_predictor_corrector = abs(y_true - y_pred_corrected);


subplot(2, 1, 2);
plot(t, error_modified_euler, '-*', 'Color', 'Red', 'LineWidth', 1, 'DisplayName', 'Implicit Midpoint Error');
hold on;
plot(t, error_adams_moulton_implicit, '-ob', 'LineWidth', 1, 'DisplayName', '3-step Implicit Adams-Moulton Error');
plot(t, error_adams_moulton_explicit, '-*g', 'LineWidth', 1, 'DisplayName', '3-step Explicit Adams-Bashforth Error');
plot(t, error_rk4, '-*m', 'LineWidth', 1, 'DisplayName', 'RK4 Error');
plot(t, error_predictor_corrector, '-oc', 'LineWidth', 1, 'DisplayName', 'Predictor-Corrector (Adams-Bashforth-Moulton) Error');
xlabel('Time [sec]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
ylabel('Local Error [m]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
legend('Location', 'Best');
%title('Absolute Error: Implicit Midpoint vs Adams-Moulton vs RK4 vs Predictor-Corrector', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 12);
grid on;