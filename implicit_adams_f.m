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

% Adams-Moulton Explicit Methods Initialization
y_adams_moulton_2step = zeros(1, length(t));
y_adams_moulton_2step(1:2) = y_true(1:2);  % Use true solution for initial steps
y_adams_moulton_3step = zeros(1, length(t));
y_adams_moulton_3step(1:3) = y_true(1:3);  % Use true solution for initial steps
y_adams_moulton_4step = zeros(1, length(t));
y_adams_moulton_4step(1:4) = y_true(1:4);  % Use true solution for initial steps
y_adams_moulton_5step = zeros(1, length(t));
y_adams_moulton_5step(1:5) = y_true(1:5);  % Use true solution for initial steps

% Adams-Moulton 2-step Method
for i = 2:n
    % Use fsolve to solve the implicit equation iteratively
    implicit_eq = @(x) x - y_adams_moulton_2step(i) - ...
        (h / 2) * (y_dot(t(i + 1), x) + y_dot(t(i), y_adams_moulton_2step(i)));
    y_adams_moulton_2step(i + 1) = fsolve(implicit_eq, y_adams_moulton_2step(i));
end

% Adams-Moulton 3-step Method
for i = 3:n
    % Use fsolve to solve the implicit equation iteratively
    implicit_eq = @(x) x - y_adams_moulton_3step(i) - ...
        (h / 12) * (5 * y_dot(t(i + 1), x) + 8 * y_dot(t(i), y_adams_moulton_3step(i)) - ...
        y_dot(t(i - 1), y_adams_moulton_3step(i - 1)));
    y_adams_moulton_3step(i + 1) = fsolve(implicit_eq, y_adams_moulton_3step(i));
end

% Adams-Moulton 4-step Method
for i = 4:n
    % Use fsolve to solve the implicit equation iteratively
    implicit_eq = @(x) x - y_adams_moulton_4step(i) - ...
        (h / 24) * (9 * y_dot(t(i + 1), x) + 19 * y_dot(t(i), y_adams_moulton_4step(i)) - ...
        5 * y_dot(t(i - 1), y_adams_moulton_4step(i - 1)) + y_dot(t(i - 2), y_adams_moulton_4step(i - 2)));
    y_adams_moulton_4step(i + 1) = fsolve(implicit_eq, y_adams_moulton_4step(i));
end

% Adams-Moulton 5-step Method
for i = 5:n
    % Use fsolve to solve the implicit equation iteratively
    implicit_eq = @(x) x - y_adams_moulton_5step(i) - ...
        (h / 720) * (251 * y_dot(t(i + 1), x) + 646 * y_dot(t(i), y_adams_moulton_5step(i)) - ...
        264 * y_dot(t(i - 1), y_adams_moulton_5step(i - 1)) + 106 * y_dot(t(i - 2), y_adams_moulton_5step(i - 2)) - ...
        19 * y_dot(t(i - 3), y_adams_moulton_5step(i - 3)));
    y_adams_moulton_5step(i + 1) = fsolve(implicit_eq, y_adams_moulton_5step(i));
end

% Plotting the results
figure;
subplot(2,1,1);
hold on;
 
plot(t, y_true, '-o', 'Color', 'Blue', 'LineWidth', 1, 'DisplayName', 'True Solution');
plot(t, y_adams_moulton_2step, '-s', 'Color', 'Red', 'LineWidth', 1, 'DisplayName', '2-step Implicit Adams-Moulton');
plot(t, y_adams_moulton_3step, '-^', 'Color', 'Green', 'LineWidth', 1, 'DisplayName', '3-step Implicit Adams-Moulton');
plot(t, y_adams_moulton_4step, '-d', 'Color', 'Magenta', 'LineWidth', 1, 'DisplayName', '4-step Implicit Adams-Moulton');
plot(t, y_adams_moulton_5step, 'm-', 'LineWidth', 2, 'DisplayName', '5-step Implicit Adams-Moulton');
xlabel('Time [sec.]');
xlabel('Time [sec]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
ylabel('Level of water [m]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
legend('Location', 'Best');
grid on;
hold off;

% Calculate errors
error_adams_moulton_2step = abs(y_true - y_adams_moulton_2step);
error_adams_moulton_3step = abs(y_true - y_adams_moulton_3step);
error_adams_moulton_4step = abs(y_true - y_adams_moulton_4step);
error_adams_moulton_5step = abs(y_true - y_adams_moulton_5step);

% Plotting the errors
subplot(2, 1, 2);
hold on;
plot(t, error_adams_moulton_2step, '-o', 'Color', 'Blue', 'LineWidth', 2, 'DisplayName', '2-step Implicit Adams-Moulton Error');
plot(t, error_adams_moulton_3step, '-s', 'Color', 'Red', 'LineWidth', 2, 'DisplayName', '3-step Implicit Adams-Moulton Error');
plot(t, error_adams_moulton_4step, '-^', 'Color', 'Green','LineWidth', 2, 'DisplayName', '4-step Implicit Adams-Moulton Error');
plot(t, error_adams_moulton_5step, '-d', 'Color', 'Magenta', 'LineWidth', 2, 'DisplayName', '5-step Implicit Adams-Moulton Error');
xlabel('Time [sec]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
ylabel('Local Error [m]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
title('Error Comparison of Adams-Moulton Methods');
legend;
grid on;
hold off;
