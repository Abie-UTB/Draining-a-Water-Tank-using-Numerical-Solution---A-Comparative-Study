clear all
close all
clc

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

% 2-step Explicit Adams-Moulton method
y_adams_moulton_explicit_2step = zeros(1, length(t));
y_adams_moulton_explicit_2step(1:2) = y_true(1:2);  % Use true solution for initial steps
for i = 2:n
    % Calculate the next value using the explicit formula
    y_adams_moulton_explicit_2step(i + 1) = y_adams_moulton_explicit_2step(i) + ...
        h * (3 * y_dot(t(i), y_adams_moulton_explicit_2step(i)) - y_dot(t(i - 1), y_adams_moulton_explicit_2step(i - 1))) / 2;
end

% 3-step Explicit Adams-Moulton method
y_adams_moulton_explicit_3step = zeros(1, length(t));
y_adams_moulton_explicit_3step(1:3) = y_true(1:3);  % Use true solution for initial steps
for i = 3:n
    % Calculate the next value using the explicit formula
    y_adams_moulton_explicit_3step(i + 1) = y_adams_moulton_explicit_3step(i) + ...
        h * (23 * y_dot(t(i), y_adams_moulton_explicit_3step(i)) - 16 * y_dot(t(i - 1), y_adams_moulton_explicit_3step(i - 1)) ...
        + 5 * y_dot(t(i - 2), y_adams_moulton_explicit_3step(i - 2))) / 12;
end

% 4-step Explicit Adams-Moulton method
y_adams_moulton_explicit_4step = zeros(1, length(t));
y_adams_moulton_explicit_4step(1:4) = y_true(1:4);  % Use true solution for initial steps
for i = 4:n
    % Calculate the next value using the explicit formula
    y_adams_moulton_explicit_4step(i + 1) = y_adams_moulton_explicit_4step(i) + ...
        h * (55 * y_dot(t(i), y_adams_moulton_explicit_4step(i)) - 59 * y_dot(t(i - 1), y_adams_moulton_explicit_4step(i - 1)) ...
        + 37 * y_dot(t(i - 2), y_adams_moulton_explicit_4step(i - 2)) - 9 * y_dot(t(i - 3), y_adams_moulton_explicit_4step(i - 3))) / 24;
end

% 5-step Explicit Adams-Moulton method
y_adams_moulton_explicit_5step = zeros(1, length(t));
y_adams_moulton_explicit_5step(1:5) = y_true(1:5);  % Use true solution for initial steps
for i = 5:n
    % Calculate the next value using the explicit formula
    y_adams_moulton_explicit_5step(i + 1) = y_adams_moulton_explicit_5step(i) + ...
        h * (1901 * y_dot(t(i), y_adams_moulton_explicit_5step(i)) - 2774 * y_dot(t(i - 1), y_adams_moulton_explicit_5step(i - 1)) ...
        + 2616 * y_dot(t(i - 2), y_adams_moulton_explicit_5step(i - 2)) - 1274 * y_dot(t(i - 3), y_adams_moulton_explicit_5step(i - 3)) ...
        + 251 * y_dot(t(i - 4), y_adams_moulton_explicit_5step(i - 4))) / 720;
end
% Change the display format to long
format long;
% Calculate absolute errors
error_adams_moulton_explicit_2step = abs(y_true - y_adams_moulton_explicit_2step);
error_adams_moulton_explicit_3step = abs(y_true - y_adams_moulton_explicit_3step);
error_adams_moulton_explicit_4step = abs(y_true - y_adams_moulton_explicit_4step);
error_adams_moulton_explicit_5step = abs(y_true - y_adams_moulton_explicit_5step);

% Plotting
figure;
subplot(2, 1, 1);
plot(t, y_true, '-k', 'LineWidth', 2, 'DisplayName', 'True Solution');
hold on;
plot(t, y_adams_moulton_explicit_2step, '-o', 'Color', 'Blue', 'LineWidth', 1, 'DisplayName', '2-step Adams-Bashforth(Explicit)');
plot(t, y_adams_moulton_explicit_3step, '-s', 'Color', 'Red', 'LineWidth', 1, 'DisplayName', '3-step Adams-Bashforth(Explicit)');
plot(t, y_adams_moulton_explicit_4step, '-^', 'Color', 'Green', 'LineWidth', 1, 'DisplayName', '4-step Adams-Bashforth(Explicit)');
plot(t, y_adams_moulton_explicit_5step, '-d', 'Color', 'Magenta', 'LineWidth', 1, 'DisplayName', '5-step Adams-Bashforth(Explicit)');
xlabel('Time [sec]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
ylabel('Level of water [m]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
legend('Location', 'Best');
%title('Comparison of Solutions: True vs Adams-Moulton (2 to 5 steps, Explicit)', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 12);
grid on;

subplot(2, 1, 2);
plot(t, error_adams_moulton_explicit_2step, '-o', 'Color', 'Blue', 'LineWidth', 1, 'DisplayName', '2-step Adams-Bashforth(Explicit)');
hold on;
plot(t, error_adams_moulton_explicit_3step, '-s', 'Color', 'Red', 'LineWidth', 1, 'DisplayName', '3-step Adams-Bashforth(Explicit)');
plot(t, error_adams_moulton_explicit_4step, '-^', 'Color', 'Green', 'LineWidth', 1, 'DisplayName', '4-step Adams-Bashforth(Explicit)');
plot(t, error_adams_moulton_explicit_5step, '-d', 'Color', 'Magenta', 'LineWidth', 1, 'DisplayName', '5-step Adams-Bashforth(Explicit)');
xlabel('Time [sec]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
ylabel('Local Error [m]', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 10);
legend('Location', 'Best');
%title('Comparison of Absolute Error: Adams-Moulton (2 to 5 steps, Explicit)', 'fontweight', 'bold', 'color', 'Black', 'fontsize', 12);
grid on;

% Calculate global errors (L2 norm)
global_error_explicit_2step = norm(error_adams_moulton_explicit_2step, 2) / sqrt(length(y_true));
global_error_explicit_3step = norm(error_adams_moulton_explicit_3step, 2) / sqrt(length(y_true));
global_error_explicit_4step = norm(error_adams_moulton_explicit_4step, 2) / sqrt(length(y_true));
global_error_explicit_5step = norm(error_adams_moulton_explicit_5step, 2) / sqrt(length(y_true));

% Display or use the global errors as needed
disp(['Global Error (Explicit 2-step): ' num2str(global_error_explicit_2step)]);
disp(['Global Error (Explicit 3-step): ' num2str(global_error_explicit_3step)]);
disp(['Global Error (Explicit 4-step): ' num2str(global_error_explicit_4step)]);
disp(['Global Error (Explicit 5-step): ' num2str(global_error_explicit_5step)]);
% Calculate maximum local errors
% Change the display format to long
format long;
max_error_explicit_2step = max(error_adams_moulton_explicit_2step);
max_error_explicit_3step = max(error_adams_moulton_explicit_3step);
max_error_explicit_4step = max(error_adams_moulton_explicit_4step);
max_error_explicit_5step = max(error_adams_moulton_explicit_5step);

% Display the results
fprintf('Maximum Local Error (Explicit 2-step): %f\n', max_error_explicit_2step);
fprintf('Maximum Local Error (Explicit 3-step): %f\n', max_error_explicit_3step);
fprintf('Maximum Local Error (Explicit 4-step): %f\n', max_error_explicit_4step);
fprintf('Maximum Local Error (Explicit 5-step): %f\n', max_error_explicit_5step);
