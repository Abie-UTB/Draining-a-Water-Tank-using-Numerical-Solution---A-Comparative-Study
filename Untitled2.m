clear all; close all; clc; % Constants
g = 9.81;  % Gravitational acceleration
D = 2;     % Diameter of the cylinder
d = 0.2;   % Diameter of the small hole at the bottom of the cylinder
h = 5;     % Step size
tend = 127; % End time

% Initialization
t = 0:h:tend;  % Time interval to be evaluated
n = length(t) - 1;
y_true = zeros(1, length(t));  % True solution
y_true(1) = 8;   % True solution for the initial level of water
y_adams_bashforth_1 = zeros(1, length(t));  % Water level for Adams-Bashforth method (1st order)
y_adams_bashforth_2 = zeros(1, length(t));  % Water level for Adams-Bashforth method (2nd order)
y_adams_bashforth_3 = zeros(1, length(t));  % Water level for Adams-Bashforth method (3rd order)
y_adams_bashforth_4 = zeros(1, length(t));  % Water level for Adams-Bashforth method (4th order)
y_adams_bashforth_5 = zeros(1, length(t));  % Water level for Adams-Bashforth method (5th order)

% Function to compute y_dot = dy/dt
y_dot = @(t, y) -sqrt(2 * g) * (d / D)^2 * sqrt(y);

% True solution (for comparison)
for i = 1:n
    y_true(i + 1) = (sqrt(8) - sqrt(g/2) * (d / D)^2 * t(i+1)).^2;
end

% Apply First-Order Adams-Bashforth method
y_adams_bashforth_1(1) = y_true(1);
for i = 1:n
    y_adams_bashforth_1(i + 1) = y_adams_bashforth_1(i) + h * y_dot(t(i), y_adams_bashforth_1(i));
end

% Apply Second-Order Adams-Bashforth method
y_adams_bashforth_2(1:2) = y_true(1:2);

for i = 2:n
    y_adams_bashforth_2(i+1) = y_adams_bashforth_2(i) + h/2 * (3 * y_dot(t(i), y_adams_bashforth_2(i)) - 1 * y_dot(t(i-1), y_adams_bashforth_2(i-1) ));
end

% Apply Third-Order Adams-Bashforth method
y_adams_bashforth_3(1:3) = y_true(1:3);
for i = 3:n
    y_adams_bashforth_3(i+1) = y_adams_bashforth_3(i) + h/12 * (23 * y_dot(t(i), y_adams_bashforth_3(i)) - 16 * y_dot(t(i-1), y_adams_bashforth_3(i-1)) + 5 * y_dot(t(i-2), y_adams_bashforth_3(i-2)) );
end

y_adams_bashforth_4(1:4) = y_true(1:4); % Initial condition
for i = 4:n
    y_adams_bashforth_4(i+1) = y_adams_bashforth_4(i) + h/24 * (55 * y_dot(t(i), y_adams_bashforth_4(i)) - 59 * y_dot(t(i-1), y_adams_bashforth_4(i-1)) + 37 * y_dot(t(i-2), y_adams_bashforth_4(i-2)) - 9 * y_dot(t(i-3), y_adams_bashforth_4(i-3)));
end

% Apply Fifth-Order Adams-Bashforth method
y_adams_bashforth_5(1:5) = y_true(1:5);
for i = 5:n
    y_adams_bashforth_5(i+1) = y_adams_bashforth_5(i) + h/720 * (1901 * y_dot(t(i), y_adams_bashforth_5(i)) - 2774 * y_dot(t(i-1), y_adams_bashforth_5(i-1)) + 2616 * y_dot(t(i-2), y_adams_bashforth_5(i-2)) - 1274 * y_dot(t(i-3), y_adams_bashforth_5(i-3))+251 * y_dot(t(i-3), y_adams_bashforth_5(i-4)));
end
figure(1)
% Plotting the results
plot(t, y_true, 'k-', 'LineWidth', 2);
hold on;
% plot(t, y_adams_bashforth_1, '-o', 'LineWidth', 1.5);
plot(t, y_adams_bashforth_2, '-s', 'LineWidth', 1.5);
plot(t, y_adams_bashforth_3, '-^', 'LineWidth', 1.5);
plot(t, y_adams_bashforth_4, '-d', 'LineWidth', 1.5);
plot(t, y_adams_bashforth_5, '-+', 'LineWidth', 1.5);
hold off;

% Add labels and legend
title('Water Drainage Using Different Orders of Adams-Bashforth Methods');
xlabel('Time (sec)');
ylabel('Water Level (m)');
legend('True Solution', '2nd Order', '3rd Order', '4th Order', '5th Order');
grid on;

% Calculate errors
error_adams_bashforth_1 = abs(y_true - y_adams_bashforth_1);
error_adams_bashforth_2 = abs(y_true - y_adams_bashforth_2);
error_adams_bashforth_3 = abs(y_true - y_adams_bashforth_3);
error_adams_bashforth_4 = abs(y_true - y_adams_bashforth_4);
error_adams_bashforth_5 = abs(y_true - y_adams_bashforth_5);

% Plotting the errors
figure(2)
% plot(t, error_adams_bashforth_1, '-o', 'LineWidth', 1.5);

plot(t, error_adams_bashforth_2, '-s', 'LineWidth', 1.5);
hold on;
plot(t, error_adams_bashforth_3, '-^', 'LineWidth', 1.5);
plot(t, error_adams_bashforth_4, '-d', 'LineWidth', 1.5);
plot(t, error_adams_bashforth_5, '-+', 'LineWidth', 1.5);
hold off;

% Add labels and legend for errors plot
title('Error Comparison of Adams-Bashforth Methods');
xlabel('Time (sec)');
ylabel('Absolute Error');
legend('2nd Order', '3rd Order', '4th Order', '5th Order');
grid on;