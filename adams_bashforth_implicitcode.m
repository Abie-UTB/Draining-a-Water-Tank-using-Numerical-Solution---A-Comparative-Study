% Given parameters
g = 9.81;  % Gravitational acceleration
D = 2;     % Diameter of the cylinder
d = 0.2;   % Diameter of the small hole at the bottom of the cylinder
h = 5;     % Step size
t = 0:h:127;  % Time interval to be evaluated
n = length(t) - 1;

% Function to compute the derivative
y_dot = @(x, y) (-sqrt(2 * g) * (d / D)^2 * sqrt(y));

% True solution (for comparison)
y_true = zeros(1, length(t));
y_true(1) = 8;   % True solution for the initial level of water
for i = 1:n
    y_true(i + 1) = (sqrt(8) - sqrt(g/2) * (d / D)^2 * t(i+1)).^2;
end

% Implement Implicit Euler method
y_implicit_euler = zeros(1, length(t));
y_implicit_euler(1) = 8; % Initial water level
for i = 1:n
    % Use fzero to find the next value of y
    y_implicit_euler(i + 1) = fzero(@(y) y - y_implicit_euler(i) - h * y_dot(t(i), y), y_implicit_euler(i));
end

% Implement First-order Adams-Moulton method
y_adams_moulton = zeros(1, length(t));
y_adams_moulton(1) = 8; % Initial water level
for i = 1:n
    % Use fzero to find the next value of y
    y_adams_moulton(i + 1) = fzero(@(y) y - y_adams_moulton(i) - h/2 * (y_dot(t(i+1), y) + y_dot(t(i), y_adams_moulton(i))), y_adams_moulton(i));
end

% Plotting the results and errors for comparison
figure;

% Plotting water level comparison
subplot(2,1,1);
plot(t, y_true, 'k-', t, y_implicit_euler, 'b--', t, y_adams_moulton, 'r-.');
legend('Exact Solution', 'Implicit Euler Method', 'First-order Adams-Moulton Method');
xlabel('Time');
ylabel('Water Level');
title('Water Drainage System: Comparison of Methods');

% Calculate errors
errors_implicit_euler = abs(y_true - y_implicit_euler);
errors_adams_moulton = abs(y_true - y_adams_moulton);

% Plotting error comparison
subplot(2,1,2);
plot(t, errors_implicit_euler, 'b-', t, errors_adams_moulton, 'r--');
legend('Implicit Euler Method', 'First-order Adams-Moulton Method');
xlabel('Time');
ylabel('Error');
title('Error Comparison between Methods');
