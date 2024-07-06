clc; 
clear all; 
clf;

% Given parameters
g = 9.81;  % Gravitational acceleration
D = 2;     % Diameter of the cylinder
d = 0.2;   % Diameter of the small hole at the bottom of the cylinder
h = 5;     % Step size
steps = round(127.7 / h); % Number of steps based on the final time

% Initialize arrays
t = zeros(steps+1, 1);
y = zeros(steps+1, 1);

% Initial condition
t(1) = 0;
y(1) = 8; % Initial water level

% Function to compute y_dot = dy/dt
y_dot = @(y) (-sqrt(2 * g) * (d / D)^2 * sqrt(y));

% True solution (for comparison)
y_true = zeros(steps+1, 1);
y_true(1) = 8;

% Compute true solution
for i = 1:steps
    t(i + 1) = t(i) + h;
    y_true(i + 1) = (sqrt(8) - sqrt(g/2) * (d / D)^2 * t(i+1)).^2;
end
figure;
subplot(2,1,1);
plot(t, y_true, '-b', 'LineWidth', 1.5);
hold on
% Compute using Euler's method (1st order Adams-Bashforth-Moulton)
for k = 2:3
     t(k) = t(k-1) + h;
    y(k) = y(k-1) + h * y_dot(y(k-1));
end

% Adams-Bashforth-Moulton method (2nd order)
for k = 3:steps
    t(k) = t(k-1) + h;
    % Adams-Bashforth predictor
    p = y(k-1) + (h/2) * (3*y_dot(y(k-1)) - y_dot(y(k-2)));
    
    % Adams-Moulton corrector
    y(k+1) = y(k) + (h/2) * (y_dot(p) + y_dot(y(k-1)));
end
err1 = abs(y_true - y);
plot(t, y, '-y', 'LineWidth', 1.5);
% Compute the first three steps using a different method (e.g., Euler's method)
for k = 2:4
    t(k) = t(k-1) + h;
    y(k) = y(k-1) + h * y_dot(y(k-1));
end

% Adams-Bashforth-Moulton method (3rd order)
for k = 4:steps
    t(k+1) = t(k) + h;
   
    % Predictor step
    p = y(k) + (h/12) * (23*y_dot(y(k)) - 16*y_dot(y(k-1)) + 5*y_dot(y(k-2)));
    
    % Corrector step
    y(k+1) = y(k) + (h/12) * (5*y_dot(p) + 8*y_dot(y(k)) - y_dot(y(k-1)));
end
err2 = abs(y_true - y);
plot(t, y, '-g', 'LineWidth', 1.5);
% Compute the first four steps using a different method (e.g., Euler's method)
for k = 2:5
    t(k) = t(k-1) + h;
    y(k) = y(k-1) + h * y_dot(y(k-1));
end

% Adams-Bashforth-Moulton method (4th order)
for k = 5:steps
    t(k+1) = t(k) + h;
   
    % Predictor step
    p = y(k) + (h/24) * (55*y_dot(y(k)) - 59*y_dot(y(k-1)) + 37*y_dot(y(k-2)) - 9*y_dot(y(k-3)));
    
    % Corrector step
    y(k+1) = y(k) + (h/24) * (9*y_dot(p) + 19*y_dot(y(k)) - 5*y_dot(y(k-1)) + y_dot(y(k-2)));
end
err3 = abs(y_true - y);
plot(t, y, '-r', 'LineWidth', 1.5);
% Compute the first four steps using a different method (e.g., Euler's method)
for k = 2:6
    t(k) = t(k-1) + h;
    y(k) = y(k-1) + h * y_dot(y(k-1));
end

% Adams-Bashforth-Moulton method (5th order)
for k = 6:steps
    t(k+1) = t(k) + h;
   
    % Predictor step
    p = y(k) + (h/720) * (1901*y_dot(y(k)) - 2774*y_dot(y(k-1)) + 2616*y_dot(y(k-2)) - 1274*y_dot(y(k-3)) + 251*y_dot(y(k-4)));
    
    % Corrector step
    y(k+1) = y(k) + (h/720) * (251*y_dot(p) + 646*y_dot(y(k)) - 264*y_dot(y(k-1)) + 106*y_dot(y(k-2)) - 19*y_dot(y(k-3)));
end
err4 = abs(y_true - y);
plot(t, y, '-k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Water Level [m]');
%title('Water Drainage System (Adams-Bashforth-Moulton Method)');
legend('True Value','2nd Adams-Bashforth-Moulton Method','3rd Adams-Bashforth-Moulton Method','4th Adams-Bashforth-Moulton Method','5th Adams-Bashforth-Moulton Method');
grid on;

% Plotting the errors
subplot(2,1,2);
plot(t, err1, '-g', 'LineWidth', 1.5);
hold on;
plot(t, err2, '-b', 'LineWidth', 1.5);
plot(t, err3, '-r', 'LineWidth', 1.5);
plot(t, err4, '-k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Local Error [m]');
%title('Local Error Comparison');
legend('2nd Order', '3rd Order', '4th Order', '5th Order');
grid on;
