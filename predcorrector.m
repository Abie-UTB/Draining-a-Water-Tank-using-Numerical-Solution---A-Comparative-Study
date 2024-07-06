clc;
clear all;

% Constants
g = 9.81;  % Gravitational acceleration
D = 2;     % Diameter of the cylinder
d = 0.2;   % Diameter of the small hole at the bottom of the cylinder
end_time = 127;  % Simulation time

% Time array for h = 1
h1 = 1;    
t1 = 0:h1:end_time;

% Time array for h = 2.5
h2_5 = 2.5;
t2_5 = 0:h2_5:end_time;

% Time array for h = 5
h5 = 5;
t5 = 0:h5:end_time;

% Initialize exact solution for all step sizes
y_exact1 = (sqrt(8) - sqrt(g/2) * (d / D)^2 * t1).^2;
y_exact2_5 = (sqrt(8) - sqrt(g/2) * (d / D)^2 * t2_5).^2;
y_exact5 = (sqrt(8) - sqrt(g/2) * (d / D)^2 * t5).^2;

% Adams-Bashforth-Moulton Method for h = 1
y_abm1 = zeros(size(t1));
y_abm1(1) = y_exact1(1);
y_abm1(2) = y_abm1(1) + h1 * (-sqrt(2 * g) * (d / D)^2 * sqrt(y_abm1(1)));

for i = 2:length(t1)-1
    predictor = y_abm1(i) + h1/2 * (3 * (-sqrt(2 * g) * (d / D)^2 * sqrt(y_abm1(i))) - (-sqrt(2 * g) * (d / D)^2 * sqrt(y_abm1(i-1))));
    y_abm1(i+1) = y_abm1(i) + h1/2 * ((-sqrt(2 * g) * (d / D)^2 * sqrt(predictor)) + (-sqrt(2 * g) * (d / D)^2 * sqrt(y_abm1(i))));
end

% Adams-Bashforth-Moulton Method for h = 2.5
y_abm2_5 = zeros(size(t2_5));
y_abm2_5(1) = y_exact2_5(1);
y_abm2_5(2) = y_abm2_5(1) + h2_5 * (-sqrt(2 * g) * (d / D)^2 * sqrt(y_abm2_5(1)));

for i = 2:length(t2_5)-1
    predictor = y_abm2_5(i) + h2_5/2 * (3 * (-sqrt(2 * g) * (d / D)^2 * sqrt(y_abm2_5(i))) - (-sqrt(2 * g) * (d / D)^2 * sqrt(y_abm2_5(i-1))));
    y_abm2_5(i+1) = y_abm2_5(i) + h2_5/2 * ((-sqrt(2 * g) * (d / D)^2 * sqrt(predictor)) + (-sqrt(2 * g) * (d / D)^2 * sqrt(y_abm2_5(i))));
end

% Adams-Bashforth-Moulton Method for h = 5
y_abm5 = zeros(size(t5));
y_abm5(1) = y_exact5(1);
y_abm5(2) = y_abm5(1) + h5 * (-sqrt(2 * g) * (d / D)^2 * sqrt(y_abm5(1)));

for i = 2:length(t5)-1
    predictor = y_abm5(i) + h5/2 * (3 * (-sqrt(2 * g) * (d / D)^2 * sqrt(y_abm5(i))) - (-sqrt(2 * g) * (d / D)^2 * sqrt(y_abm5(i-1))));
    y_abm5(i+1) = y_abm5(i) + h5/2 * ((-sqrt(2 * g) * (d / D)^2 * sqrt(predictor)) + (-sqrt(2 * g) * (d / D)^2 * sqrt(y_abm5(i))));
end

% Plot the water levels for all step sizes
figure;
subplot(2,1,1);
plot(t1, y_abm1, 'LineWidth', 2, 'DisplayName', 'Adams-Bashforth-Moulton (h = 1)');
hold on;
plot(t2_5, y_abm2_5, '--', 'LineWidth', 2, 'DisplayName', 'Adams-Bashforth-Moulton (h = 2.5)');
plot(t5, y_abm5, '-.', 'LineWidth', 2, 'DisplayName', 'Adams-Bashforth-Moulton (h = 5)');
xlabel('Time [sec]');
ylabel('Water Level [m]');
%title('Comparison of Different Step Sizes');
legend('Location', 'Best');
grid on;

% Plot the exact solution for comparison
%plot(t1, y_exact1, '-.', 'LineWidth', 2, 'DisplayName', 'Exact Solution (h = 1)');
%plot(t2_5, y_exact2_5, '-.', 'LineWidth', 2, 'DisplayName', 'Exact Solution (h = 2.5)');
plot(t5, y_exact5, '-.', 'LineWidth', 2, 'DisplayName', 'Exact Solution');
legend('Location', 'Best');

% Plot the errors for all step sizes
subplot(2,1,2);
error_abm1 = abs(y_exact1 - y_abm1);
error_abm2_5 = abs(y_exact2_5 - y_abm2_5);
error_abm5 = abs(y_exact5 - y_abm5);
plot(t1, error_abm1, 'LineWidth', 2, 'DisplayName', 'Error (h = 1)');
hold on;
plot(t2_5, error_abm2_5, '--', 'LineWidth', 2, 'DisplayName', 'Error (h = 2.5)');
plot(t5, error_abm5, '-.', 'LineWidth', 2, 'DisplayName', 'Error (h = 5)');
xlabel('Time [sec]');
ylabel('Local Error');
%title('Comparison of Errors for Different Step Sizes');
legend('Location', 'Best');
grid on;

% Calculate and display global error for each case
global_error1 = norm(error_abm1, 2) / sqrt(length(y_exact1));
max_error1 = max(error_abm1);
disp(['Global Error (h = 1): ' num2str(global_error1)]);
fprintf('Maximum Error (h = 1): %f\n', max_error1);

global_error2_5 = norm(error_abm2_5, 2) / sqrt(length(y_exact2_5));
max_error2_5 = max(error_abm2_5);
disp(['Global Error (h = 2.5): ' num2str(global_error2_5)]);
fprintf('Maximum Error (h = 2.5): %f\n', max_error2_5);

global_error5 = norm(error_abm5, 2) / sqrt(length(y_exact5));
max_error5 = max(error_abm5);
disp(['Global Error (h = 5): ' num2str(global_error5)]);
fprintf('Maximum Error (h = 5): %f\n', max_error5);
