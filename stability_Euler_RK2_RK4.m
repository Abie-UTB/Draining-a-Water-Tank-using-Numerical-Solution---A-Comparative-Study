% Define colors
clc; clear all; clf;
c = {'b','r','g'};

% Define lines for plotting the axes
x = [0 0]; 
y_ax = [-8 8]; % Define y-axis
K = 'k'; 

% Plot the x and y axes
plot(y_ax, x, K, 'LineWidth', 1.5), hold on
plot(x, y_ax, K, 'LineWidth', 1.5)

% Parameterize the circle
t = linspace(0, 2*pi, 100);
z = exp(1i*t); 

clf
plot(y_ax, x, K, 'LineWidth', 1.5), hold on
plot(x, y_ax, K, 'LineWidth', 1.5)

w = z - 1;
fill(real(w), imag(w), c{1}, 'FaceAlpha', 0.2) % Shade inside stability region for order 1
plot(w, c{1}, 'LineWidth', 2) % Plot the boundary for order 1

% Add text labels for order 1
gtext('Explicit Euler')

for i = 1:3
   w = w - (1 + w + 0.5*w.^2 - z.^2) ./ (1 + w);
end
fill(real(w), imag(w), c{2}, 'FaceAlpha', 0.2) % Shade inside stability region for order 2
plot(w, c{2}, 'LineWidth', 2) % Plot the boundary for order 2

% Add text labels for order 2
gtext('Modified Euler')

for i = 1:4
   w = w - (1 + w + 0.5*w.^2 + w.^3/6 - z.^3) ./ (1 + w + 0.5*w.^2);
end
fill(real(w), imag(w), c{3}, 'FaceAlpha', 0.2) % Shade inside stability region for order 4
plot(w, c{3}, 'LineWidth', 2) % Plot the boundary for order 4

% Add text labels for order 4
gtext('RK4')

axis([-5 2 -3.5 3.5]), axis square, grid on
xlabel('Real (z)');ylabel('Im (z)');
%title('Runge-Kutta Stability Regions')
