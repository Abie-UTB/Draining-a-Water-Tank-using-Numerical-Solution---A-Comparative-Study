% Define the complex range for plotting the stability region
[x, y] = meshgrid(-3:0.01:3, -3:0.01:3);
z = x + 1i*y;

% Stability region for Forward Euler (explicit)
FE_stability = 1 + z;

% Stability region for Backward Euler (implicit)
BE_stability = 1 ./ (1 - z);

% Plotting
figure;

subplot(1, 2, 1);
contourf(x, y, abs(FE_stability), [0, 1], 'LineColor', 'none'); % Plot stability regions for FE
colormap([1 1 1; 1 1 0]); % Set colormap to white and yellow
hold on;
contour(x, y, abs(FE_stability), [0, 1], 'LineColor', 'k'); % Add outlines for unstable region
xlabel('Re(z)');
ylabel('Im(z)');
title('Stability region for Explicit Euler');
axis square;
hold off;

subplot(1, 2, 2);
contourf(x, y, abs(BE_stability), [0, 1], 'LineColor', 'none'); % Plot stability regions for BE
colormap([1 1 0; 1 1 1]); % Set colormap to white and yellow
hold on;
contour(x, y, abs(BE_stability), [0, 1], 'LineColor', 'k'); % Add outlines for unstable region
xlabel('Re(z)');
ylabel('Im(z)');
title('Stability region for Implicit Euler');
axis square;
hold off;
