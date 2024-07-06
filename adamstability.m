clc; clear all; clf;% Define colors
c = {'b','r','g','m','y','c'};

% Define lines for plotting the axes
x = [0 0]; 
y = [-25 25]; 
K = 'k'; 

% Plot the x and y axes
plot(8*y, x, K), hold on
plot(x, 8*y, K)

% Parameterize the circle
t = linspace(0, 2*pi, 100);
z = exp(1i*t); 
r = 0;

% Plot stability regions for different orders
for i = 1:5
    d = 1 - 1./z;
    for j = 1:i
        r = r + d.^j / j;
    end
    fill(real(r), imag(r), c{i}, 'FaceAlpha', 0.2); % Fill the region inside the stability curve
    plot(real(r), imag(r), c{i}); % Plot the stability curve
    r = 0;
end

% Set up axis limits, make the axis square, and display the grid
axis([-15 35 -25 25]), axis square, grid on

% Add text labels inside each region using gtext (may not work in all MATLAB environments)
gtext('1st Order', 'Color', c{5});
gtext('2nd Order', 'Color', c{4});
gtext('3rd Order', 'Color', c{3});
gtext('4th Order', 'Color', c{2});
gtext('5th Order', 'Color', c{1});
xlabel('Real (z)');ylabel('Im (z)');
% Add title
%title('Backward Differentiation Orders 1-6')
