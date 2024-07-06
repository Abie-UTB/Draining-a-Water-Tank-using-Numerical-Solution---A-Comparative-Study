% Program to plot stability plots for Adam Bashforth methods

% Setting grid points
axisbox = [-3.5 0 -1 0];
xa = axisbox(1); xb = axisbox(2); ya = axisbox(3); yb = axisbox(4);
npts  = 60;
theta = linspace(0, 2*pi, 2*npts+1);
z = exp(1i*theta);

% 1st order Adam Bashforth
nu1=(z - 1);

% 2nd order Adam Bashforth
nu2 = (z.^2 - z)./((3*z - 1)/2);

% 3rd order Adam Bashforth
nu3 = (z.^3 - z.^2)./((5-16*z + 23*z.^2)/12);

% 4th order Adam Bashforth
nu4 = (z.^4 -z.^3)./((55*z.^3 -59*z.^2 +37*z -9)/24);
nu5 = (z.^5 -z.^4)./((1901*z.^4 -2774*z.^3 +2616*z.^2 -1274*z + 251)/720);

% chop off intersecting loops for 4th order AB
for k = 1: length(nu5)-1
    z = [real(nu5); imag(nu5)];
    iloop = 0;
    for j = k+2 : length(nu5)-1
        % find first place where z(k) intersects the rest of the path
        lam = inv([z(:,k)-z(:,k+1), z(:,j+1)-z(:,j)])*(z(:,j+1)-z(:,k+1));
        if (lam >=0 & lam <=1)
            iloop = k+1:j;
            zint = lam(1)*z(:,k)+(1-lam(1))*z(:,k+1);
            break
        end
    end
    if (iloop ~= 0)
        % chop out nu4(iloop) and replace with single interpolated value zint
        zcp = complex(zint(1),zint(2));
        nu5(iloop(1)) = zcp;
        iloop(1) = [];
        nu5(iloop) = [];
    end
end

figure();
clf;
hold on;
fill(real(nu1), imag(nu1), 'm', 'FaceAlpha', 0.2); % Fill with magenta color for 1st order
fill(real(nu2), imag(nu2), 'r', 'FaceAlpha', 0.2); % Fill with red color for 2nd order
fill(real(nu3), imag(nu3), 'b', 'FaceAlpha', 0.2); % Fill with blue color for 3rd order
fill(real(nu4), imag(nu4), 'g', 'FaceAlpha', 0.2); % Fill with green color for 4th order
fill(real(nu5), imag(nu5), 'y', 'FaceAlpha', 0.2); % Fill with yellow color for 5th order
plot([xa xb],[0 0],'k-','LineWidth',0.5); % Very thin x-axis line
plot([0 0],[ya yb],'k-','LineWidth',0.5); % Very thin y-axis line
xlim([-3.5 0]);  % Limiting the x-axis range
legend({'Adams-Bashforth (1st Order)','Adams-Bashforth (2nd Order)','Adams-Bashforth (3rd Order)','Adams-Bashforth (4th Order)','Adams-Bashforth (5th Order)'},'FontSize',12);
set(gca,'FontSize',15);
xlabel('Real (z)');ylabel('Im (z)');
grid on;
