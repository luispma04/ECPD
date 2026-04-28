x1min = -2; x1max = 2;
x2min = -1; x2max = 3;
N1 = 100;
N2 = 100;

xv1 = linspace(x1min, x1max, N1);
xv2 = linspace(x2min, x2max, N2);
[xx1, xx2] = meshgrid(xv1, xv2);

for ii = 1:N1
    for jj = 1:N2
        x = [xx1(ii,jj); xx2(ii,jj)];
        ff(ii,jj) = RosenbrockFunction(x);
    end
end

x0 = [-1; 1];

options = optimoptions('fminunc', 'Algorithm', 'quasi-newton');

xopt = fminunc(@RosenbrockFunction, x0, options);

Ac = [1 0];
Bc = 0.5;
xoptconstr = fmincon(@RosenbrockFunction, x0, Ac, Bc);

%--------------------------------------------------------------------------
% Figure 1: level curves with initial estimate, unconstrained and
% constrained minima

Nlevel = 30;
figure(1)
contour(xv1, xv2, ff, Nlevel, 'LineWidth', 1.2)
colorbar
axis([x1min x1max x2min x2max])
axis square
hold on

plot([0.5 0.5], [x2min x2max], 'k--', 'LineWidth', 1.5)

% Plot initial estimate marked with "o"
gg = plot(x0(1), x0(2), 'or');
set(gg, 'LineWidth', 2, 'MarkerSize', 12)

% Plot unconstrained minimum marked with "x"
gg = plot(xopt(1), xopt(2), 'xb');
set(gg, 'LineWidth', 2, 'MarkerSize', 12)

% Plot constrained minimum marked with "*"
gg = plot(xoptconstr(1), xoptconstr(2), '*g');
set(gg, 'LineWidth', 2, 'MarkerSize', 12)

gg = xlabel('x_1');
set(gg, 'FontSize', 14)
gg = ylabel('x_2');
set(gg, 'FontSize', 14)
title('Rosenbrock function - level curves', 'FontSize', 14)

hold off

%--------------------------------------------------------------------------
% Figure 2: 3D view of the Rosenbrock function

figure(2)
surf(xx1, xx2, ff)
shading interp

gg = xlabel('x_1');
set(gg, 'FontSize', 14)
gg = ylabel('x_2');
set(gg, 'FontSize', 14)
gg = zlabel('f(x)');
set(gg, 'FontSize', 14)
title('Rosenbrock function - 3D view', 'FontSize', 14)

fprintf('Unconstrained minimum: x1 = %.4f, x2 = %.4f\n', xopt(1), xopt(2))
fprintf('Constrained minimum:   x1 = %.4f, x2 = %.4f\n', xoptconstr(1), xoptconstr(2))