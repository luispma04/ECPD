% Illustrates the unconstrained and constrained minimization of the 
% Rosenbrock function.
%
% Matlab Toolboxes: Optimization Toolbox
% Functions: fminunc, fmincon, Rosenbrock
%--------------------------------------------------------------------------

%% 1. Setup Mesh Grid and Function Surface
% Range of independent variables
x1_min = -2;   x1_max = 2;
x2_min = -0.5; x2_max = 3.5;

% Grid discretization
n1 = 100; % Points along x1 axis
n2 = 100; % Points along x2 axis

x_vec1 = linspace(x1_min, x1_max, n1);
x_vec2 = linspace(x2_min, x2_max, n2);
[xx1, xx2] = meshgrid(x_vec1, x_vec2);

% Evaluate Rosenbrock function across the grid
ff = zeros(n1, n2);
for ii = 1:n1
    for jj = 1:n2
        x_point = [xx1(ii,jj); xx2(ii,jj)];
        ff(ii,jj) = Rosenbrock(x_point);
    end
end

%% 2. Plot Level Curves
% Visualization parameters
n_level = 20;     % Number of contour lines
lw      = 1.2;    % Line width
fs      = 14;     % Font size

figure(1);
contour(x_vec1, x_vec2, ff, n_level, 'linewidth', lw);
colorbar;
axis([x1_min x1_max x2_min x2_max]); 
axis square;
hold on;

%% 3. Compute Minima
% Initial guess
x0 = [-1; 1];

% Unconstrained optimization
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton');
x_opt = fminunc(@Rosenbrock, x0, options);

% Constrained optimization (x1 <= 0.5)
% Formulated as A_c * x <= B_c
A_c = [1 0];
B_c = 0.5;
x_opt_constr = fmincon(@Rosenbrock, x0, A_c, B_c);

%% 4. Annotate Plot
% Plot markers
plot(x0(1), x0(2), 'or', 'LineWidth', 1.5);           % Initial guess
plot(x_opt(1), x_opt(2), 'xr', 'LineWidth', 1.5);     % Unconstrained min
plot(x_opt_constr(1), x_opt_constr(2), '*r', 'LineWidth', 1.5); % Constrained min

% Plot constraint boundary (x1 = 0.5)
z2_line = linspace(x2_min, x2_max, 100);
z1_line = 0.5 * ones(size(z2_line));
plot(z1_line, z2_line, 'k', 'LineWidth', 1.5);

% Axis labels
xlabel('x_1', 'FontSize', fs);
ylabel('x_2', 'FontSize', fs);
hold off;

%% 5. 3D Surface Visualization
figure(2);
surf(xx1, xx2, ff); 
hold on;

% Add constraint line projected on the surface
n_c = 200;
z2_3d = linspace(x2_min, x2_max, n_c);
z1_3d = 0.5 * ones(size(z2_3d));

f_constr = zeros(1, n_c);
for kk = 1:n_c
    f_constr(kk) = Rosenbrock([z1_3d(kk); z2_3d(kk)]);
end

plot3(z1_3d, z2_3d, f_constr, 'r', 'LineWidth', 2);
legend('Rosenbrock surface', 'Constraint x_1 = 0.5', 'Location', 'best');

% Label axes
xlabel('x_1', 'FontSize', fs);
ylabel('x_2', 'FontSize', fs);
zlabel('f(x)', 'FontSize', fs);

% Set viewing angle
view([-45, 30]); 
hold off;