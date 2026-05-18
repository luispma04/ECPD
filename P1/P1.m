% Illustrates the unconstrained and constrained minimization of the 
% Rosenbrock function following the structure of ProbBasic.m.
%
% Matlab sw required: optimization toolbox
% Functions called
%    fminunc    - Matlab function for unconstrained minimization
%    fmincon    - Matlab function for constrained minimization
%    Rosenbrock - user defined function to be minimized
%
% IST, MEEC, Distributed Predictive Control and Estimation
% Afonso Botelho, Joao Miranda Lemos, 2025
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plots function level curves

% Range of independent variables to consider (Adjusted for Rosenbrock)
x1min = -2;
x1max = 2;
x2min = -0.5;
x2max = 3.5;

% Number of intervals in the mesh grid
N1 = 100;
N2 = 100;

xv1 = linspace(x1min, x1max, N1);
xv2 = linspace(x2min, x2max, N2);
[xx1, xx2] = meshgrid(xv1, xv2);

% Computes the function at the different points of the mesh grid
for ii = 1:N1
    for jj = 1:N2
        x = [xx1(ii,jj); xx2(ii,jj)];
        ff(ii,jj) = Rosenbrock(x);
    end
end

% Plots the level curves using the Matlab function contour
Nlevel = 20;  
LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';
figure(1), contour(xv1, xv2, ff, Nlevel, LW, 1.2), colorbar
axis([x1min x1max x2min x2max]), axis square
hold on

%--------------------------------------------------------------------------
% Compute the minimum

% Initial estimate of the minimum
x0 = [-1; 1];

% Define the options to be used with the fminunc solver
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton');

% Uses the solver fminunc to compute the unconstrained minimum
xopt = fminunc(@Rosenbrock, x0, options);

%--------------------------------------------------------------------------
% Computes the constrained minimum associated to the constraint x1 <= 0.5
% Written as Ac*x <= Bc

Ac = [1 0];
Bc = 0.5;
xoptconstr = fmincon(@Rosenbrock, x0, Ac, Bc);

%--------------------------------------------------------------------------
% Plots markers (using the same red color style as ProbBasic.m)

% Initial point as a red circle
gg = plot(x0(1), x0(2), 'or');
set(gg, 'Linewidth', 1.5);

% Unconstrained minimum as a red cross
gg = plot(xopt(1), xopt(2), 'xr');
set(gg, 'Linewidth', 1.5);

% Constrained minimum as a red star
gg = plot(xoptconstr(1), xoptconstr(2), '*r');
set(gg, 'Linewidth', 1.5);

% Plots the constraint boundary (x1 = 0.5)
z2c = linspace(x2min, x2max, 100);
z1c = 0.5 * ones(size(z2c));
gg = plot(z1c, z2c, 'k');
set(gg, 'Linewidth', 1.5);

% Identifies axis
gg = xlabel('x_1'); set(gg, 'FontSize', 14);
gg = ylabel('x_2'); set(gg, 'FontSize', 14);

hold off

%--------------------------------------------------------------------------
% Plots the 3d view of the function with the square grid

figure(2)
surf(xx1, xx2, ff); % Standard surf creates the 'squares grid'
hold on

% --- NEW: Create and plot the 3D constraint plane ---
% Determine the vertical bounds based on the evaluated function
z_min = min(ff(:)); 
z_max = max(ff(:)); 

% Create a meshgrid for the plane spanning x2 and z
[P_x2, P_z] = meshgrid([x2min, x2max], [z_min, z_max]);

% The constraint is fixed at x1 = 0.5
P_x1 = 0.5 * ones(size(P_x2));

% Plot the plane (red, 50% transparent, no grid lines on the plane)
surf(P_x1, P_x2, P_z, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Identifies axis
gg = xlabel('x_1'); set(gg, 'FontSize', 14);
gg = ylabel('x_2'); set(gg, 'FontSize', 14);
gg = zlabel('f(x)'); set(gg, 'FontSize', 14);

% Optional: Set a good default 3D viewing angle to see the intersection
view([-45, 30]); 

hold off
%--------------------------------------------------------------------------

