% Simulation of the TCLab linear model previously identified
%
% Loads the model identified in the TCLab_identification script, creates
% the h1 and T1C functions that mimick the TCLab interface, and performs a
% simulation starting at ambient temperature.
% You will be developing and testing your MPC controller and Kalman filter
% in this simulation environment.
%
% Q2: Unconstrained closed-loop MPC — study effect of H (fixed R) and R (fixed H)
% Q3: Constrained closed-loop MPC   — add control limits [0,100]%, tune R for saturation
%
% Afonso Botelho and J. Miranda Lemos, IST, May 2023
%__________________________________________________________________________

% Initialization
clear all
close all
clc

% ── Paths ─────────────────────────────────────────────────────────────────
base = fileparts(mfilename('fullpath'));   % folder where this script lives
addpath(fullfile(base, '../P3'));          % singleheater_model.mat
addpath(fullfile(base, '../P4'));          % mpc_solve.m

% Load model
load('singleheater_model.mat','A','B','C','Ke','e_var','y_ss','u_ss','Ts');
n = size(A,1);
e_std = sqrt(e_var); %choose to activate or not noise, guide recommends to keep it off to test the controller in ideal conditions

% Build the functions for applying the control and reading the temperature,
% mimicking the TCLab interface
x_ss = [eye(n)-A; C]\[B*u_ss; y_ss];
c1 = ((eye(n)-A)*x_ss - B*u_ss);
c2 = (y_ss - C*x_ss);
h1  = @(x,u) A*x + B*u + Ke*e_std*randn + c1;  % apply control
T1C = @(x) C*x + e_std*randn + c2;              % read temperature

% Simulation parameters
T = 4000;       % Experiment duration [s]
N = T/Ts;       % Number of samples to collect

% Initial conditions (start at ambient temperature, i.e. equilibrium for u = 0)
Dx0Dy0 = [eye(n)-A, zeros(n,1); C, -1]\[-B*u_ss; 0];
Dx0 = Dx0Dy0(1:n);
x0  = Dx0 + x_ss;

%% ═══════════════════════════════════════════════════════════════════════
%  Q2 — UNCONSTRAINED MPC
%  No lb/ub passed to mpc_solve. No clipping of u: the unconstrained
%  behaviour is intentionally observed (may exceed [0,100]%).
% ═══════════════════════════════════════════════════════════════════════

%% ── Q2 Study 1: Effect of H (R fixed) ───────────────────────────────────
% R is fixed at a neutral value so only H changes and its effect is isolated.
R_fixed = 1;
H_list  = [3, 5, 10, 20, 50];      % H >= n = 3 (model order)

colors_H = lines(length(H_list));

figure('Units','normalized','Position',[0.05 0.5 0.55 0.4])
subplot(2,1,1), hold on, grid on
title(sprintf('Q2 — Effect of H on \\Deltay  (R = %.1f fixed)', R_fixed))
ylabel('\Delta y [°C]'), xlabel('Time [s]')
yline(0,'k--','LineWidth',1)
subplot(2,1,2), hold on, grid on
title('Corresponding \Deltau  (unconstrained)')
ylabel('\Delta u [%]'), xlabel('Time [s]')

for i = 1:length(H_list)
    H_i = H_list(i);
    fprintf('Q2 Study 1 — H=%d, R=%.1f ...\n', H_i, R_fixed)

    t  = nan(1,N);  y  = nan(1,N);
    Dy = nan(1,N);  Du = nan(1,N);
    u  = nan(1,N);
    x  = nan(n,N+1);
    x(:,1) = x0;

    for k = 1:N
        t(k)    = (k-1)*Ts;
        y(k)    = T1C(x(:,k));
        Dy(k)   = y(k) - y_ss;
        Dx_k    = x(:,k) - x_ss;

        % Unconstrained: no lb/ub, no clipping
        Du_k    = mpc_solve(Dx_k, H_i, R_fixed, A, B, C);
        u(k)    = u_ss + Du_k;
        Du(k)   = u(k) - u_ss;

        x(:,k+1) = h1(x(:,k), u(k));
    end

    subplot(2,1,1)
    plot(t, Dy, 'Color', colors_H(i,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('H = %d', H_i))
    subplot(2,1,2)
    stairs(t, Du, 'Color', colors_H(i,:), 'LineWidth', 1.5, ...
           'DisplayName', sprintf('H = %d', H_i))
end
subplot(2,1,1), legend('Location','best')
subplot(2,1,2), legend('Location','best')

%% ── Q2 Study 2: Effect of R (H fixed at chosen value) ───────────────────
% After analysing Study 1, set H_chosen to the value where the response
% has converged. R is then varied to study its effect on aggressiveness.
H_chosen = 10;      % <-- update after analysing Study 1
R_list   = [0.1, 1, 10, 100];

colors_R = lines(length(R_list));

figure('Units','normalized','Position',[0.55 0.5 0.55 0.4])
subplot(2,1,1), hold on, grid on
title(sprintf('Q2 — Effect of R on \\Deltay  (H = %d fixed)', H_chosen))
ylabel('\Delta y [°C]'), xlabel('Time [s]')
yline(0,'k--','LineWidth',1)
subplot(2,1,2), hold on, grid on
title('Corresponding \Deltau  (unconstrained)')
ylabel('\Delta u [%]'), xlabel('Time [s]')

for i = 1:length(R_list)
    R_i = R_list(i);
    fprintf('Q2 Study 2 — H=%d, R=%.2f ...\n', H_chosen, R_i)

    t  = nan(1,N);  y  = nan(1,N);
    Dy = nan(1,N);  Du = nan(1,N);
    u  = nan(1,N);
    x  = nan(n,N+1);
    x(:,1) = x0;

    for k = 1:N
        t(k)    = (k-1)*Ts;
        y(k)    = T1C(x(:,k));
        Dy(k)   = y(k) - y_ss;
        Dx_k    = x(:,k) - x_ss;

        % Unconstrained: no lb/ub, no clipping
        Du_k    = mpc_solve(Dx_k, H_chosen, R_i, A, B, C);
        u(k)    = u_ss + Du_k;
        Du(k)   = u(k) - u_ss;

        x(:,k+1) = h1(x(:,k), u(k));
    end

    subplot(2,1,1)
    plot(t, Dy, 'Color', colors_R(i,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('R = %.1f', R_i))
    subplot(2,1,2)
    stairs(t, Du, 'Color', colors_R(i,:), 'LineWidth', 1.5, ...
           'DisplayName', sprintf('R = %.1f', R_i))
end
subplot(2,1,1), legend('Location','best')
subplot(2,1,2), legend('Location','best')

%% ═══════════════════════════════════════════════════════════════════════
%  Q3 — CONSTRAINED MPC
%  H is fixed from Q2. lb and ub encode [0,100]% in incremental form.
%  No clipping after mpc_solve: the solver guarantees the bounds.
%  R is tuned until the control saturates at the limits.
% ═══════════════════════════════════════════════════════════════════════

% Control bounds in incremental form (fixed, computed once before the loop)
lb = -u_ss       * ones(H_chosen, 1);   % Du >= -u_ss  =>  u >= 0
ub = (100-u_ss)  * ones(H_chosen, 1);   % Du <= 100-u_ss  =>  u <= 100

% Choose R small enough that the control saturates
% Decrease R until saturation is visible in the plot
R_Q3 = 0.1;    % <-- tune this: try 1, 0.5, 0.1, 0.01

fprintf('Q3 — H=%d, R=%.3f (constrained) ...\n', H_chosen, R_Q3)

t  = nan(1,N);  y  = nan(1,N);
Dy = nan(1,N);  Du = nan(1,N);
u  = nan(1,N);
x  = nan(n,N+1);
x(:,1) = x0;

for k = 1:N
    t(k)    = (k-1)*Ts;
    y(k)    = T1C(x(:,k));
    Dy(k)   = y(k) - y_ss;
    Dx_k    = x(:,k) - x_ss;

    % Constrained: lb and ub passed, no clipping needed
    Du_k    = mpc_solve(Dx_k, H_chosen, R_Q3, A, B, C, lb, ub);
    u(k)    = u_ss + Du_k;
    Du(k)   = u(k) - u_ss;

    x(:,k+1) = h1(x(:,k), u(k));
end

% Plot absolute variables
figure('Units','normalized','Position',[0.2 0.05 0.3 0.4])
subplot(2,1,1), hold on, grid on
title(sprintf('Q3 — Absolute  (H=%d, R=%.3f, constrained)', H_chosen, R_Q3))
plot(t, y, '.', 'MarkerSize', 5)
yl = yline(y_ss,'k--');
xlabel('Time [s]'), ylabel('y [°C]')
legend(yl,'$\bar{y}$','Interpreter','latex','Location','best')
subplot(2,1,2), hold on, grid on
stairs(t, u, 'LineWidth', 2)
yl = yline(u_ss,'k--');
yline(0,  'r--')
yline(100,'r--')
xlabel('Time [s]'), ylabel('u [%]')
legend(yl,'$\bar{u}$','Interpreter','latex','Location','best')

% Plot incremental variables
figure('Units','normalized','Position',[0.5 0.05 0.3 0.4])
subplot(2,1,1), hold on, grid on
title(sprintf('Q3 — Incremental  (H=%d, R=%.3f, constrained)', H_chosen, R_Q3))
plot(t, Dy, '.', 'MarkerSize', 5)
yline(0,'k--')
xlabel('Time [s]'), ylabel('\Delta y [°C]')
subplot(2,1,2), hold on, grid on
stairs(t, Du, 'LineWidth', 2)
yline(-u_ss,     'r--')
yline(100-u_ss,  'r--')
xlabel('Time [s]'), ylabel('\Delta u [%]')

%--------------------------------------------------------------------------
% End of File