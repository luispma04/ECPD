%% Closed-Loop Experiment for Data Collection in TCLab (P5)
%--------------------------------------------------------------------------
% Applies the Kalman filter and MPC designed in P4 to the REAL plant.
% Matches the specification of P4 Question 7: the reference commutes across
% multiple temperature levels with feedforward and disturbance estimation 
% (zero-offset MPC).
%
% Recommended workflow: test the Kalman filter first with an open-loop
% profile (set USE_MPC = false), then enable the MPC (USE_MPC = true).
%
% Required variables in .mat: A, B, C, Ke, e_var, y_ss, u_ss, Ts
%
% Matlab Toolboxes: Control System Toolbox, Optimization Toolbox
% Functions: tclab, mpc_solve, mpc_solve_sparse_regularized
%
% IST - MEEC - Distributed Predictive Control and Estimation
% Afonso Botelho, Joao Miranda Lemos, 2025
%--------------------------------------------------------------------------

%% 1. Initialization and Setup

clear; clc; close all;

% Add required folders to path
base_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(base_dir, '..', 'TCLab_software_test'));
addpath(fullfile(base_dir, '..', 'P4'));

% Initialize TCLab hardware
tclab;

% Load identified model parameters
load(fullfile(base_dir, '..', 'P4', 'P4_4', 'singleheater_model2D.mat'), ...
    'A', 'B', 'C', 'Ke', 'e_var', 'y_ss', 'u_ss', 'Ts');

n = size(A, 1); % State vector dimension

%% 2. Experiment Parameters

% Same specification as P4 Q7: reference levels held for seg_time each.
% r_levels = [50 40 65 45];                    % absolute reference [degC]
% r_levels = [50];
% seg_time = 400;                              % seconds held at each level

r_levels = [50 40 65 45 30 80 50 65 35 20 50]; % absolute reference [degC]
seg_time = 600;                                % seconds held at each level

% r_levels = [50 60 70 80 90 30 60 30 70 30 80 30 90]; 
% seg_time = 400;                              

T_total = seg_time * numel(r_levels);          % Total experiment duration [s]
N       = round(T_total / Ts);                 % Number of samples to collect

%% 3. Controller / Observer Options

USE_MPC    = true; % false = open-loop test of the Kalman filter only
mode       = 1;    % 0 = dense mpc_solve, 1 = sparse mpc_solve_regularized
const_type = 1;    % 0 = hard output constraint, 1 = soft

% MPC tuning (fixed from P4)
H = 50;            % Prediction horizon
R = 0.04;          % Control weight

% Safety constraints
y_max  = 55;       % Hard safety cap [degC]
margin = 0.5;      % Clamp reference this far BELOW the cap [degC]

%% 4. Kalman Filter Design

% Augmented state matrices (xd = [Dx; d])
Ad = [A,          B;
      zeros(1,n), 1];
Bd = [B; 0];
Cd = [C, 0];

Q_E     = Ke * e_var * Ke'; % Process-noise covariance (identified)
R_E     = e_var;            % Measurement-noise covariance (identified)
delta_e = 0.1;              % Disturbance-state covariance (Tuning knob)

Q_Ed = [Q_E,         zeros(n, 1);
        zeros(1, n), delta_e];

% Compute steady-state Kalman gain
L = dlqe(Ad, eye(n+1), Cd, Q_Ed, R_E);

% Initial conditions (Start at ambient temperature, i.e., equilibrium for u=0)
Dx0Dy0 = [eye(n) - A, zeros(n, 1); 
          C,          -1] \ [-B * u_ss; 0];
      
Dx0    = Dx0Dy0(1:n); % Used to initialize the filter

%% 5. Reference Trajectory Generation

% Build the absolute reference signal r(k) segmented by time
r = zeros(1, N);
for s = 1:numel(r_levels)
    idx = (1:N) > (s-1)*seg_time/Ts & (1:N) <= s*seg_time/Ts;
    r(idx) = r_levels(s);
end

Dr = r - y_ss; % Incremental reference

% Steady-state matrix (constant -- reused every step in the loop)
Mss = [eye(n) - A, -B;
       C,           0];

%% 6. Real-Time Plot Setup and Pre-allocation

% Real-time plot flag (true = real-time updating, false = print to console)
rt_plot = true;

if rt_plot
    figure;
    drawnow;
end

% Pre-allocate signal arrays
t        = nan(1, N);
u        = zeros(1, N);
y        = zeros(1, N);
Dy       = nan(1, N);
Du       = zeros(1, N);
xd_est   = nan(n+1, N+1);
exitflag = nan(1, N);

% Kalman filter initialization
xd_est(:, 1) = [Dx0; 0]; 

% Date string for uniquely saving results
time_str = char(datetime('now', 'Format', 'yyMMdd_HHmmSS'));

%% 7. Hardware Execution Loop

% Signal the start of the experiment by lighting the LED
led(1);
disp('Temperature test started.');

for k = 1:N
    tic;

    % Compute elapsed time
    t(k) = (k - 1) * Ts;

    % Read the sensor temperatures
    y(1, k) = T1C();

    % Compute incremental variables
    Dy(:, k) = y(:, k) - y_ss;

    % --- Kalman Filter CORRECTION Step ---
    % xd_est(:, k) currently holds the PREDICTION made last iteration.
    % Correct it using the new measurement Dy(:, k).
    xd_est(:, k) = xd_est(:, k) + L * (Dy(:, k) - Cd * xd_est(:, k));

    Dx_hat = xd_est(1:n, k); % Estimated incremental state
    d_hat  = xd_est(end, k); % Estimated input disturbance

    if USE_MPC
        % --- Feedforward: Steady-state solve using estimated d_hat ---
        
        % Clamp the reference to a reachable target (bypassed based on code)
        % Dr_track = min(Dr(k), (y_max - margin) - y_ss);
        Dr_track = Dr(k);

        % Solve: (I-A)*Dx_bar = B*(Du_bar + d_hat) AND C*Dx_bar = Dr_track
        sol    = Mss \ [B * d_hat; Dr_track];
        Dx_bar = sol(1:n);
        Du_bar = sol(end);

        % Shifted control limits: du = Du - Du_bar, Du in [-u_ss, 100-u_ss]
        lb = (-u_ss        - Du_bar) * ones(H, 1);
        ub = ( 100 - u_ss  - Du_bar) * ones(H, 1);

        % Shifted output cap: y_hat(ii) <= y_max  ->  dy_hat(ii) <= y_max_inc
        y_max_inc = (y_max - y_ss) - Dr_track;

        % --- MPC: Solve in shifted coordinates using estimated state ---
        dx_k = Dx_hat - Dx_bar;

        if mode == 0
            [du_k, exitflag(k)] = mpc_solve(dx_k, H, R, A, B, C, ...
                                            lb, ub, y_max_inc, const_type);
        else
            [du_k, exitflag(k)] = mpc_solve_sparse_regularized(dx_k, H, R, ...
                                            A, B, C, lb, ub, y_max_inc, const_type);
        end

        Du(:, k) = du_k + Du_bar; % Incremental control
    else
        % Open-loop test of the Kalman filter: keep the equilibrium input
        Du(:, k) = 0;
    end

    % Compute the absolute control variable to apply
    u(:, k) = Du(:, k) + u_ss;

    % Saturate to the physical heater range before applying (safety)
    u(:, k) = min(max(u(:, k), 0), 100);

    % --- Kalman Filter PREDICTION Step ---
    % Propagate to the next sample using the control just applied
    xd_est(:, k+1) = Ad * xd_est(:, k) + Bd * Du(:, k);

    % Apply the control variable to the plant
    h1(u(1, k));

    % --- Real-Time Plotting ---
    if rt_plot
        clf;
        
        % Temperature subplot
        subplot(2, 1, 1); hold on; grid on;
        plot(t(1:k), y(1, 1:k), '.', 'MarkerSize', 10);
        stairs(t, r, 'g--');
        yline(y_max, 'r--');
        xlabel('Time [s]');
        ylabel('y [°C]');
        
        % Control subplot
        subplot(2, 1, 2); hold on; grid on;
        stairs(t(1:k), u(1, 1:k), 'LineWidth', 2);
        xlabel('Time [s]');
        ylabel('u [%]');
        ylim([0 100]);
        
        drawnow;
    else
        fprintf('t = %d, y1 = %.1f C, u1 = %.1f\n', t(k), y(1, k), u(1, k)); %#ok<UNRCH>
    end

    % Check for timing violations
    if toc > Ts
        warning('Computation time exceeded sampling time by %.2f s at sample %d.', toc - Ts, k);
    end
    
    % Wait for the remainder of the sampling interval
    pause(max(0, Ts - toc));
end

%% 8. Hardware Shutdown

% Turn off both heaters at the end of the experiment
h1(0);
h2(0);

% Signal the end of the experiment by shutting off the LED
led(0);

disp('Temperature test complete.');

%% 9. Final Diagnostic Plots

% Reconstruct the estimated output for plotting / discussion
y_hat     = y_ss + C * xd_est(1:n, 1:N);
d_hat_log = xd_est(end, 1:N);

figure('Units', 'normalized', 'Position', [0.15 0.3 0.4 0.55]);

% Plot 1: Output Tracking
subplot(3, 1, 1); hold on; grid on;
plot(t(1:k), y(1, 1:k), '.', 'MarkerSize', 10);
plot(t(1:k), y_hat(1:k), '-', 'LineWidth', 1.2);
stairs(t, r, 'g--');
yline(y_max, 'r--');
xlabel('Time [s]'); 
ylabel('y [°C]');
legend('measured y', 'estimated $\hat{y}$', 'reference r', '$y_{max}$', ...
       'Interpreter', 'latex', 'Location', 'best');
title('P5 — Closed-loop experiment on the real plant');

% Plot 2: Control Effort
subplot(3, 1, 2); hold on; grid on;
stairs(t(1:k), u(1, 1:k), 'LineWidth', 2);
yline(0, 'r--'); 
yline(100, 'r--'); 
yline(u_ss, 'k--');
xlabel('Time [s]'); 
ylabel('u [%]');
ylim([0 100]);

% Plot 3: Disturbance Estimation
subplot(3, 1, 3); hold on; grid on;
plot(t(1:k), d_hat_log(1:k), 'LineWidth', 1.5);
xlabel('Time [s]'); 
ylabel('$\hat{d}$ [%]', 'Interpreter', 'latex');
legend('estimated $\hat{d}$', 'Location', 'best');

%% 10. Save Results

% Save figure and experiment data to file
exportgraphics(gcf, ['closedloop_plot_', time_str, '.png'], 'Resolution', 300);
save('closedloop_data_1.mat', 'y', 'u', 't', 'r', 'y_hat', 'd_hat_log', 'exitflag');

%--------------------------------------------------------------------------