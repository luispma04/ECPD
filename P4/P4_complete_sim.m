%% Complete MPC Tracker with Feedforward and Kalman Filter
%--------------------------------------------------------------------------
% Simulates the TCLab linear model (P4 Question 7). 
%
% The MPC controller receives the ESTIMATED state (Dx_hat) from an augmented 
% Kalman filter rather than the real state. The feedforward references 
% (Dx_bar, Du_bar) are dynamically computed using the ESTIMATED input 
% disturbance (d_hat) to achieve zero-offset tracking.
% 
% The reference signal steps through: 50 -> 40 -> 60 (unreachable) -> 45 °C.
%
% Required variables in .mat: A, B, C, Ke, e_var, y_ss, u_ss, Ts
%
% Matlab Toolboxes: Control System Toolbox, Optimization Toolbox
% Functions: dlqe, mpc_solve, mpc_solve_sparse_regular
%--------------------------------------------------------------------------

%% 1. Initialization and Setup

clear; clc; close all;

% Load identified model parameters
load('singleheater_model2D.mat', 'A', 'B', 'C', 'Ke', 'e_var', 'y_ss', 'u_ss', 'Ts');
n = size(A, 1);

% Measurement noise setup (Keep ON for Q7 as per guide)
e_std = sqrt(e_var);
% e_std = 0; % Uncomment only for clean debugging

%% 2. MPC Tuning Parameters

H          = 50;   % Prediction horizon
R          = 0.04; % Control weight
mode       = 1;    % Solver mode: 0 = dense, 1 = sparse
const_type = 0;    % Constraint type: 0 = hard, 1 = soft (soft rec. for 60°C leg)

% Safety constraints
y_max  = 55;       % Hard safety cap [°C]
margin = 0.5;      % Margin below the cap [°C] (Currently disabled in tracking logic)

% Plant input disturbance (Unknown to the MPC model)
d_true = 1.15;     % Constant heater offset applied to the real plant [%]

%% 3. Reference Schedule Definition

r_levels = [50 40 60 45];             % Absolute reference targets [°C]
seg_time = 400;                       % Time held at each level [s]
T_total  = seg_time * numel(r_levels);% Total simulation duration [s]
N        = round(T_total / Ts);       % Total number of samples
t        = (0:N-1) * Ts;              % Time vector [s]

% Build absolute reference signal
r_abs = zeros(1, N);
for s = 1:numel(r_levels)
    idx = (1:N) > (s-1)*seg_time/Ts & (1:N) <= s*seg_time/Ts;
    r_abs(idx) = r_levels(s);
end

%% 4. Plant Bias Terms and Real System Handles

% Compute state and output biases to match absolute simulator with incremental model
x_ss = [eye(n) - A; C] \ [B * u_ss; y_ss];
c1   = (eye(n) - A) * x_ss - B * u_ss;
c2   = y_ss - C * x_ss;

% Plant functions (h1 advances REAL state, T1C reads REAL temp)
% Note: d_true enters here. The MPC model never sees it directly.
h1  = @(x, u) A * x + B * (u + d_true) + Ke * e_std * randn + c1;
T1C = @(x)    C * x + e_std * randn + c2;

%% 5. Augmented Kalman Filter Setup (xd = [Dx; d])

Ad = [A,          B;
      zeros(1,n), 1];
Bd = [B; 0];
Cd = [C, 0];

Q_E     = Ke * e_var * Ke'; % Process-noise covariance (identified)
R_E     = e_var;            % Measurement-noise covariance (identified)
delta_e = 1;                % Disturbance-state covariance (Tuning knob)

Q_Ed = [Q_E,          zeros(n, 1);
        zeros(1, n),  delta_e];

% Compute steady-state Kalman gain
L = dlqe(Ad, eye(n+1), Cd, Q_Ed, R_E);

%% 6. Pre-allocation and Initial Conditions

% Signal storage arrays
x        = nan(n,   N+1); % Real plant state (absolute)
y        = nan(1,   N);   % Measured output
Du       = nan(1,   N);   % Incremental control effort
u        = nan(1,   N);   % Absolute control effort
xd_hat   = nan(n+1, N);   % Augmented estimate [Dx_hat; d_hat]
y_hat    = nan(1,   N);   % Estimated output
exitflag = nan(1,   N);   % Solver status flag
solve_ms = nan(1,   N);   % MPC solver execution time [ms]

% Real plant starts at ambient equilibrium
x(:, 1) = x_ss; 

% Estimator initialized with an ~5°C output error (as requested)
Dx_hat0        = pinv(C) * 5; 
xd_hat(:, 1)   = [Dx_hat0; 0]; % Initial disturbance estimate is zero

%% 7. Closed-Loop Simulation

fprintf('Running Q7 simulation ...\n');

for k = 1:N

    % --- Step 1: SENSE (Read real plant) ---
    y(k) = T1C(x(:, k));
    Dy_k = y(k) - y_ss;

    % --- Step 2: ESTIMATE (Kalman Correction) ---
    if k == 1
        xd_pred = xd_hat(:, 1); % First step corrects initial guess directly
    end
    
    xd_hat(:, k) = xd_pred + L * (Dy_k - Cd * xd_pred);

    Dx_hat   = xd_hat(1:n, k);  % Estimated incremental state
    d_hat    = xd_hat(end, k);  % Estimated input disturbance
    y_hat(k) = y_ss + C * Dx_hat;

    % --- Step 3: FEEDFORWARD (Steady-state solve using d_hat) ---
    Dr_raw   = r_abs(k) - y_ss;
    Dr_track = Dr_raw; % (margin limit logic preserved but bypassed as requested)
    
    % Solve: (I-A)*Dx_bar = B*(Du_bar + d_hat) AND C*Dx_bar = Dr_track
    Mss = [eye(n) - A, -B;
           C,           0];
    bss = [B * d_hat; Dr_track];
    
    sol    = Mss \ bss;
    Dx_bar = sol(1:n);
    Du_bar = sol(end);

    % Shifted bounds for standard coordinates: Du in [-u_ss, 100-u_ss]
    lb = (-u_ss        - Du_bar) * ones(H, 1);
    ub = ( 100 - u_ss  - Du_bar) * ones(H, 1);

    % Shifted output cap: y_hat(ii) <= y_max  -->  dy_hat(ii) <= y_max_inc
    y_max_inc = (y_max - y_ss) - Dr_track;

    % --- Step 4: COMPUTE MPC (Solve in shifted coordinates) ---
    dx_k = Dx_hat - Dx_bar; 

    tic;
    if mode == 0
        [du_k, exitflag(k)] = mpc_solve(dx_k, H, R, A, B, C, lb, ub, y_max_inc);
    else
        % [du_k, exitflag(k)] = mpc_solve_sparse_regularized( ...
        %                        dx_k, H, R, A, B, C, lb, ub, y_max_inc, const_type);
        [du_k, exitflag(k)] = mpc_solve_sparse_regular( ...
                                dx_k, H, R, A, B, C, lb, ub, y_max_inc, const_type);
    end
    solve_ms(k) = toc * 1e3;

    % Recover absolute and incremental control signals
    Du(k) = du_k + Du_bar; 
    u(k)  = u_ss + Du(k);  

    % --- Step 5: ACT (Advance real plant) ---
    x(:, k+1) = h1(x(:, k), u(k));

    % --- Step 6: PREDICT (Kalman Prediction for next step) ---
    xd_pred = Ad * xd_hat(:, k) + Bd * Du(k);
    
end

fprintf(' Done.\n');

%% 8. Diagnostics

fprintf('\nInfeasible MPC steps : %d / %d\n', sum(exitflag ~= 1), N);
fprintf('MPC solve time       : mean %.2f ms, max %.2f ms\n', mean(solve_ms), max(solve_ms));
fprintf('10%% of Ts            : %.1f ms (rule of thumb: solve < this)\n', 0.10 * Ts * 1e3);
fprintf('Final d_hat          : %.4f (true %.4f)\n', xd_hat(end, N), d_true);

%% 9. Plot Results

figure('Units', 'normalized', 'Position', [0.15 0.3 0.4 0.55]);

% Plot 1: Output Tracking
subplot(3, 1, 1); hold on; grid on;
plot(t, y, '.', 'MarkerSize', 6);
plot(t, y_hat, '-', 'LineWidth', 1.2);
stairs(t, r_abs, 'g--', 'LineWidth', 1.5);
yline(y_max, 'r--', 'LineWidth', 1.2);
xlabel('Time [s]'); 
ylabel('y [°C]');
legend('measured y', 'estimated $\hat{y}$', 'reference r', '$y_{max}$', ...
       'Interpreter', 'latex', 'Location', 'best');
title('Q7 — Output tracking with Kalman filter');

% Plot 2: Control Effort
subplot(3, 1, 2); hold on; grid on;
stairs(t, u, 'LineWidth', 1.5);
yline(0, 'r--'); 
yline(100, 'r--'); 
yline(u_ss, 'k--');
xlabel('Time [s]'); 
ylabel('u [%]');
legend('u', '', '', '$\bar{u}$', 'Interpreter', 'latex', 'Location', 'best');

% Plot 3: Disturbance Estimation
subplot(3, 1, 3); hold on; grid on;
plot(t, xd_hat(end, :), 'LineWidth', 1.5);
yline(d_true, 'k--', 'LineWidth', 1.2);
xlabel('Time [s]'); 
ylabel('$\hat{d}$ [%]', 'Interpreter', 'latex');
legend('estimated $\hat{d}$', 'true d', 'Interpreter', 'latex', 'Location', 'best');

%--------------------------------------------------------------------------