%% Kalman Filter for Disturbance Estimation (TCLab)
%--------------------------------------------------------------------------
% Simulates the real plant with a constant input disturbance and 
% measurement noise, and uses an augmented Kalman filter to jointly 
% estimate the system states and the unknown disturbance.
%
% Expected variables in the .mat file:
%   A, B, C : State-space matrices
%   Ke      : Kalman gain from system identification
%   e_var   : Measurement noise variance
%   u_ss    : Steady-state input
%   y_ss    : Steady-state output
%   Ts      : Sampling time
%
% Matlab Toolboxes: Control System Toolbox
% Functions: dlqe
%--------------------------------------------------------------------------

%% 1. Initialization and Setup

clear; clc; close all;

% Load identified model parameters
load('singleheater_model_n2.mat');

n = size(A, 1);        % State vector dimension
N = 800;               % Number of simulation steps
t = (0:N-1) * Ts;      % Time vector [s]

% Open-loop input sequence (constant steady-state)
u  = u_ss * ones(1, N);
Du = u - u_ss;         % Incremental input (0 in this case)

% Real plant simulation parameters
d_true = 1.15;         % True input disturbance (bias)
e_std  = sqrt(e_var);  % Measurement noise standard deviation

%% 2. Augmented Model Definition
% Augment the state vector to include the disturbance model: x_d = [x; d]
% Assumes a constant disturbance model: d(k+1) = d(k) + noise

Ad = [A, B;
      zeros(1, n), 1];

Bd = [B;
      0];

Cd = [C, 0];

%% 3. Kalman Filter Design

% Base covariances from system identification
Q_E = Ke * e_var * Ke';
R_E = e_var;

% Tuning parameter for the disturbance estimation speed
delta_e = 1; 

% Augmented process noise covariance matrix
Q_Ed = [Q_E, zeros(n, 1);
        zeros(1, n), delta_e];

% Compute steady-state Kalman gain for the augmented system
L = dlqe(Ad, eye(n+1), Cd, Q_Ed, R_E);

%% 4. Pre-allocation and Initial Conditions

% Real plant state variables
Dx = zeros(n, N); 
Dy = zeros(1, N); 
y  = zeros(1, N); 

% Estimated variables
xd_hat = zeros(n+1, N); % Augmented state estimate
Dx_hat = zeros(n, N);   % State estimate
d_hat  = zeros(1, N);   % Disturbance estimate
Dy_hat = zeros(1, N);   % Incremental output estimate
y_hat  = zeros(1, N);   % Total output estimate

% Apply an initial estimation error (approx. 5 degrees Celsius)
xd_hat(1:n, 1) = pinv(C) * 5;

% Initialize estimated arrays for k = 1
Dx_hat(:, 1) = xd_hat(1:n, 1);
d_hat(1)     = xd_hat(end, 1);
Dy_hat(1)    = C * Dx_hat(:, 1);
y_hat(1)     = y_ss + Dy_hat(1);

% Initialize real plant arrays for k = 1 with noise
e_initial = e_std * randn;
Dy(1)     = C * Dx(:, 1) + e_initial;
y(1)      = y_ss + Dy(1);

%% 5. Simulation Loop

for k = 1:N-1

    % Generate measurement noise for the current step
    e = e_std * randn;

    % 1. Simulate real plant step
    Dx(:, k+1) = A * Dx(:, k) + B * (Du(k) + d_true) + Ke * e;
    Dy(k+1)    = C * Dx(:, k+1) + e;
    y(k+1)     = y_ss + Dy(k+1);

    % 2. Kalman prediction step (using previous estimate and known input)
    xd_pred = Ad * xd_hat(:, k) + Bd * Du(k);

    % 3. Kalman correction step (using new measurement)
    xd_hat(:, k+1) = xd_pred + L * (Dy(k+1) - Cd * xd_pred);

    % 4. Extract and store estimates
    Dx_hat(:, k+1) = xd_hat(1:n, k+1);
    d_hat(k+1)     = xd_hat(end, k+1);

    Dy_hat(k+1) = C * Dx_hat(:, k+1);
    y_hat(k+1)  = y_ss + Dy_hat(k+1);

end

%% 6. Plot Results

figure('Units', 'normalized', 'Position', [0.2 0.1 0.6 0.75]);

% Output estimation plot
subplot(2, 1, 1); hold on; grid on;
plot(t, y, 'LineWidth', 1.5, 'DisplayName', 'Measured (y)');
plot(t, y_hat, '--', 'LineWidth', 1.5, 'DisplayName', 'Estimated (\hat{y})');
xlabel('Time [s]');
ylabel('Temperature [°C]');
legend('Location', 'best');
title('Kalman Filter Output Estimation');

% Disturbance estimation plot
subplot(2, 1, 2); hold on; grid on;
plot(t, d_hat, 'LineWidth', 1.5, 'DisplayName', 'Estimated (\hat{d})');
yline(d_true, '--k', 'LineWidth', 1.5, 'DisplayName', 'True (d)');
xlabel('Time [s]');
ylabel('Disturbance');
legend('Location', 'best');
title('Estimated Input Disturbance');

sgtitle(sprintf('Augmented Kalman Filter Performance (\\delta_E = %g)', delta_e));