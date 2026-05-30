% Simulation of the TCLab linear model -- P4 Question 7
%
% Zero-offset MPC tracker with feedforward + Kalman filter (Figure 7).
% The MPC receives the ESTIMATED state Dx_hat (not the real Dx), and the
% feedforward terms Dx_bar, Du_bar are computed using the ESTIMATED input
% disturbance d_hat.
%
% Reference r commutes through 50, 40, 60 (unreachable), 45 degrees.
%
% Afonso Botelho and J. Miranda Lemos, IST  --  Q7 version
%__________________________________________________________________________
clear all; close all; clc

% ── Load model ────────────────────────────────────────────────────────────
load('singleheater_model2D.mat','A','B','C','Ke','e_var','y_ss','u_ss','Ts');
n = size(A,1);

% Noise: keep ON for Q7 (the guide requires it).
e_std = sqrt(e_var);
% e_std = 0;   % uncomment only for clean debugging

% ── MPC tuning (fixed from earlier questions) ─────────────────────────────
H = 50;          % prediction horizon
R = 0.04;        % control weight
mode       = 1;  % 0 = dense, 1 = sparse
const_type = 0;  % 0 = hard, 1 = soft  (soft recommended for the 60 deg leg)

% ── Safety constraint ─────────────────────────────────────────────────────
y_max  = 55;             % hard safety cap [degC]
margin = 0.5;            % clamp reference this far BELOW the cap [degC]

% ── Plant input disturbance (unknown to the MPC model) ────────────────────
% Applied to the plant as B*(u + d_true). The Kalman filter must estimate it.
d_true = 1.15;           % constant heater offset [%]

% ── Reference schedule: 50 -> 40 -> 60 -> 45 degC ─────────────────────────
r_levels  = [50 40 60 45];          % absolute reference [degC]
seg_time  = 400;                    % seconds held at each level
T  = seg_time*numel(r_levels);      % total duration [s]
N  = round(T/Ts);                   % number of samples
t  = (0:N-1)*Ts;

% r(k) as an absolute signal
r_abs = zeros(1,N);
for s = 1:numel(r_levels)
    idx = (1:N) > (s-1)*seg_time/Ts & (1:N) <= s*seg_time/Ts;
    r_abs(idx) = r_levels(s);
end

% ── Bias terms so the absolute simulator matches the incremental model ────
x_ss = [eye(n)-A; C] \ [B*u_ss; y_ss];
c1   = (eye(n)-A)*x_ss - B*u_ss;
c2   = y_ss - C*x_ss;

% Plant handles (real system): h1 advances the REAL state, T1C reads temp.
% d_true enters here -- the MPC model never sees it.
h1  = @(x,u) A*x + B*(u + d_true) + Ke*e_std*randn + c1;
T1C = @(x)   C*x + e_std*randn + c2;

% ═════════════════════════════════════════════════════════════════════════
%  Kalman filter setup -- augmented state xd = [Dx ; d]
% ═════════════════════════════════════════════════════════════════════════
Ad = [A,            B;
      zeros(1,n),   1];
Bd = [B; 0];
Cd = [C, 0];

QE  = Ke * e_var * Ke';   % process-noise covariance (identified)
RE  = e_var;              % measurement-noise covariance (identified)
deltaE = 1;               % disturbance-state covariance -- TUNING KNOB

QEd = [QE,          zeros(n,1);
       zeros(1,n),  deltaE];

L = dlqe(Ad, eye(n+1), Cd, QEd, RE);   % Kalman gain

% ── Signal storage ────────────────────────────────────────────────────────
x        = nan(n,   N+1);   % real plant state (absolute)
y        = nan(1,   N);
Du       = nan(1,   N);
u        = nan(1,   N);
xd_hat   = nan(n+1, N);     % augmented estimate [Dx_hat ; d_hat]
y_hat    = nan(1,   N);
exitflag = nan(1,   N);
solve_ms = nan(1,   N);     % MPC solve time per step [ms]

% ── Initial conditions ────────────────────────────────────────────────────
% Real plant starts at ambient equilibrium (u = 0 -> Dy = 0).
x(:,1) = x_ss; % + pinv(C)*(-20)

% Estimator initialised WITH AN ERROR equivalent to ~5 degC of output error,
% as the guide requests, and with zero initial disturbance estimate.
Dx_hat0   = pinv(C) * 5;          % C*Dx_hat0 = 5 degC
xd_hat(:,1) = [Dx_hat0; 0];

% ═════════════════════════════════════════════════════════════════════════
%  Closed-loop simulation
% ═════════════════════════════════════════════════════════════════════════
fprintf('Running Q7 simulation ...\n');
for k = 1:N

    % ── 1. SENSE: read the real plant ────────────────────────────────────
    y(k)  = T1C(x(:,k));
    Dy_k  = y(k) - y_ss;

    % ── 2. ESTIMATE: Kalman correction with the latest measurement ───────
    % Prediction was done at the end of the previous iteration (xd_pred).
    % On the first step we correct the initial guess directly.
    if k == 1
        xd_pred = xd_hat(:,1);
    end
    xd_hat(:,k) = xd_pred + L*(Dy_k - Cd*xd_pred);

    Dx_hat = xd_hat(1:n, k);     % estimated incremental state
    d_hat  = xd_hat(end,  k);    % estimated input disturbance
    y_hat(k) = y_ss + C*Dx_hat;

    % ── 3. FEEDFORWARD: steady-state solve using d_hat ───────────────────
    % Reference for this step (incremental), clamped to a reachable target.
    Dr_raw   = r_abs(k) - y_ss;
    %Dr_track = min(Dr_raw, (y_max - margin) - y_ss);   % <-- key fix
    Dr_track = Dr_raw;
    % Solve  (I-A)Dx_bar = B(Du_bar + d_hat),  C Dx_bar = Dr_track
    Mss = [eye(n)-A, -B;
           C,         0];
    bss = [B*d_hat; Dr_track];
    sol = Mss \ bss;
    Dx_bar = sol(1:n);
    Du_bar = sol(end);

    % Shifted control limits  du = Du - Du_bar,  Du in [-u_ss, 100-u_ss]
    lb = (-u_ss       - Du_bar) * ones(H,1);
    ub = ( 100 - u_ss - Du_bar) * ones(H,1);

    % Shifted output cap   y_hat(i) <= y_max  ->  dy_hat(i) <= y_max_inc
    y_max_inc = (y_max - y_ss) - Dr_track;

    % ── 4. MPC: solve in shifted coordinates using the ESTIMATED state ───
    dx_k = Dx_hat - Dx_bar;                 % dx = Dx_hat - Dx_bar

    tic;
    if mode == 0
        [du_k, exitflag(k)] = mpc_solve(dx_k, H, R, A, B, C, lb, ub, y_max_inc);
    else
        %[du_k, exitflag(k)] = mpc_solve_sparse_regularized( ...
        %                        dx_k, H, R, A, B, C, lb, ub, y_max_inc, const_type);
        [du_k, exitflag(k)] = mpc_solve_sparse_regular( ...
                                dx_k, H, R, A, B, C, lb, ub, y_max_inc, const_type);
    end
    solve_ms(k) = toc*1e3;

    Du(k) = du_k + Du_bar;       % incremental control
    u(k)  = u_ss + Du(k);        % absolute control commanded to the plant

    % ── 5. ACT: advance the real plant ───────────────────────────────────
    x(:,k+1) = h1(x(:,k), u(k));

    % ── 6. KALMAN PREDICTION for the next step (uses Du(k)) ──────────────
    xd_pred = Ad*xd_hat(:,k) + Bd*Du(k);
end
fprintf(' Done.\n');

% ═════════════════════════════════════════════════════════════════════════
%  Diagnostics
% ═════════════════════════════════════════════════════════════════════════
fprintf('Infeasible MPC steps : %d / %d\n', sum(exitflag ~= 1), N);
fprintf('MPC solve time       : mean %.2f ms, max %.2f ms\n', ...
        mean(solve_ms), max(solve_ms));
fprintf('10%% of Ts            : %.1f ms  (rule of thumb: solve < this)\n', ...
        0.10*Ts*1e3);
fprintf('Final d_hat          : %.4f  (true %.4f)\n', xd_hat(end,N), d_true);

% ═════════════════════════════════════════════════════════════════════════
%  Plots
% ═════════════════════════════════════════════════════════════════════════
figure('Units','normalized','Position',[0.15 0.3 0.4 0.55])

subplot(3,1,1); hold on; grid on
plot(t, y,     '.', 'MarkerSize', 6)
plot(t, y_hat, '-', 'LineWidth', 1.2)
stairs(t, r_abs, 'g--', 'LineWidth', 1.5)
yline(y_max, 'r--', 'LineWidth', 1.2)
xlabel('Time [s]'); ylabel('y [°C]')
legend('measured y','estimated $\hat{y}$','reference r','$y_{max}$', ...
       'Interpreter','latex','Location','best')
title('Q7 — Output tracking with Kalman filter')

subplot(3,1,2); hold on; grid on
stairs(t, u, 'LineWidth', 1.5)
yline(0,  'r--'); yline(100,'r--'); yline(u_ss,'k--')
xlabel('Time [s]'); ylabel('u [%]')
legend('u','','','$\bar{u}$','Interpreter','latex','Location','best')

subplot(3,1,3); hold on; grid on
plot(t, xd_hat(end,:), 'LineWidth', 1.5)
yline(d_true, 'k--', 'LineWidth', 1.2)
xlabel('Time [s]'); ylabel('$\hat{d}$ [%]','Interpreter','latex')
legend('estimated $\hat{d}$','true d','Interpreter','latex','Location','best')

%--------------------------------------------------------------------------
% End of File