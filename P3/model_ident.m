%% System Identification and Validation (TCLab)
%  -----------------------------------------------------------------------
%  Identifies state-space models of varying orders from open-loop data
%  and evaluates their performance against a separate validation dataset.
%  -----------------------------------------------------------------------

clear; clc; close all;

%% 1. General Settings
n_save = 3;  % Model order to save to disk
n_min  = 1;  % Minimum model order to evaluate
n_max  = 5;  % Maximum model order to evaluate

%% 2. Load and Prepare Identification Data
load('openloop_data_1.mat', 'y', 'u', 't');

% Extract relevant vectors
u1_raw = u(1, :);
y1_raw = y(1, :);
t1_raw = t;

% Compute steady-state equilibrium values
k_ss_begin = 201;
k_ss_end   = 400;
y_ss = mean(y1_raw(:, k_ss_begin:k_ss_end), 2);
u_ss = u1_raw(:, k_ss_begin);

% Truncate initial transients
t1 = t1_raw(k_ss_begin:end-1);
u1 = u1_raw(:, k_ss_begin:end-1);
y1 = y1_raw(:, k_ss_begin:end-1);
Ts = t1(2) - t1(1); % Sampling time

% Compute incremental variables (deviations from steady-state)
Dy1 = y1 - y_ss;
Du1 = u1 - u_ss;
N1  = length(t1);

%% 3. Load and Prepare Validation Data
load('openloop_data_2.mat', 'y', 'u', 't');

u2_raw = u(1, :);
y2_raw = y(1, :);
t2     = t;

% Compute incremental variables for the validation set
Dy2 = y2_raw - y_ss;
Du2 = u2_raw - u_ss;
N2  = length(t2);

%% 4. Pre-allocate Storage for Metrics
n_orders    = n_max - n_min + 1;
color_map   = lines(n_orders);

mse_all     = zeros(1, n_orders); % Mean Squared Error
rmse_db_all = zeros(1, n_orders); % De-biased Root Mean Squared Error
bias_all    = zeros(1, n_orders); % Mean error (offset)
Dy2_sims    = cell(1, n_orders);  % Simulated outputs for validation

%% 5. Model Identification and Evaluation Loop
for n = n_min:n_max
    
    idx = n - n_min + 1;
    
    % --- System Identification ---
    sys = ssest(Du1', Dy1', n, 'Ts', Ts);
    [A, B, C, ~, Ke] = idssdata(sys);
    e_var = sys.NoiseVariance;
    
    % Save the designated model order
    if n == n_save
        save('singleheater_model.mat', 'A', 'B', 'C', 'Ke', 'e_var', 'y_ss', 'u_ss', 'Ts');
    end
    
    % --- Simulate on Identification Dataset ---
    Dy1_sim = nan(1, N1);
    Dx1_sim = nan(n, N1);
    
    % Estimate initial states
    Dx0 = findstates(sys, iddata(Dy1', Du1', Ts));
    Dx1_sim(:, 1) = Dx0;
    Dy1_sim(:, 1) = Dy1(:, 1);
    
    for k = 1:N1-1
        Dx1_sim(:, k+1) = A * Dx1_sim(:, k) + B * Du1(:, k);
        Dy1_sim(:, k+1) = C * Dx1_sim(:, k+1);
    end
    
    % --- Simulate on Validation Dataset ---
    Dy2_sim = nan(1, N2);
    Dx2_sim = nan(n, N2);
    
    Dx02 = findstates(sys, iddata(Dy2', Du2', Ts));
    Dx2_sim(:, 1) = Dx02;
    Dy2_sim(:, 1) = Dy2(:, 1);
    
    for k = 1:N2-1
        Dx2_sim(:, k+1) = A * Dx2_sim(:, k) + B * Du2(:, k);
        Dy2_sim(:, k+1) = C * Dx2_sim(:, k+1);
    end
    
    % --- Compute Validation Metrics ---
    err    = Dy2_sim - Dy2; % Signed error
    bias   = mean(err);     % Vertical offset
    err_db = err - bias;    % De-biased error (isolates shape mismatch)
    
    mse     = mean(err.^2);
    rmse_db = sqrt(mean(err_db.^2));
    
    mse_all(idx)     = mse;
    rmse_db_all(idx) = rmse_db;
    bias_all(idx)    = bias;
    Dy2_sims{idx}    = Dy2_sim;
    
    % --- Plot Individual Model Performance ---
    figure('Units', 'normalized', 'Position', [0.05 0.1 0.88 0.75], ...
           'Name', sprintf('Order n=%d', n));
       
    % Identification Output
    subplot(2, 2, 1); hold on; grid on;
    title(sprintf('n=%d | Identification dataset – output', n));
    plot(t1, Dy1, '.', 'MarkerSize', 4, 'Color', [0.3 0.5 1]);
    plot(t1, Dy1_sim, '-', 'Color', color_map(idx, :), 'LineWidth', 1.5);
    xlabel('Time [s]'); ylabel('\Delta y [°C]');
    xlim([t1(1), t1(end)]);
    legend('Experimental', 'Model', 'Location', 'best');
    
    % Identification Input
    subplot(2, 2, 3); hold on; grid on;
    title(sprintf('n=%d | Identification dataset – input', n));
    stairs(t1, Du1, 'LineWidth', 2, 'Color', [0.2 0.2 0.2]);
    xlabel('Time [s]'); ylabel('\Delta u [%]');
    xlim([t1(1), t1(end)]);
    
    % Validation Output
    subplot(2, 2, 2); hold on; grid on;
    title(sprintf('n=%d | Validation – output (MSE=%.4f, RMSE_{db}=%.4f, bias=%.4f)', ...
                  n, mse, rmse_db, bias));
    plot(t2, Dy2, '.', 'MarkerSize', 4, 'Color', [0.3 0.5 1]);
    plot(t2, Dy2_sim, '-', 'Color', color_map(idx, :), 'LineWidth', 1.5);
    xlabel('Time [s]'); ylabel('\Delta y [°C]');
    xlim([t2(1), t2(end)]);
    legend('Experimental', 'Model', 'Location', 'best');
    
    % Validation Error
    subplot(2, 2, 4); hold on; grid on;
    title(sprintf('n=%d | Validation dataset – error (raw and de-biased)', n));
    plot(t2, err, '-', 'Color', color_map(idx, :), 'LineWidth', 1);
    plot(t2, err_db, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    yline(bias, ':k', sprintf('bias = %.3f', bias), 'LineWidth', 1);
    xlabel('Time [s]'); ylabel('\Delta y error [°C]');
    xlim([t2(1), t2(end)]);
    legend('Raw error', 'De-biased error', 'Location', 'best');
    
    % Console Output
    fprintf('n=%d | MSE = %.4f | RMSE_db = %.4f | bias = %+.4f\n', ...
            n, mse, rmse_db, bias);
end

%% 6. Comparative Validation Plot (All Orders)
figure('Units', 'normalized', 'Position', [0.15 0.15 0.65 0.5], ...
       'Name', 'Comparative validation – all orders');
hold on; grid on;
title('Validation output: comparison across model orders');

plot(t2, Dy2, '.', 'MarkerSize', 10, 'Color', [0 0 0], ...
     'DisplayName', 'Experimental data');

for n = n_min:n_max
    idx = n - n_min + 1;
    plot(t2, Dy2_sims{idx}, '-', 'Color', color_map(idx, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('n=%d (MSE=%.4f, RMSE_{db}=%.4f)', ...
                                n, mse_all(idx), rmse_db_all(idx)));
end

xlabel('Time [s]');
ylabel('\Delta y [°C]');
xlim([t2(1), t2(end)]);
legend('Location', 'best');

%% 7. Summary Bar Plots for Metrics
order_labels = arrayfun(@(k) sprintf('n=%d', k), n_min:n_max, 'UniformOutput', false);

% Define bar colors (highlighting the saved model order)
bar_colors_mse = repmat([0.30 0.55 0.85], n_orders, 1);
bar_colors_rms = repmat([0.85 0.45 0.30], n_orders, 1);

sel_idx = n_save - n_min + 1;
bar_colors_mse(sel_idx, :) = [0.10 0.70 0.30];
bar_colors_rms(sel_idx, :) = [0.10 0.70 0.30];

figure('Units', 'normalized', 'Position', [0.20 0.20 0.65 0.45], ...
       'Name', 'Validation metrics – bar plots');

% MSE Bar Plot
subplot(1, 2, 1);
b1 = bar(mse_all, 'FaceColor', 'flat'); 
b1.CData = bar_colors_mse;
grid on;
set(gca, 'XTickLabel', order_labels);
ylabel('MSE [°C^{2}]');
title('Metric 1: Mean Squared Error');
ylim([0, max(mse_all) * 1.18]);

for kk = 1:n_orders
    text(kk, mse_all(kk), sprintf('%.3f', mse_all(kk)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 10);
end

% RMSE Bar Plot
subplot(1, 2, 2);
b2 = bar(rmse_db_all, 'FaceColor', 'flat'); 
b2.CData = bar_colors_rms;
grid on;
set(gca, 'XTickLabel', order_labels);
ylabel('RMSE_{db} [°C]');
title('Metric 2: Bias-removed RMSE (Shape Only)');
ylim([0, max(rmse_db_all) * 1.18]);

for kk = 1:n_orders
    text(kk, rmse_db_all(kk), sprintf('%.3f', rmse_db_all(kk)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 10);
end

sgtitle('Model performance on validation dataset – selected order in green');