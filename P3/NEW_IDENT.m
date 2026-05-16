clear
close all

%% Settings
n_save  = 3;
n_min   = 1;
n_max   = 5;

%% Load identification data (use distinct variable names)
load('openloop_data_1.mat','y','u','t');
u1_raw = u(1,:);
y1_raw = y(1,:);
t1_raw = t;

% Steady-state equilibrium
k_ss_begin = 201;
k_ss_end   = 400;
y_ss = mean(y1_raw(:,k_ss_begin:k_ss_end),2);
u_ss = u1_raw(:,k_ss_begin);

% Truncate initial transient
t1  = t1_raw(k_ss_begin:end-1);
u1  = u1_raw(:,k_ss_begin:end-1);
y1  = y1_raw(:,k_ss_begin:end-1);
Ts  = t1(2) - t1(1);

% Incremental variables for identification
Dy1 = y1 - y_ss;
Du1 = u1 - u_ss;
N1  = length(t1);

%% Load validation data (completely separate variables)
load('openloop_data_2.mat','y','u','t');
u2_raw = u(1,:);
y2_raw = y(1,:);
t2     = t;

% Incremental variables for validation
Dy2 = y2_raw - y_ss;
Du2 = u2_raw - u_ss;
N2  = length(t2);

%% Pre-allocate storage
n_orders     = n_max - n_min + 1;
colors       = lines(n_orders);
mse_all      = zeros(1, n_orders);   % Metric 1: Mean Squared Error
rmse_db_all  = zeros(1, n_orders);   % Metric 2: bias-removed (de-biased) RMSE
bias_all     = zeros(1, n_orders);   % auxiliary: mean error (offset)
Dy2_sims     = cell(1,  n_orders);

%% Loop over model orders
for n = n_min:n_max

    idx = n - n_min + 1;

    % ── Identify model ───────────────────────────────────────────────────
    sys = ssest(Du1', Dy1', n, 'Ts', Ts);
    [A,B,C,~,Ke] = idssdata(sys);
    e_var = sys.NoiseVariance;

    if n == n_save
        save('singleheater_model.mat','A','B','C','Ke','e_var','y_ss','u_ss','Ts');
    end

    % ── Simulate on identification dataset ───────────────────────────────
    Dy1_sim = nan(1,N1);
    Dx1_sim = nan(n,N1);
    Dx0     = findstates(sys, iddata(Dy1',Du1',Ts));
    Dx1_sim(:,1) = Dx0;
    Dy1_sim(:,1) = Dy1(:,1);
    for k = 1:N1-1
        Dx1_sim(:,k+1) = A*Dx1_sim(:,k) + B*Du1(:,k);
        Dy1_sim(:,k+1) = C*Dx1_sim(:,k+1);
    end

    % ── Simulate on validation dataset ───────────────────────────────────
    Dy2_sim = nan(1,N2);
    Dx2_sim = nan(n,N2);
    Dx02    = findstates(sys, iddata(Dy2',Du2',Ts));
    Dx2_sim(:,1) = Dx02;
    Dy2_sim(:,1) = Dy2(:,1);
    for k = 1:N2-1
        Dx2_sim(:,k+1) = A*Dx2_sim(:,k) + B*Du2(:,k);
        Dy2_sim(:,k+1) = C*Dx2_sim(:,k+1);
    end

    % ── Metrics on validation dataset ────────────────────────────────────
    err      = Dy2_sim - Dy2;          % signed error
    bias     = mean(err);              % vertical offset between curves
    err_db   = err - bias;             % de-biased error (shape mismatch only)

    mse      = mean(err.^2);           % Metric 1: MSE  (penalizes both shape and offset)
    rmse_db  = sqrt(mean(err_db.^2));  % Metric 2: bias-removed RMSE (shape only)

    mse_all(idx)     = mse;
    rmse_db_all(idx) = rmse_db;
    bias_all(idx)    = bias;
    Dy2_sims{idx}    = Dy2_sim;

    % ── Individual figure ─────────────────────────────────────────────────
    figure('Units','normalized','Position',[0.05 0.1 0.88 0.75], ...
           'Name', sprintf('Order n=%d',n));

    subplot(2,2,1), hold on, grid on
    title(sprintf('n=%d | Identification dataset – output',n))
    plot(t1, Dy1, '.', 'MarkerSize',4, 'Color',[0.3 0.5 1])
    plot(t1, Dy1_sim, '-', 'Color', colors(idx,:), 'LineWidth',1.5)
    xlabel('Time [s]'), ylabel('\Delta y [°C]')
    xlim([t1(1),t1(end)])
    legend('Experimental','Model','Location','best')

    subplot(2,2,3), hold on, grid on
    title(sprintf('n=%d | Identification dataset – input',n))
    stairs(t1, Du1, 'LineWidth',2, 'Color',[0.2 0.2 0.2])
    xlabel('Time [s]'), ylabel('\Delta u [%]')
    xlim([t1(1),t1(end)])

    subplot(2,2,2), hold on, grid on
    title(sprintf(['n=%d | Validation – output  ' ...
                   '(MSE=%.4f, RMSE_{db}=%.4f, bias=%.4f)'], ...
                   n, mse, rmse_db, bias))
    plot(t2, Dy2, '.', 'MarkerSize',4, 'Color',[0.3 0.5 1])
    plot(t2, Dy2_sim, '-', 'Color', colors(idx,:), 'LineWidth',1.5)
    xlabel('Time [s]'), ylabel('\Delta y [°C]')
    xlim([t2(1),t2(end)])
    legend('Experimental','Model','Location','best')

    subplot(2,2,4), hold on, grid on
    title(sprintf('n=%d | Validation dataset – error (raw and de-biased)',n))
    plot(t2, err,    '-',  'Color', colors(idx,:), 'LineWidth',1)
    plot(t2, err_db, '--', 'Color', [0.4 0.4 0.4], 'LineWidth',1)
    yline(bias, ':k', sprintf('bias = %.3f',bias), 'LineWidth',1)
    xlabel('Time [s]'), ylabel('\Delta y error [°C]')
    xlim([t2(1),t2(end)])
    legend('Raw error','De-biased error','Location','best')

    fprintf('n=%d | MSE = %.4f | RMSE_db = %.4f | bias = %+.4f\n', ...
            n, mse, rmse_db, bias);
end

%% Comparative validation plot
figure('Units','normalized','Position',[0.15 0.15 0.65 0.5], ...
       'Name','Comparative validation – all orders');
hold on, grid on
title('Validation output: comparison across model orders')

plot(t2, Dy2, '.', 'MarkerSize',10, 'Color',[0 0 0], ...
     'DisplayName','Experimental data')

for n = n_min:n_max
    idx = n - n_min + 1;
    plot(t2, Dy2_sims{idx}, '-', 'Color', colors(idx,:), 'LineWidth',1.5, ...
         'DisplayName', sprintf('n=%d  (MSE=%.4f, RMSE_{db}=%.4f)', ...
                                 n, mse_all(idx), rmse_db_all(idx)))
end

xlabel('Time [s]')
ylabel('\Delta y [°C]')
xlim([t2(1),t2(end)])
legend('Location','best')

%% Bar plots — one per metric
order_labels = arrayfun(@(k) sprintf('n=%d',k), n_min:n_max, ...
                        'UniformOutput', false);

% Highlight the chosen order in a different colour
bar_colors_mse = repmat([0.30 0.55 0.85], n_orders, 1);
bar_colors_rms = repmat([0.85 0.45 0.30], n_orders, 1);
sel = n_save - n_min + 1;
bar_colors_mse(sel,:) = [0.10 0.70 0.30];
bar_colors_rms(sel,:) = [0.10 0.70 0.30];

figure('Units','normalized','Position',[0.20 0.20 0.65 0.45], ...
       'Name','Validation metrics – bar plots');

% ── MSE ──
subplot(1,2,1)
b1 = bar(mse_all, 'FaceColor','flat'); b1.CData = bar_colors_mse;
grid on
set(gca,'XTickLabel',order_labels)
ylabel('MSE  [°C^{2}]')
title('Metric 1: Mean Squared Error (validation)')
% annotate numerical values on top of each bar
for k = 1:n_orders
    text(k, mse_all(k), sprintf('%.3f', mse_all(k)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',10)
end
ylim([0, max(mse_all)*1.18])

% ── bias-removed RMSE ──
subplot(1,2,2)
b2 = bar(rmse_db_all, 'FaceColor','flat'); b2.CData = bar_colors_rms;
grid on
set(gca,'XTickLabel',order_labels)
ylabel('RMSE_{db}  [°C]')
title('Metric 2: Bias-removed RMSE (validation, shape only)')
for k = 1:n_orders
    text(k, rmse_db_all(k), sprintf('%.3f', rmse_db_all(k)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',10)
end
ylim([0, max(rmse_db_all)*1.18])

sgtitle('Model performance on validation dataset – selected order in green')