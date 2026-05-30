%% Receding-Horizon Control Basics
%  -----------------------------------------------------------------------
% Compares the Receding-Horizon (RH) state-feedback gain with the
% infinite-horizon LQ gain for two first-order plants:
%   Plant 1 (unstable): x(t+1) = 1.2 x(t) + u(t)
%   Plant 2 (stable):   x(t+1) = 0.8 x(t) + u(t)
%
% Matlab Toolboxes: Control System Toolbox
% Functions: compute_KRH, dlqr
%--------------------------------------------------------------------------

clear; clc; close all;

%% 1. User-Defined Parameters
% Maximum horizon limit
h_max = 50; 

% Representative control weights (R)
r_values = [0.01, 0.1, 1, 10]; 

% State weight (Q), fixed
q_cost = 1; 

%% 2. System Definitions
% Define unstable plant
sys(1).name = 'Unstable plant (A = 1.2)';
sys(1).a    = 1.2;
sys(1).b    = 1;
sys(1).c    = 1;

% Define stable plant
sys(2).name = 'Stable plant (A = 0.8)';
sys(2).a    = 0.8;
sys(2).b    = 1;
sys(2).c    = 1;

%% 3. Main Loop Over Plants
for pp = 1:numel(sys)
    
    % Extract parameters for current plant
    sys_a = sys(pp).a;
    sys_b = sys(pp).b;
    sys_c = sys(pp).c;
    
    fprintf('\n========================================\n');
    fprintf('  %s\n', sys(pp).name);
    fprintf('========================================\n');
    
    %% Q1: Compute LQ Gain
    fprintf('\n--- Q1: LQ gain for each R (Q = %g) ---\n', q_cost);
    
    q_lq    = sys_c' * q_cost * sys_c;
    k_lq_vec = zeros(size(r_values));
    
    for rr = 1:numel(r_values)
        r_val = r_values(rr);
        [k_lq_temp, ~, ~] = dlqr(sys_a, sys_b, q_lq, r_val);
        k_lq_vec(rr) = k_lq_temp;
        fprintf('  R = %5.2f  ->  K_LQ = %8.5f,  CL eigenvalue = %+.5f\n', ...
            r_val, k_lq_temp, sys_a - sys_b * k_lq_temp);
    end
    
    %% Q2 & Q3: Compute RH Gains and Eigenvalues
    h_range    = 1:h_max;
    color_map  = lines(numel(r_values));
    
    k_rh_data  = zeros(numel(r_values), h_max);
    eig_data   = zeros(numel(r_values), h_max);
    
    for rr = 1:numel(r_values)
        r_val = r_values(rr);
        for hh = h_range
            k_rh_temp         = compute_KRH(sys_a, sys_b, sys_c, hh, r_val);
            k_rh_data(rr, hh) = k_rh_temp;
            eig_data(rr, hh)  = abs(sys_a - sys_b * k_rh_temp);
        end
    end
    
    %% Plotting: Feedback Gain vs Horizon (Q2)
    figure('Name', sprintf('%s – Gain vs H', sys(pp).name), ...
           'Position', [100 + 300*(pp-1), 400, 700, 420]);
    hold on;
    
    for rr = 1:numel(r_values)
        r_val = r_values(rr);
        
        % Plot RH gain
        plot(h_range, k_rh_data(rr, :), '-', ...
            'Color', color_map(rr, :), 'LineWidth', 1.8, ...
            'DisplayName', sprintf('K_{RH},  R = %.2g', r_val));
        
        % Plot LQ gain reference line
        plot(h_range, k_lq_vec(rr) * ones(1, h_max), '--', ...
            'Color', color_map(rr, :), 'LineWidth', 1.4, ...
            'DisplayName', sprintf('K_{LQ},  R = %.2g', r_val));
    end
    
    xlabel('Horizon H', 'FontSize', 14);
    ylabel('Feedback gain K', 'FontSize', 14);
    title([sys(pp).name ' – RH gain vs H (Q = 1)'], 'FontSize', 13);
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 13);
    xlim([1, h_max]);
    
    %% Plotting: Closed-Loop Eigenvalue vs Horizon (Q3)
    figure('Name', sprintf('%s – |CL eigenvalue| vs H', sys(pp).name), ...
           'Position', [100 + 300*(pp-1), 50, 700, 400]);
    hold on;
    
    for rr = 1:numel(r_values)
        plot(h_range, eig_data(rr, :), '-', ...
            'Color', color_map(rr, :), 'LineWidth', 1.8, ...
            'DisplayName', sprintf('R = %.2g', r_values(rr)));
    end
    
    % Stability boundary line
    yline(1, 'k-.', 'LineWidth', 2, 'DisplayName', 'Stability boundary |λ| = 1');
    
    xlabel('Horizon H', 'FontSize', 14);
    ylabel('|CL eigenvalue| |A – bK|', 'FontSize', 14);
    title([sys(pp).name ' – Closed-loop eigenvalue vs H (Q = 1)'], 'FontSize', 13);
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 13);
    xlim([1, h_max]);
    
    % Adjust y-axis limits dynamically based on data
    max_eig = max(eig_data(:));
    ylim([-0.05, max(1.3, 1.1 * max_eig)]);
    
    %% Q4: Discussion & Results Printout
    fprintf('\n--- Q4: Discussion (%s) ---\n', sys(pp).name);
    
    fprintf('\n  Minimum H required for closed-loop stability (|A - b*KRH| < 1):\n');
    for rr = 1:numel(r_values)
        r_val = r_values(rr);
        stable_idx = find(eig_data(rr, :) < 1, 1, 'first');
        
        if isempty(stable_idx)
            fprintf('    R = %5.2f  ->  NOT stabilised within H = 1..%d\n', r_val, h_max);
        else
            fprintf('    R = %5.2f  ->  minimum stable H = %d  (|lambda| = %.4f)\n', ...
                r_val, stable_idx, eig_data(rr, stable_idx));
        end
    end
end

%% Q5: Comparison Summary
fprintf('\n========================================\n');
fprintf('  Q5: Comparison summary\n');
fprintf('========================================\n');
fprintf('  Unstable plant requires a minimum H for stability.\n');
fprintf('  Stable plant is stable even at H = 1.\n');
fprintf('  Enlarging H is essential for the unstable plant.\n');
fprintf('\nDone.\n');