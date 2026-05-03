%% P2_solution.m
%  Basics on Receding-Horizon Control
%  -----------------------------------------------------------------------
%  Compares the Receding-Horizon (RH) state-feedback gain with the
%  infinite-horizon LQ gain for two first-order plants:
%    Plant 1 (unstable):  x(t+1) = 1.2 x(t) + u(t)
%    Plant 2 (stable):    x(t+1) = 0.8 x(t) + u(t)
%
%  For each plant the script:
%    Q1) Computes the optimal LQ state-feedback gain with dlqr.
%    Q2) Plots KRH vs horizon H, overlaid with the LQ gain.
%    Q3) Plots |eigenvalue of (A - b*K)| vs H with the stability boundary.
%    Q4) Discusses the results (printed to command window).
%    Q5) Repeats for the stable plant and compares.
%
%  Per the lab guide: Q = 1 (fixed), R is swept over representative
%  small and large values (explicitly required by the problem statement).
%
%  Requires: compute_KRH.m  (in the same folder)
%            MATLAB Control System Toolbox
%  -----------------------------------------------------------------------

clear; clc; close all;

%% ── User-defined parameters ─────────────────────────────────────────────
H_max  = 50;                    % maximum horizon to sweep
R_vals = [0.01, 0.1, 1, 10];   % representative small and large R values
                                % (required by the lab guide)
Q_cost = 1;                     % Q = 1 as specified by the lab guide

%% ── Plant definitions ───────────────────────────────────────────────────
plants(1).name = 'Unstable plant  (A = 1.2)';
plants(1).A    = 1.2;
plants(1).b    = 1;
plants(1).C    = 1;

plants(2).name = 'Stable plant  (A = 0.8)';
plants(2).A    = 0.8;
plants(2).b    = 1;
plants(2).C    = 1;

%% ─────────────────────────────────────────────────────────────────────────
%  Loop over the two plants
%% ─────────────────────────────────────────────────────────────────────────
for pp = 1:numel(plants)

    A = plants(pp).A;
    b = plants(pp).b;
    C = plants(pp).C;

    fprintf('\n========================================\n');
    fprintf('  %s\n', plants(pp).name);
    fprintf('========================================\n');

    % ── Q1  LQ gain ──────────────────────────────────────────────────────
    fprintf('\n--- Q1: LQ gain for each R (Q = %g) ---\n', Q_cost);

    Q_lq    = C' * Q_cost * C;   % scalar plant: Q_lq = 1
    KLQ_vec = zeros(size(R_vals));

    for rr = 1:numel(R_vals)
        R = R_vals(rr);
        [KLQ_temp, ~, ~] = dlqr(A, b, Q_lq, R);
        KLQ_vec(rr) = KLQ_temp;
        fprintf('  R = %5.2f  ->  KLQ = %8.5f,  CL eigenvalue = %+.5f\n', ...
            R, KLQ_temp, A - b*KLQ_temp);
    end

    % ── Compute RH gains and eigenvalues for all R and H ─────────────────
    H_range   = 1:H_max;
    colors    = lines(numel(R_vals));
    all_KRH_H = zeros(numel(R_vals), H_max);
    all_eig_H = zeros(numel(R_vals), H_max);

    for rr = 1:numel(R_vals)
        R = R_vals(rr);
        for hh = H_range
            KRH_temp         = compute_KRH(A, b, C, hh, R);
            all_KRH_H(rr,hh) = KRH_temp;
            all_eig_H(rr,hh) = abs(A - b*KRH_temp);
        end
    end

    % ── Figure Q2: gain vs H ─────────────────────────────────────────────
    figure('Name', sprintf('%s – Gain vs H', plants(pp).name), ...
           'Position', [100+300*(pp-1), 400, 700, 420]);
    hold on;

    for rr = 1:numel(R_vals)
        % RH gain curve
        plot(H_range, all_KRH_H(rr,:), '-', ...
            'Color', colors(rr,:), 'LineWidth', 1.8, ...
            'DisplayName', sprintf('K_{RH},  R = %.2g', R_vals(rr)));

        % FIX: each LQ dashed line gets its own legend entry in the same
        % colour as its RH curve, so the pairing is immediately clear.
        % Previously all dashed lines shared one dummy black legend entry.
        plot(H_range, KLQ_vec(rr)*ones(1,H_max), '--', ...
            'Color', colors(rr,:), 'LineWidth', 1.4, ...
            'DisplayName', sprintf('K_{LQ},  R = %.2g  (H\\to\\infty)', R_vals(rr)));
    end

    xlabel('Horizon H',       'FontSize', 14);
    ylabel('Feedback gain K', 'FontSize', 14);
    title([plants(pp).name '  – RH gain vs H  (Q = 1)'], 'FontSize', 13);
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 13);
    xlim([1, H_max]);

    % ── Figure Q3: |CL eigenvalue| vs H ──────────────────────────────────
    figure('Name', sprintf('%s – |CL eigenvalue| vs H', plants(pp).name), ...
           'Position', [100+300*(pp-1), 50, 700, 400]);
    hold on;

    for rr = 1:numel(R_vals)
        plot(H_range, all_eig_H(rr,:), '-', ...
            'Color', colors(rr,:), 'LineWidth', 1.8, ...
            'DisplayName', sprintf('R = %.2g', R_vals(rr)));
    end

    yline(1, 'k-.', 'LineWidth', 2, 'DisplayName', 'Stability boundary |λ| = 1');

    xlabel('Horizon H',                 'FontSize', 14);
    ylabel('|CL eigenvalue|  |A – bK|', 'FontSize', 14);
    title([plants(pp).name '  – Closed-loop eigenvalue vs H  (Q = 1)'], 'FontSize', 13);
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 13);
    xlim([1, H_max]);

    % FIX: ylim derived from the actual maximum across ALL R curves so no
    % curve is clipped. Previously the upper limit was computed from
    % abs(A - b*0) = |A|, which ignored what KRH actually produces.
    eig_max = max(all_eig_H(:));
    ylim([-0.05, max(1.3, 1.1 * eig_max)]);

    % ── Q4 Discussion ─────────────────────────────────────────────────────
    fprintf('\n--- Q4: Discussion (%s) ---\n', plants(pp).name);
    fprintf('  Q = %g  |  R values: %s\n', Q_cost, ...
        strjoin(arrayfun(@(r) sprintf('%.2g',r), R_vals, 'UniformOutput', false), ', '));

    fprintf('\n  KLQ (infinite-horizon target) for each R:\n');
    for rr = 1:numel(R_vals)
        fprintf('    R = %5.2f  ->  KLQ = %.5f\n', R_vals(rr), KLQ_vec(rr));
    end

    fprintf(['\n  * As H increases, KRH converges to KLQ for every R value.\n', ...
             '  * Small R -> larger, more aggressive gain; faster convergence with H.\n', ...
             '  * Large R -> smaller, gentler gain; slower convergence with H.\n']);

    % FIX: explicitly report the minimum H for stability for each R.
    % This directly answers Q3 and the Q4 discussion requirement.
    fprintf('\n  Minimum H required for closed-loop stability (|A - b*KRH| < 1):\n');
    for rr = 1:numel(R_vals)
        stable_idx = find(all_eig_H(rr,:) < 1, 1, 'first');
        if isempty(stable_idx)
            fprintf('    R = %5.2f  ->  NOT stabilised within H = 1..%d\n', ...
                R_vals(rr), H_max);
        else
            fprintf('    R = %5.2f  ->  minimum stable H = %d  (|lambda| = %.4f)\n', ...
                R_vals(rr), stable_idx, all_eig_H(rr, stable_idx));
        end
    end

    if abs(A) > 1
        fprintf(['\n  * UNSTABLE plant: a minimum H is NECESSARY for stability.\n', ...
                 '    With large R the gain is small, so more prediction steps\n', ...
                 '    are needed before KRH grows large enough to stabilise the loop.\n', ...
                 '    Enlarging H is therefore both a performance and a stability issue.\n']);
    else
        fprintf(['\n  * STABLE plant: the plant contracts on its own (|A| < 1), so\n', ...
                 '    even H = 1 gives a stable closed loop for all tested R values.\n', ...
                 '    Enlarging H refines the gain toward KLQ and reduces cost,\n', ...
                 '    but stability is not at risk.\n']);
    end

end  % plant loop

%% ── Q5 comparison summary ────────────────────────────────────────────────
fprintf('\n========================================\n');
fprintf('  Q5: Comparison summary\n');
fprintf('========================================\n');
fprintf([...
    '  Unstable plant (A = 1.2):\n', ...
    '    A sufficiently large H is NECESSARY: for small H (especially\n', ...
    '    with large R) the closed-loop eigenvalue exceeds 1 (unstable).\n', ...
    '    The gain must grow large enough to overcome the open-loop\n', ...
    '    instability, and this only happens once H is large enough.\n\n', ...
    '  Stable plant (A = 0.8):\n', ...
    '    Stability is achieved even at H = 1 for all R values tested.\n', ...
    '    The plant''s natural contraction (|A| < 1) keeps the loop stable\n', ...
    '    regardless of how small the gain is.\n', ...
    '    Enlarging H improves performance (gain converges to KLQ, cost\n', ...
    '    decreases) but is not required for stability.\n\n', ...
    '  Conclusion: enlarging H is far more advantageous -- and often\n', ...
    '  essential -- for the UNSTABLE plant.\n']);

fprintf('\nDone.\n');
