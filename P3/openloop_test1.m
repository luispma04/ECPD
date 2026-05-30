%% First open-loop experiment for data collection in TCLab.
%--------------------------------------------------------------------------
% This script applies a sequence of open-loop control signals to Heater 1
% and records the resulting temperature evolution. The test begins with a
% brief stabilization phase, followed by a series of step changes with
% varying amplitudes and durations. This specific sequence is designed to
% excite a wide range of system dynamics for robust system identification.
%
% Matlab Toolboxes: (None required beyond base MATLAB)
% Functions: tclab, T1C, T2C, h1, h2, led
%--------------------------------------------------------------------------
%% Initialization

clear; clc; close all;

% Initialize TCLab hardware
tclab;

%% Experiment Parameters

t_total = 8000;         % Experiment duration [s] (approx. 2h 13m)
Ts      = 5;            % Sampling period [s]
N       = t_total / Ts; % Number of samples to collect (1600 samples)

%% Control Sequence Definition

% Pre-allocate input matrix (2 heaters x N samples)
u = zeros(2, N);

% Phase A: Wait 2000 seconds (400 samples * 5s) to ensure steady state
u(1, 1:400) = 45;

% Phase B: 5 incremental steps, each held for ~200 samples (~1000s)
u(1, 401:600)   = 48;   % +3
u(1, 601:800)   = 45;   % -3
u(1, 801:1000)  = 43;   % -2
u(1, 1001:1200) = 51;   % +8
u(1, 1201:1400) = 40;   % -11
u(1, 1401:1600) = 45;   % Return to baseline

%% Experiment Setup

% Real-time plot flag
% true  = plots input and measured temperature in real time
% false = plots at the end and prints results in command window
rt_plot = true;

if rt_plot
    figure;
    drawnow;
end

% Pre-allocate arrays for time and temperature measurements
t = nan(1, N);
y = nan(2, N);

% String with date for uniquely saving results
time_str = char(datetime('now', 'Format', 'yyMMdd_HHmmSS'));

%% Experiment Execution

% Signal the start of the experiment by lighting the LED
led(1);
disp('Temperature test started.');

for k = 1:N
    tic;

    % Compute elapsed time
    t(k) = (k - 1) * Ts;

    % Read sensor temperatures
    y(1, k) = T1C();
    y(2, k) = T2C();

    % Apply control variables to the plant
    h1(u(1, k));
    h2(u(2, k));

    if rt_plot
        % Update real-time plots
        clf;
        
        % Temperature subplot
        subplot(2, 1, 1); hold on; grid on;
        plot(t(1:k), y(1, 1:k), '.', 'MarkerSize', 10);
        plot(t(1:k), y(2, 1:k), '.', 'MarkerSize', 10);
        legend('Temperature 1', 'Temperature 2', 'Location', 'northwest');
        xlabel('Time [s]');
        ylabel('Temperature [°C]');
        
        % Control input subplot
        subplot(2, 1, 2); hold on; grid on;
        stairs(t(1:k), u(1, 1:k), 'LineWidth', 2);
        stairs(t(1:k), u(2, 1:k), 'LineWidth', 2);
        legend('Heater 1', 'Heater 2', 'Location', 'northwest');
        xlabel('Time [s]');
        ylabel('Heater [%]');
        ylim([0 100]);
        
        drawnow;
    else
        % Print to command window instead of plotting
        fprintf('t = %d, y1 = %.1f C, y2 = %.1f C, u1 = %.1f, u2 = %.1f\n', ...
                t(k), y(1, k), y(2, k), u(1, k), u(2, k));
    end

    % Check for timing violations
    if toc > Ts
        warning('Computation time exceeded sampling time by %.2f s at sample %d.', toc - Ts, k);
    end
    
    % Wait for the remainder of the sampling interval
    pause(max(0, Ts - toc));
end

%% Finalization and Plotting

% Turn off both heaters
h1(0);
h2(0);

% Signal the end of the experiment
led(0);
disp('Temperature test complete.');

% Plot final results if real-time plotting was disabled
if ~rt_plot
    figure;
    
    subplot(2, 1, 1); hold on; grid on;
    plot(t, y(1, :), '.', 'MarkerSize', 10);
    plot(t, y(2, :), '.', 'MarkerSize', 10);
    legend('Temperature 1', 'Temperature 2', 'Location', 'best');
    xlabel('Time [s]');
    ylabel('Temperature [°C]');
    
    subplot(2, 1, 2); hold on; grid on;
    stairs(t, u(1, :), 'LineWidth', 2);
    stairs(t, u(2, :), 'LineWidth', 2);
    legend('Heater 1', 'Heater 2', 'Location', 'best');
    xlabel('Time [s]');
    ylabel('Heater control [%]');
    ylim([0 100]);
end

%% Save Results

% Save figure and experiment data to file
exportgraphics(gcf, ['openloop_plot_', time_str, '.png'], 'Resolution', 300);
save(['openloop_data_1_', time_str, '.mat'], 'y', 'u', 't');

%--------------------------------------------------------------------------