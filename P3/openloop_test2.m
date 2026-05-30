%% Second open-loop experiment for data collection in TCLab.
%--------------------------------------------------------------------------
% This script applies a sequence of open-loop control signals to Heater 1
% and records the resulting temperature evolution. The test begins with a
% brief stabilization phase, followed by a rapid series of 10 step changes
% with varying amplitudes and durations. This specific sequence is designed
% to strongly excite the system's transient dynamics, generating a dataset
% that is highly suitable for system identification or model validation.
%
% Matlab Toolboxes: (None required beyond base MATLAB)
% Functions: tclab, T1C, T2C, h1, h2, led
%--------------------------------------------------------------------------

%% Initialization

clear; clc; close all;

% Initialize TCLab hardware
tclab;

%% Experiment Parameters

t_total = 1000;         % Experiment duration [s] (approx. 16m 40s)
Ts      = 5;            % Sampling period [s]
N       = t_total / Ts; % Number of samples to collect (200 samples)

%% Control Sequence Definition

% Pre-allocate input matrix (2 heaters x N samples)
u = zeros(2, N);

% Phase A: Wait 100 seconds (20 samples * 5s) to ensure steady state
u(1, 1:20) = 30;

% Phase B: Sequence of step changes with varying durations
u(1, 21:30)   = 10;   % -20 
u(1, 30:42)   = 75;   % +65 
u(1, 43:53)   = 30;   % -45 
u(1, 54:74)   = 50;   % +20 
u(1, 75:92)   = 20;   % -30 
u(1, 93:110)  = 90;   % +70 
u(1, 111:130) = 10;   % -80 
u(1, 131:150) = 60;   % +50 
u(1, 151:170) = 50;   % -10 
u(1, 171:200) = 30;   % -20 

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
exportgraphics(gcf, ['openloop_plot_2_', time_str, '.png'], 'Resolution', 300);
save(['openloop_data_2_', time_str, '.mat'], 'y', 'u', 't');

%--------------------------------------------------------------------------