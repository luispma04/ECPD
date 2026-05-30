%% Quick inspection of a saved open-loop dataset.
%--------------------------------------------------------------------------
% This script loads a previously saved open-loop dataset (for example,
% 'openloop_data_1.mat') and generates side-by-side subplots of the 
% measured temperatures and applied heater control signals over time. 
% It serves as a quick visualization tool to verify the integrity and 
% quality of the acquired experimental data before proceeding to the 
% system identification phase.
%
% Expected variables in the .mat file:
%   t : 1xN time vector            [s]
%   u : 2xN heater command history [%]
%   y : 2xN measured temperatures  [°C]
%
% Matlab Toolboxes: (None required beyond base MATLAB)
% Functions: (None custom required)
%--------------------------------------------------------------------------

%% Initialization

clear; clc; close all;

%% Load Data

% Define and load data file
file_name = 'openloop_data_1.mat';
load(file_name); % Expects variables: t (1xN), u (2xN), y (2xN)

%% Plot Results

figure;

% Temperature subplot
subplot(2, 1, 1); hold on; grid on;
plot(t, y(1, :), '.', 'MarkerSize', 10);
plot(t, y(2, :), '.', 'MarkerSize', 10);

legend('Temperature 1', 'Temperature 2', 'Location', 'best');
xlabel('Time [s]');
ylabel('Temperature [°C]');
title('Identification run - Temperature evolution');

% (Optional) Add vertical lines without affecting the legend
% xline(1000, '--r', 'HandleVisibility', 'off');
% xline(2000, '--r', 'HandleVisibility', 'off');

% Control input subplot
subplot(2, 1, 2); hold on; grid on;
stairs(t, u(1, :), 'LineWidth', 2);
stairs(t, u(2, :), 'LineWidth', 2);

legend('Heater 1', 'Heater 2', 'Location', 'best');
xlabel('Time [s]');
ylabel('Heater control [%]');
ylim([0 100]);
title('Identification run - Control input evolution');