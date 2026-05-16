% Quick inspection of a saved open-loop dataset, using the same plot style
% as the TCLab open-loop acquisition script.
%
% Expects the .mat file to contain:
%   t : 1xN time vector              [s]
%   u : 2xN heater command history   [%]
%   y : 2xN measured temperatures    [°C]
%__________________________________________________________________________

clear all
close all
clc

% Load data file
filename = 'openloop_data_1.mat';
load(filename);   % expects variables: t (1xN), u (2xN), y (2xN)

% Plot in the same style as TCLab_openloop.m
figure

subplot(2,1,1), hold on, grid on
plot(t, y(1,:), '.', 'MarkerSize', 10)
plot(t, y(2,:), '.', 'MarkerSize', 10)

% The legend will only look at the plot commands above
legend('Temperature 1', 'Temperature 2', 'Location', 'best')
xlabel('Time [s]')
ylabel('Temperature [°C]')
title('Identification run - Temperature evolution')

% Add 'HandleVisibility', 'off' to exclude them from the legend
%xline(1000, '--r', 'HandleVisibility', 'off');
%xline(2000, '--r', 'HandleVisibility', 'off');

subplot(2,1,2), hold on, grid on
stairs(t, u(1,:), 'LineWidth', 2)
stairs(t, u(2,:), 'LineWidth', 2)
legend('Heater 1', 'Heater 2', 'Location', 'best')
xlabel('Time [s]')
ylabel('Heater control [%]')
ylim([0 100])
title('Identification run - Control input evolution')

%--------------------------------------------------------------------------
% End of File