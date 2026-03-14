clearvars; close all;

q = tf('q');
F = (0.7157 + 1.4315*q^-1 + 0.7157*q^-2)/(1 + 1.3490*q^-1 + 0.5140*q^-2);

%% Part 1: Understanding saturation and Butterworth filter
% 1.1 Determine saturation parameter M by experiment
% 100 second input sweep from -5 to 5
r = linspace(-,5,100);

[u,y] = assignment_sys_18(r,'open loop');

M = max(u); % == -1 * min(u)

% The combination of actuator saturation and Butterworth filter limits the
% amount of energy that can be applied to the system in the high frequency
% range. This will reduce the signal-to-noise ratio in that range.

%% Part 2: Nonparametric identification
% 2.1 Frequency behavior of G_0
% Reference that does