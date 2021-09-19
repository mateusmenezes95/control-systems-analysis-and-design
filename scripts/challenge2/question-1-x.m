clear all, close all, clc
addpath(genpath("../../lib")) % Add lib path to Octave script file search paths

run common_functions_script
run graphs_functions_script

clear s;

% =============================================================================
% Simulation parameters
% =============================================================================

sim = get_sim_time (dt = 0.01, end_time = 60);

% =============================================================================
% Input signals
% =============================================================================

reference = get_signal (sim.time, amplitude = 1,
                        start_time = 2, end_time = inf);

step_output_disturbance = get_signal (sim.time, amplitude = -0.2,
                                      start_time = 20, end_time = inf,
                                      signal = 'step');
step_input_disturbance = get_signal (sim.time, amplitude = -0.2,
                                     start_time = 40, end_time = inf,
                                     signal = 'step');

sine_output_disturbance.signal = get_sine_signal (sim.time,
                                                  amplitude = -0.2, natural_freq = 2,
                                                  start_time = 20, end_time = inf);
sine_input_disturbance.signal = get_sine_signal (sim.time,
                                                 amplitude = -0.2, natural_freq = 2,
                                                 start_time = 40, end_time = inf);

% =============================================================================
% Transfer Functions Definions
% =============================================================================

w = logspace(-2, 4, 1e4);
s = j*w;

l = [0.9 0.7 0.6 0.4];
k = [1.3 0.9 1.2 0.8];
tau = [1.2 1.1 0.8 0.9];

ln = mean(l);
kn = mean(k);
taun = mean(tau);

gn = kn .* exp(-s * ln) ./ (taun * s + 1);
g = zeros(length(l), length(gn));
delta = zeros(length(l), length(gn));

tauc = 0.5;
kc = taun/(kn * (tauc + ln));
ti = min(taun, 4 * (tauc + ln))

c = kc * ((s*ti + 1) ./ s*ti);

% =============================================================================
% Main of the script
% =============================================================================

for i = 1:length(l)
    g(i, :) = k(i) .* exp(-s * l(i)) ./ (tau(i) * s + 1);
    delta(i, :) = (abs(gn) - abs(g(i, :))) ./ abs(gn);
    % delta(i, :) = abs((gn - g(i, :)) ./ gn);
    figure(i);
    semilogx(w, delta(i, :));
endfor

figure(i+1)
delta_th = max(delta);
semilogx(w, delta_th)

figure(i+2)

comp_sensibility = (c .* gn) ./ (1 + c .* gn);
cdelta_th = abs(comp_sensibility) .* delta_th;
semilogx(w, 20*log10(cdelta_th))
