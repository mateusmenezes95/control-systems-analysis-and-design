clear all, close all, clc
addpath(genpath("../../lib")) % Add lib path to Octave script file search paths

run common_functions_script
run graphs_functions_script

clear s;

% =============================================================================
% Simulation parameters
% =============================================================================
dt = 0.01;
sim = get_sim_time (dt, end_time = 60);

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

num_plants = length(l);

ln = mean(l);
kn = mean(k);
taun = mean(tau);

gn = kn .* exp(-s * ln) ./ (taun * s + 1);
g = zeros(num_plants, length(gn));
delta = zeros(num_plants, length(gn));

tauc = 0.5;
kc = taun/(kn * (tauc + ln));
ti = min(taun, 4 * (tauc + ln))

c = kc * ((s*ti + 1) ./ s*ti);

f = 1;

% =============================================================================
% Main of the script
% =============================================================================

axs = zeros(1,num_plants);

for i = 1:num_plants
    g(i, :) = k(i) .* exp(-s * l(i)) ./ (tau(i) * s + 1);
    delta(i, :) = abs((g(i, :) - gn) ./ gn);
    axs(i) = subplot(num_plants,1,i);
    semilogx(w, delta(i, :));
endfor


lm = max(delta);
linkaxes(axs);
ylim([min(lm) 2.5])

figure("name", "Análise da Incerteza Multiplicativa")
subplot(3,1,1)
semilogx(w, lm)

comp_sensibility = (c .* gn) ./ (1 + c .* gn);
subplot(3,1,2)
semilogx(w, abs(comp_sensibility))

comp_sensibility_x_lm = abs(comp_sensibility) .* lm;
subplot(3,1,3)
semilogx(w, comp_sensibility_x_lm)

s = tf('s');
sampling_period = 1 / 100;

c = kc * ((s*ti + 1) / s*ti);

for i = 1:num_plants
    g = k(i) / (tau(i) * s + 1)
    comp_sensibility = (c * g)/(1 + (c * g));
endfor

[u, y, error_] = simulate_sys(sim.time, dt,
                      g, 0.4, c, f,
                      reference.signal,
                      step_output_disturbance.signal,
                      step_input_disturbance.signal);

subplot(3,1,1)
plot(sim.time, y)
subplot(3,1,2)
plot(sim.time, u)
subplot(3,1,3)
plot(sim.time, error_)

