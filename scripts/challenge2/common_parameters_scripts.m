clear all, close all, clc
addpath(genpath("../../lib")) % Add lib path to Octave script file search paths

run common_functions_script
run graphs_functions_script

% =============================================================================
% Simulation parameters
% =============================================================================
dt = 0.01;
sim = get_sim_time (dt, end_time = 60);
sim_discrete_size = length(sim.time) 

% =============================================================================
% Input signals
% =============================================================================

reference = get_signal (sim.time, amplitude = 1,
                        start_time = 2, end_time = inf);

step_output_disturbance = get_signal (sim.time, amplitude = 0.2,
                                      start_time = 25, end_time = inf,
                                      signal = 'step');
step_input_disturbance.signal = zeros(1, length(reference.signal));

% =============================================================================
% Transfer Functions Definions
% =============================================================================

l = [0.9 0.7 0.6 0.4];
k = [1.3 0.9 1.2 0.8];
tau = [1.2 1.1 0.8 0.9];

num_plants = length(l)

ln = mean(l)
kn = mean(k)
taun = mean(tau)

tauc = 0.5
kc = taun/(kn * (tauc + ln))
ti = min(taun, 4 * (tauc + ln))
