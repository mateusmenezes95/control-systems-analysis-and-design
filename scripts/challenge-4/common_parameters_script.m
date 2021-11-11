clear all, close all, clc
addpath(genpath("../../lib")) % Add lib path to Octave script file search paths

run common_functions_script
run graphs_functions_script

pkg load miscellaneous
format short

minreal_precision = 0.01;

gn_of_s = minreal((0.2*(10-s))/(s+1)^2)

wb = bandwidth_lti(gn_of_s)
fs = 30*(wb/(2*pi()))
sampling_period = truncate(1/fs, -2)

gn_of_z = filt(c2d(gn_of_s, sampling_period))

[b, a, sampling_period] = filtdata(gn_of_z);  
az = filt(a{1}, 1, sampling_period)
bz = filt(b{1}, 1, sampling_period)

desirable_s_poles = [-2 -2 -20 -20];
desirable_z_poles = e.^(desirable_s_poles.*sampling_period)
controllable_poles = filt(zpk(desirable_z_poles(1:2), [0 0], 1,sampling_period));

% =============================================================================
% Simulation parameters
% =============================================================================
integration_step_ratio = 20;
dt = sampling_period/integration_step_ratio;
sim = get_sim_time (dt, end_time = 20);

% =============================================================================
% Input signals
% =============================================================================

reference = get_signal (sim.time, amplitude = 1,
                        start_time = 1, end_time = inf);
output_disturbance = get_signal (sim.time, amplitude = 0.2,
                                 start_time = 7, end_time = inf);
input_disturbance = get_signal (sim.time, amplitude = 0.2,
                                start_time = 12, end_time = inf);
