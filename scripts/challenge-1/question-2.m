clear all, close all, clc
addpath(genpath("../../lib")) % Add lib path to Octave script file search paths

run common_functions_script
run graphs_functions_script

% =============================================================================
% Simulation parameters
% =============================================================================

sim = get_sim_time (dt = 0.01, end_time = 50);

% =============================================================================
% Input signals
% =============================================================================

reference = get_signal (sim.time, amplitude = 1,
                        start_time = 2, end_time = inf);

step_output_disturbance = get_signal (sim.time, amplitude = -0.2,
                                      start_time = 15, end_time = inf,
                                      signal = 'step');
step_input_disturbance = get_signal (sim.time, amplitude = -0.2,
                                     start_time = 25, end_time = inf,
                                     signal = 'step');

ramp_output_disturbance = get_signal (sim.time, slope = -0.2,
                                      start_time = 15, end_time = inf,
                                      signal = 'ramp');
ramp_input_disturbance = get_signal (sim.time, slope = -0.2,
                                     start_time = 25, end_time = inf,
                                     signal = 'ramp');

% =============================================================================
% Transfer Functions Definions
% =============================================================================

tau_n = 0;
tau_d = 2;
F = ((tau_n * s) + 1)/((tau_d * s) + 1);
K = 1;
z = 0.5;
G = 2 / s;
C = K * ((s + z) / s);

% =============================================================================
% Main of the script
% =============================================================================

y_to_r = get_output_to_reference_tf (G, C, F)
y_to_qy = get_output_to_output_disturbance_tf (G, C, F)
y_to_qu = get_output_to_input_disturbance_tf (G, C, F)

u_to_r = get_controller_output_to_reference_tf (G, C, F)
u_to_qy = get_controller_output_to_output_disturbance_tf (G, C, F)
u_to_qu = get_controller_output_to_input_disturbance_tf (G, C, F)

reference.response = lsim(y_to_r, reference.signal, sim.time);

step_output_disturbance.response = lsim(y_to_qy, step_output_disturbance.signal,sim.time);
step_input_disturbance.response = lsim(y_to_qu, step_input_disturbance.signal, sim.time);

ramp_output_disturbance.response = lsim(y_to_qy, ramp_output_disturbance.signal, sim.time);
ramp_input_disturbance.response = lsim(y_to_qu, ramp_input_disturbance.signal, sim.time);

reference.controller_output = lsim(u_to_r, reference.signal, sim.time);

step_output_disturbance.controller_output = lsim(u_to_qy, step_output_disturbance.signal, sim.time);
step_input_disturbance.controller_output = lsim(u_to_qu, step_input_disturbance.signal, sim.time);

ramp_output_disturbance.controller_output = lsim(u_to_qy, ramp_output_disturbance.signal, sim.time);
ramp_input_disturbance.controller_output = lsim(u_to_qu, ramp_input_disturbance.signal, sim.time);

control_loop_response_to_steps = (reference.response +
                                  step_output_disturbance.response +
                                  step_input_disturbance.response);

controller_output_signal_to_steps = (reference.controller_output +
                                     step_output_disturbance.controller_output +
                                     step_input_disturbance.controller_output);

control_loop_response_to_ramps = (reference.response +
                                  ramp_output_disturbance.response +
                                  ramp_input_disturbance.response);

controller_output_signal_to_ramps = (reference.controller_output +
                                     ramp_output_disturbance.controller_output +
                                     ramp_input_disturbance.controller_output);

% =============================================================================
% Plot Graphs
% =============================================================================

plot_responses_of_disturbances_signals (sim.time, figure_num = 1,
                                        reference.signal,
                                        control_loop_response_to_steps,
                                        controller_output_signal_to_steps,
                                        step_input_disturbance.signal,
                                        step_output_disturbance.signal)

plot_responses_of_disturbances_signals (sim.time, figure_num = 2,
                                        reference.signal,
                                        control_loop_response_to_ramps,
                                        controller_output_signal_to_ramps,
                                        ramp_input_disturbance.signal,
                                        ramp_output_disturbance.signal)
