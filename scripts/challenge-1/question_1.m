clear all, close all, clc
addpath(genpath("../../lib")) % Add lib path to Octave script file search paths

run common_functions_script
run graphs_functions_script

% ===============================================================================
% Simulation parameters
% ===============================================================================

sim = get_sim_time (dt = 0.01, end_time = 50);

% ===============================================================================
% Input signals
% ===============================================================================

reference = get_signal (sim.time, start_time = 2, end_time = inf, amplitude = 1);
output_disturbance = get_signal (sim.time, start_time = 15, end_time = inf, amplitude = -0.2);
input_disturbance = get_signal (sim.time, start_time = 25, end_time = inf, amplitude = -0.2);

% ===============================================================================
% Transfer Functions Definions
% ===============================================================================

F = 1;
K = 0.5;
G = 2 / s;
C = K;

% ===============================================================================
% Main of the script
% ===============================================================================

y_to_r = get_output_to_reference_tf (G, C, F)
y_to_qy = get_output_to_output_disturbance_tf (G, C, F)
y_to_qu = get_output_to_input_disturbance_tf (G, C, F)

u_to_r = get_controller_output_to_reference_tf (G, C, F)
u_to_qy = get_controller_output_to_output_disturbance_tf (G, C, F)
u_to_qu = get_controller_output_to_input_disturbance_tf (G, C, F)

reference.response = lsim(y_to_r, reference.signal, sim.time);
output_disturbance.response = lsim(y_to_qy, output_disturbance.signal, sim.time);
input_disturbance.response = lsim(y_to_qu, input_disturbance.signal, sim.time);

reference.controller_output = lsim(u_to_r, reference.signal, sim.time);
output_disturbance.controller_output = lsim(u_to_qy, output_disturbance.signal, sim.time);
input_disturbance.controller_output = lsim(u_to_qu, input_disturbance.signal, sim.time);

control_loop_response = (reference.response +
                         output_disturbance.response +
                         input_disturbance.response);

controller_output_signal = (reference.controller_output +
                            output_disturbance.controller_output +
                            input_disturbance.controller_output);

% ===============================================================================
% Plot Graphs
% ===============================================================================

plot_responses_of_disturbances_signals (sim.time, reference.signal,
                                        control_loop_response,
                                        controller_output_signal,
                                        input_disturbance.signal,
                                        output_disturbance.signal)
