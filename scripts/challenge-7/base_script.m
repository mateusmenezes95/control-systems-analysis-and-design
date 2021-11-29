clear all, close all, clc
addpath(genpath("../../lib")) % Add lib path to Octave script file search paths

% =============================================================================
% External packages and script
% =============================================================================

run common_functions_script
run graphs_functions_script

pkg load miscellaneous
format short

% =============================================================================
% Model definition
% =============================================================================

xt0 = [-5; -5];

% =============================================================================
% Equilibrium states
% =============================================================================

xt_eq(:,1:2) = [1 4; 0 0]
ut_eq(1:2) = [1 4];
yt_eq(1:2) = xt_eq(1,:);

% =============================================================================
% Simulation parameters
% =============================================================================

integration_step_size = 0.01;
end_time = 60;
sim = get_sim_time (integration_step_size, end_time);
sim_time_length = length(sim.time);

% =============================================================================
% Input signals
% =============================================================================

reference = get_signal (sim.time, amplitude = 1,
                        start_time = 1, end_time = inf);
output_disturbance = get_signal (sim.time, amplitude = 0.2,
                                 start_time = 10, end_time = inf);
input_disturbance = get_signal (sim.time, amplitude = 0.2,
                                start_time = 20, end_time = inf);
qu = [output_disturbance.signal'; output_disturbance.signal'];
qy = [input_disturbance.signal'; input_disturbance.signal'];
unit_step = [reference.signal'; reference.signal'; reference.signal'; reference.signal'];
