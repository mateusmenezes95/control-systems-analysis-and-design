clear all, close all, clc
addpath(genpath("../../lib")) % Add lib path to Octave script file search paths

% =============================================================================
% External packages and script
% =============================================================================

run common_functions_script

pkg load miscellaneous
format short

% =============================================================================
% Model definition
% =============================================================================

g11 = 2/(10*s+1)
g12 = 0.8/((10*s+1)*(2*s+1))
g21 = 0.6/((10*s+1)*(2*s+1))
g22 = 1.5/(10*s+1)

wc11 = bandwidth_lti(g11)
wc12 = bandwidth_lti(g12)
wc21 = 0
wc22 = bandwidth_lti(g22)

gs = [g11 g12;
     g21 g22]
W = [wc11 wc12;
     wc21 wc22]

gs_ss = ss(gs)
det_A = det(gs_ss.A)

% =============================================================================
% Sampling period definition and discrete state space realization
% =============================================================================

wb_max = max(max(W))
fs = 30*(wb_max/(2*pi()))
sampling_period_computed = truncate(1/fs, -2)
sampling_period = 0.75

gz = c2d(gs, sampling_period);
gz_ss = c2d(gs_ss, sampling_period)

% =============================================================================
% Equilibrium output and control signal computation
% =============================================================================

g0 = dcgain(gs)
det_g0 = det(g0)
inv_g0 = inv(g0)
yeq = [2; 1];
ueq = inv_g0*yeq

% =============================================================================
% Equilibrium states
% =============================================================================

xeq_continue = -1*inv(gs_ss.A)*gs_ss.B*ueq
xeq_discrete = (inv(eye(size(gz_ss.A)) - gz_ss.A)*gz_ss.B)*ueq
yeq_discrete = gz_ss.C*xeq_discrete

% =============================================================================
% Simulation parameters
% =============================================================================

integration_step_ratio = 30;
integration_step_size = sampling_period/integration_step_ratio;
end_time = 140;
sim = get_sim_time (integration_step_size, ceil(end_time/sampling_period)*sampling_period);

% =============================================================================
% Input signals
% =============================================================================

reference = get_signal (sim.time, amplitude = 1,
                        start_time = 10, end_time = inf);
output_disturbance = get_signal (sim.time, amplitude = 0.2,
                                 start_time = 50, end_time = inf);
input_disturbance = get_signal (sim.time, amplitude = 0.2,
                                start_time = 100, end_time = inf);
qu = [output_disturbance.signal'; output_disturbance.signal'];
qy = [input_disturbance.signal'; input_disturbance.signal'];
unit_step = [reference.signal'; reference.signal'; reference.signal'; reference.signal'];
