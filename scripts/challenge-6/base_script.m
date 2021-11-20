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
% Simulation parameters
% =============================================================================

sampling_period = 0.1;
integration_step_ratio = 20;
dt = sampling_period/integration_step_ratio;
end_time = 20;
sim = get_sim_time (dt, round(end_time/sampling_period)*sampling_period);
minreal_precision = 0.01;

% =============================================================================
% Input signals
% =============================================================================

reference = get_signal (sim.time, amplitude = 1,
                        start_time = 1, end_time = inf);
output_disturbance = get_signal (sim.time, amplitude = 0.2,
                                 start_time = 7, end_time = inf);
input_disturbance = get_signal (sim.time, amplitude = 0.2,
                                start_time = 12, end_time = inf);

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

Gs = [g11 g12;
     g21 g22]
W = [wc11 wc12;
     wc21 wc22]

Gs_ss = ss(Gs)
det_A = det(Gs_ss.A)

% =============================================================================
% Sampling period definition and discrete state space realization
% =============================================================================

wb_max = max(max(W))
fs = 30*(wb_max/(2*pi()))
sampling_period_computed = truncate(1/fs, -2)
sampling_period = 0.75

Gz = c2d(Gs, sampling_period);
Gz_ss = ss(Gz)

% =============================================================================
% Equilibrium output and control signal computation
% =============================================================================

g0 = dcgain(G)
det_g0 = det(g0)
inv_g0 = inv(g0)
yeq = [2; 1];
ueq = inv_g0*yeq

% =============================================================================
% Equilibrium states
% =============================================================================

xeq = -1*inv(Gs_ss.A)*Gs_ss.B*ueq