clear all, close all, clc
addpath(genpath("../../lib")) % Add lib path to Octave script file search paths

run common_functions_script
run graphs_functions_script

% =============================================================================
% Transfer Functions Definions
% =============================================================================

g = 0.5 / ((s^2 + 0.6*s + 1)*(0.1*s + 1))
k = 1 / g(0)
w = logspace(-1, 2, 1e4);
