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
ts = truncate(1/fs, -2)

gn_of_z = filt(c2d(gn_of_s, ts))

[b, a, ts] = filtdata(gn_of_z);  
az = filt(a{1}, 1, ts)
bz = filt(b{1}, 1, ts)

desirable_s_poles = [-2 -2 -20 -20];
desirable_z_poles = e.^(desirable_s_poles.*ts)
controllable_poles = filt(zpk(desirable_z_poles(1:2), [0 0], 1,ts));