addpath(genpath("./")) % Add current dir path to Octave script file search paths
run common_parameters_script.m

kf = 1/dcgain(bz/controllable_poles)
fr_of_z = kf*(bz/controllable_poles)

[fr_num_coefficients, fr_den_coefficients, ts] = filtdata(fr_of_z);
nfr_of_z = filt(fr_num_coefficients{1}, 1, ts)
dfr_of_z = filt(fr_den_coefficients{1}, 1, ts)

c_of_z = (kf*az)/(dfr_of_z - nfr_of_z)

closed_loop = minreal(feedback(c_of_z*gn_of_z, 1), minreal_precision)