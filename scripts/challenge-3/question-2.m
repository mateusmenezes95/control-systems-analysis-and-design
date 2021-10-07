addpath(genpath("./")) % Add lib path to Octave script file search paths
run common_parameters_scripts.m

% =============================================================================
% Main of the script
% =============================================================================
p = k * g
p_wb = bandwidth_lti(p)
bodemag(p, w)
[p_gain, w] = bodemag(p, w);

set(legend, 'string', '$|P(jw)|$')
plot_coordinates(p_wb, -3.0, [1 1 1], '', 'r', '--', 'o')
text(p_wb + 0.2, -3, ['$|P(w_{b} \approx  ' num2str(p_wb, '%.2f') ')| \approx -3_{db}$'])
xlabel('Frequência [rad/s]')
title('')
set(legend, 'location', 'southwest')