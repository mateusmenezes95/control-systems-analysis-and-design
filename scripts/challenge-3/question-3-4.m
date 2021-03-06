addpath(genpath("./")) % Add lib path to Octave script file search paths
run common_parameters_scripts.m

% =============================================================================
% Main of the script
% =============================================================================

disp('============================ KG(jw) margins ===========================')
[gain_margin, phase_margin, w_gain_margin, w_phase_margin] = margin(k*g)
disp('=======================================================================')

disp('================= Lead compensator with slash = 12 deg ================')
[phimax, alpha, t, cl12slash] = get_lead_compensator(k*g, desirable_margin = 60, slack = 12)
openloop12slash = cl12slash * k * g;
[gain_margin, phase_margin, w_gain_margin, w_phase_margin] = margin(openloop12slash)
wb12slash = bandwidth_lti(openloop12slash)
disp('=======================================================================')

disp('================= Lead compensator with slash = 24 deg ================')
[phimax, alpha, t, cl24slash] = get_lead_compensator(k * g, desirable_margin = 60, slack = 24)
openloop24slash = cl24slash * k * g;
wb24slash = bandwidth_lti(openloop24slash)
[gain_margin, phase_margin, w_gain_margin, w_phase_margin] = margin(openloop24slash)
disp('=======================================================================')

% =============================================================================
% Plot Graphs
% =============================================================================
magnitude_subplot = [2 1 1];
fase_subplot = [2 1 2];

figure(1)
plot_bode (k * g, w, 'b', '-', "$\\overline{K}G(jw)$", magnitude_subplot, fase_subplot)
plot_phase_margin(k * g, magnitude_subplot, fase_subplot, 'b')

plot_bode (openloop12slash, w, 'r', '--', "$C(jw)G(jw)$ p/ folga $= 12^{\\circ}$", magnitude_subplot, fase_subplot)
plot_phase_margin(openloop12slash, magnitude_subplot, fase_subplot, 'r')

plot_bode (openloop24slash, w, 'k', '--', "$C(jw)G(jw)$ p/ folga $= 24^{\\circ}$", magnitude_subplot, fase_subplot)
plot_phase_margin(openloop24slash, magnitude_subplot, fase_subplot, 'k')

figure(2)
plot_bode (cl12slash, w, 'r', '-', "${C}\'(jw)$ p/ folga $= 12^{\\circ}$", [2 1 1], [2 1 2])
plot_phase_peak(cl12slash, w, [2 1 1], [2 1 2], 'r')

plot_bode (cl24slash, w, 'k', '-', "${C}\'(jw)$ p/ folga $= 24^{\\circ}$", [2 1 1], [2 1 2])
plot_phase_peak(cl24slash, w, [2 1 1], [2 1 2], 'k');
set(legend, 'location', 'southwest')

subplot(2,1,1)
set(legend, 'location', 'northwest')
