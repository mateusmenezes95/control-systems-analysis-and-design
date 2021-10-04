addpath(genpath("./")) % Add lib path to Octave script file search paths
run common_parameters_scripts.m

% =============================================================================
% Main of the script
% =============================================================================

disp('================= Lead compensator with slash = 24 deg ================')
[phimax, alpha, t, cl24slash] = get_lead_compensator(k*g, desirable_margin = 60, slack = 24)
openloop24slash = cl24slash * k * g;
[gain_margin, phase_margin, w_gain_margin, w_phase_margin] = margin(openloop24slash)

[w2dot5_idx, w2dot5_approx] = find_idx(w, 2.5, 'ge')
[openloop24slash_gains, w] = bodemag(openloop24slash, w);
openloop24slash_gain = openloop24slash_gains(w2dot5_idx)

kc = 1 / (sqrt(2) * openloop24slash_gain)

[phimax, alpha, t, cl24slash_with_kc] = get_lead_compensator(k*g, desirable_margin = 60, slack = 24, kc)
openloop24slash_with_kc = cl24slash_with_kc * k * g;
[gain_margin, phase_margin, w_gain_margin, w_phase_margin] = margin(openloop24slash_with_kc)
disp('=======================================================================')

% =============================================================================
% Plot Graphs
% =============================================================================

hold on
plot_bode (k * g, w, 'b', '-', '$\overline{K}G(jw)$', [2 2 2], [2 2 4])
plot_phase_margin(k * g, [2 2 2], [2 2 4], 'b')
hold on
plot_bode (openloop24slash, w, 'r', '--', '$C(jw)\overline{K}G(jw)$ p/ kc = 1', [2 2 2], [2 2 4])
plot_phase_margin(openloop24slash, [2 2 2], [2 2 4], 'r')
hold on
plot_bode (openloop24slash_with_kc, w, 'k', '--',
          '$C(jw)\overline{K}G(jw)$ p/ kc \approx 1.43', [2 2 2], [2 2 4])
plot_phase_margin(openloop24slash_with_kc, [2 2 2], [2 2 4], 'k')

hold on
plot_bode (cl24slash, w, 'r', '-', '$C(jw) p/ kc = 1°', [2 2 1], [2 2 3])
plot_phase_peak(cl24slash, w, [2 2 1], [2 2 3], 'r');
hold on
plot_bode (cl24slash_with_kc, w, 'k', '-', '$C(jw)$ p/ kc \appox 1.43', [2 2 1], [2 2 3])
plot_phase_peak(cl24slash_with_kc, w, [2 2 1], [2 2 3], 'k');
