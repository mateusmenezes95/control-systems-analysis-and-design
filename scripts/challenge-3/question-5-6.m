addpath(genpath("./")) % Add lib path to Octave script file search paths
run common_parameters_scripts.m

% =============================================================================
% Main of the script
% =============================================================================

disp('=================== Open loop with lead compensator ===================')
[phimax, alpha, t, cl] = get_lead_compensator(k*g,
                                              desirable_margin = 60,
                                              slack = 24,
                                              gain = 1)
openloop = cl * k * g;
[gm, pm, w_gm, w_pm] = margin(openloop)
bandwidth = bandwidth_lti(openloop)
disp('=======================================================================')

disp('============================== Compute Kc =============================')
[w2dot5_idx, w2dot5_approx] = find_idx(w, 2.5, 'ge');
[openloop_gains, w] = bodemag(openloop, w);
openloop_gain = openloop_gains(w2dot5_idx)
kc = 1 / (sqrt(2) * openloop_gain)
disp('=======================================================================')

disp('===================== Open loop with wb equal 2.5 =====================')
[phimax, alpha, t, cl_with_2dot5_wb] = get_lead_compensator(k*g,
                                                            desirable_margin = 60,
                                                            slack = 24, 
                                                            gain = kc)
openloop_with_2dot5_wb = cl_with_2dot5_wb * k * g;
[gm, pm, w_gm, w_pm] = margin(openloop_with_2dot5_wb)
bandwidth = bandwidth_lti(openloop_with_2dot5_wb)
disp('=======================================================================')

disp(['=========== Open loop with fixed Kc equal ' num2str(kc) ' ==========='])
[phimax, alpha, t, cl_with_fixed_kc] = get_lead_compensator(kc*k*g,
                                                            desirable_margin = 60,
                                                            slack = 24, 
                                                            gain = 1);
phimax, alpha, t
cl_with_fixed_kc = kc * cl_with_fixed_kc 
openloop_with_fixed_kc = cl_with_fixed_kc * k * g;
[gm, pm, w_gm, w_pm] = margin(openloop_with_fixed_kc)
bandwidth = bandwidth_lti(openloop_with_fixed_kc)
disp('=======================================================================')

% =============================================================================
% Plot Graphs
% =============================================================================

cl_str = "${C}\'(jw)$";

figure(1)

hold on
plot_bode (openloop, w, 'b', '-', "$C(jw)G(jw)$ p/ folga $= 24^{\\circ}$",
          [2 1 1], [2 1 2])
plot_phase_margin(openloop, [2 1 1], [2 1 2], 'b')
hold on
plot_bode (openloop_with_2dot5_wb, w, 'r', '-',
          [open_loop_str ' c/ $w_{b}$ = 2.5'], [2 1 1], [2 1 2])
plot_phase_margin(openloop_with_2dot5_wb, [2 1 1], [2 1 2], 'r')
hold on
plot_bode (openloop_with_fixed_kc, w, 'k', '-',
          [open_loop_str ' c/ $K_{c}$ fixo'], [2 1 1], [2 1 2])
plot_phase_margin(openloop_with_fixed_kc, [2 1 1], [2 1 2], 'k')

subplot(2,1,1)
% line(xlim(), [-3 -3], 'color', [0.15 0.15 0.15], 'linewidth', 0.1, 'linestyle', '-')

figure(2)

hold on
plot_bode (cl_with_2dot5_wb, w, 'r', '-',
          [cl_str ' c/ $w_{b}$ = 2.5'], [2 1 1], [2 1 2])
plot_phase_peak(cl_with_2dot5_wb, w, [2 1 1], [2 1 2], 'r');
hold on
plot_bode (cl_with_fixed_kc, w, 'k', '-',
          [cl_str ' c/ $K_{c}$ fixo'], [2 1 1], [2 1 2])
plot_phase_peak(cl_with_fixed_kc, w, [2 1 1], [2 1 2], 'k');
