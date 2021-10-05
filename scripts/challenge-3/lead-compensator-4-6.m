addpath(genpath("./")) % Add lib path to Octave script file search paths
run common_parameters_scripts.m

% =============================================================================
% Controllers definitions
% =============================================================================

disp('=================== Open loop with lead compensator ===================')
[phimax, alpha, t, cl] = get_lead_compensator(k*g,
                                              desirable_margin = 60,
                                              slack = 24,
                                              gain = 1)
openloop = cl * k * g;
[gm, pm, w_gm, w_pm] = margin(openloop);
bandwidth = bandwidth_lti(openloop);
f = 1 / (t*s + 1)
disp('=======================================================================')

disp('============================== Compute Kc =============================')
[w2dot5_idx, w2dot5_approx] = find_idx(w, 2.5, 'ge');
[openloop_gains, w] = bodemag(openloop, w);
openloop_gain = openloop_gains(w2dot5_idx);
kc = 1 / (sqrt(2) * openloop_gain);
disp('=======================================================================')


disp(['=========== Open loop with fixed Kc equal ' num2str(kc) ' ==========='])
[phimax, alpha, t, cl_with_fixed_kc] = get_lead_compensator(kc*k*g,
                                                            desirable_margin = 60,
                                                            slack = 24, 
                                                            gain = 1);
phimax, alpha, t
cl_with_fixed_kc = kc * cl_with_fixed_kc
openloop_with_fixed_kc = cl_with_fixed_kc * k * g;
[gm, pm, w_gm, w_pm] = margin(openloop_with_fixed_kc);
bandwidth = bandwidth_lti(openloop_with_fixed_kc);
f_with_fixed_kc = 1 / (t*s + 1)
disp('=======================================================================')

