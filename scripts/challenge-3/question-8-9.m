addpath(genpath("./")) % Add lib path to Octave script file search paths
run lead-compensator-4-6.m

% =============================================================================
% Simulation parameters
% =============================================================================

sim = get_sim_time (dt = 0.01, end_time = 200);

% =============================================================================
% Input signals
% =============================================================================

r = get_signal (sim.time, amplitude = 1, start_time = 2, end_time = inf);
qu = get_signal (sim.time, amplitude = 1, start_time = 120, end_time = inf);
qy = get_signal (sim.time, amplitude = 0, start_time = 0, end_time = inf);

% =============================================================================
% Main of the script
% =============================================================================

tlag = 5 * tlead
tlag_with_fixed_kc = 5 * tlead_with_fixed_kc

[betarvar, clag] = get_lag_compensator(cl*k*g, openloop_gain = 5, tlag)
[beta_with_fixed_kc, clag_with_fixed_kc] = get_lag_compensator(cl_with_fixed_kc*k*g,
                                                                openloop_gain = 5,
                                                                tlag_with_fixed_kc)

c = minreal(cl*k*clag)
c_with_fixed_kc = minreal(cl_with_fixed_kc*k*clag_with_fixed_kc)

[u, y, err] = simulate_sys_lsim(sim.time, g, c, f, r.signal, qy.signal, qu.signal);
[u_with_fixed_kc, y_with_fixed_kc, err_with_fixed_kc] = simulate_sys_lsim(sim.time,
    g, c_with_fixed_kc, f_with_fixed_kc, r.signal, qy.signal, qu.signal);

% =============================================================================
% Plot Graphs
% =============================================================================
unit_kc_str = ' p/ $K_{c} = 1$';
plot_response_and_control_signals (sim.time, figure_num = 1,
                                    r.signal,
                                    y, ['$y(t)$' unit_kc_str],
                                    err, ['$e(t)$' unit_kc_str],
                                    u, ['$u(t)$' unit_kc_str], 'b-')
                                    
hold on
fixed_kc_str = [' p/ $K_{c} \approx ' num2str(kc, '%.2f') '$'];
plot_response_and_control_signals (sim.time, figure_num = 1,
                                    'none',
                                    y_with_fixed_kc, ['$y(t)$' fixed_kc_str],
                                    err_with_fixed_kc, ['$e(t)$' fixed_kc_str],
                                    u_with_fixed_kc, ['$u(t)$' fixed_kc_str], 'm-')
