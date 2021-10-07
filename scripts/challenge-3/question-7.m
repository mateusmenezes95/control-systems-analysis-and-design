addpath(genpath("./")) % Add lib path to Octave script file search paths
run lead-compensator-4-6.m

% =============================================================================
% Simulation parameters
% =============================================================================

sim = get_sim_time (dt = 0.01, end_time = 50);

% =============================================================================
% Input signals
% =============================================================================

r = get_signal (sim.time, amplitude = 1, start_time = 2, end_time = inf);
qu = get_signal (sim.time, amplitude = 0, start_time = 0, end_time = inf);
qy = get_signal (sim.time, amplitude = 0, start_time = 0, end_time = inf);

% =============================================================================
% Main of the script
% =============================================================================

p = k*g;
c = cl*k;
c_with_fixed_kc = cl_with_fixed_kc*k;

[u, y, err] = simulate_sys_lsim(sim.time, p, cl, f, r.signal, qy.signal, qu.signal);
[u_with_fixed_kc, y_with_fixed_kc, err_with_fixed_kc] = simulate_sys_lsim(sim.time,
    p, cl_with_fixed_kc, f_with_fixed_kc, r.signal, qy.signal, qu.signal);

mag = bodemag(cl*k*g, w);
mag_with_fixed_kc = bodemag(cl_with_fixed_kc*k*g, w);

err_inf = 1 / (1 + c(0)*g(0))
y_inf = 1 - err_inf
err_with_fixed_kc_inf = 1 / (1 + c_with_fixed_kc(0)*g(0))
y_with_fixed_kc_inf = 1 - err_with_fixed_kc_inf

% =============================================================================
% Plot Graphs
% =============================================================================
figure(1)

plot_signal(sim.time, '', r.signal, '', '--r', '$r(t)$')
hold on

unit_kc_str = ' p/ $K_{c} = 1$';
plot_signal (sim.time, '',
             y, '',
             line_color = 'b-',
             ['$y(t)$' unit_kc_str])
                                    
hold on
fixed_kc_str = [' p/ $K_{c} = ' num2str(kc, '%.2f') '$'];
plot_signal (sim.time, 'Tempo(s)',
             y_with_fixed_kc, 'Saída $y(t)$',
             line_color = 'r-',
             ['$y(t)$' fixed_kc_str])

figure(2)
plot_semilogx(w, 20*log10(mag), 'b', '-', [open_loop_str unit_kc_str], '', '')
hold on
plot_semilogx(w, 20*log10(mag_with_fixed_kc), 'r', '-', [open_loop_str fixed_kc_str],
              'Frequência [rad/s]', 'Ganho [db]')
set(legend, 'location', 'southwest')
