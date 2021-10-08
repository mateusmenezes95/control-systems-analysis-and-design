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

[gm, pm, w_gm, w_pm] = margin(c*g);
[gm_with_fixed_kc, pm_with_fixed_kc, w_gm_with_fixed_kc, w_pm_with_fixed_kc] = margin(
    c_with_fixed_kc*g);

[u, y, err] = simulate_sys_lsim(sim.time, g, c, f, r.signal, qy.signal, qu.signal);
[u_with_fixed_kc, y_with_fixed_kc, err_with_fixed_kc] = simulate_sys_lsim(sim.time,
    g, c_with_fixed_kc, f_with_fixed_kc, r.signal, qy.signal, qu.signal);

[u_cl, y_cl, err_cl] = simulate_sys_lsim(sim.time, k*g, cl, f, r.signal, qy.signal, qu.signal);
[u_cl_with_fixed_kc, y_cl_with_fixed_kc, err_cl_with_fixed_kc] = simulate_sys_lsim(sim.time,
    k*g, cl_with_fixed_kc, f_with_fixed_kc, r.signal, qy.signal, qu.signal);

disp('======================== LeadLag with Kc equal 1 =======================')
err_inf = 1 / (1 + (c(0)*g(0)))
y_inf = 1 - err_inf

err_qu_inf = err_inf * g(0)
y_qu_inf = 1 - err_qu_inf

disp('====================== LeadLag with Kc equal 1,43 =====================')
err_with_fixed_kc_inf = 1 / (1 + c_with_fixed_kc(0)*g(0))
y_with_fixed_kc_inf = 1 - err_with_fixed_kc_inf

err_qu_with_fixed_kc_inf = err_with_fixed_kc_inf*g(0)
y_qu_with_fixed_kc_inf = 1 - err_qu_with_fixed_kc_inf
disp('=======================================================================')

% =============================================================================
% Plot Graphs
% =============================================================================

c_str = '$C(jw)$ p/ $K_{c} = 1$';
c_with_fixed_kc_str = '$C(jw)$ p/ $K_{c} = 1$';

clead_str =  "$C\'(jw)\\overline{K}G(jw)$ p/ $K_{c} = 1$";
clead_with_fixed_kc_str =  "$C\'(jw)\\overline{K}G(jw)$ p/ $K_{c} = 1,43$";

clag_str = "$C_{at}(jw)$ p/ $K_{c} = 1$";
clag_with_fixed_kc_str = "$C_{at}(jw)$ p/ $K_{c} = 1,43$";

cleadlag_str = "$C\'(jw)\\overline{K}C_{at}$ p/ $K_c = 1$";
cleadlag_with_fixed_kc_str = "$C\'(jw)\\overline{K}C_{at}$ p/ $K_c = 1,43$";

openloop_lead_str = "$C\'(jw)\\overline{K}G(jw)$ p/ $K_c = 1$";
openloop_lead_with_fixed_kc_str = "$C\'(jw)\\overline{K}G(jw)$ p/ $K_c = 1,43$";

openloop_leadlag_str = "$C(jw)G(jw)$ p/ $K_c = 1$";
openloop_leadlag_with_fixed_kc_str = "$C(jw)G(jw)$ p/ $K_c = 1,43$";


figure(1)
plot_bode (clag, w, 'r', '--', clag_str, [2 1 1], [2 1 2])
plot_bode (clag_with_fixed_kc, w, 'b', '--', clag_with_fixed_kc_str, [2 1 1], [2 1 2])

plot_bode (c, w, 'r', '-', cleadlag_str, [2 1 1], [2 1 2])
plot_bode (c_with_fixed_kc, w, 'b', '-', cleadlag_with_fixed_kc_str, [2 1 1], [2 1 2])

subplot(2,1,1)
set(legend, 'location', 'northeastoutside')
ylim([-5 30])
subplot(2,1,2)
set(legend, 'location', 'northeastoutside')

figure(2)

plot_bode (c*g, w, 'r', '-', openloop_leadlag_str, [2 1 1], [2 1 2])
plot_phase_margin(c*g, [2 1 1], [2 1 2], 'r')

plot_bode (c_with_fixed_kc*g, w, 'b', '-', openloop_leadlag_with_fixed_kc_str, [2 1 1], [2 1 2])
plot_phase_margin(c_with_fixed_kc*g, [2 1 1], [2 1 2], 'b')

plot_bode (cl*k*g, w, 'r', '--', openloop_lead_str, [2 1 1], [2 1 2])
plot_phase_margin(cl*k*g, [2 1 1], [2 1 2], 'b')

plot_phase_margin(cl_with_fixed_kc*k*g, [2 1 1], [2 1 2], 'r')
plot_bode (cl_with_fixed_kc*k*g, w, 'b', '--', openloop_lead_with_fixed_kc_str, [2 1 1], [2 1 2])
        
subplot(2,1,1)
set(legend, 'location', 'northeastoutside')
subplot(2,1,2)
set(legend, 'location', 'northeastoutside')

figure(3)

subplot(2,1,1)
plot_signal(sim.time, '', r.signal, '', '--k', '$r(t)$')
hold on

unit_kc_str = ' p/ $K_{c} = 1$';
fixed_kc_str = [' p/ $K_{c} = ' num2str(kc, '%.2f') '$'];

plot_signal (sim.time, '',
             y, '',
             line_color = 'b-',
             ['$y(t)$ c/ $C_{at}$ e' unit_kc_str])
plot_signal (sim.time, 'Tempo(jw)',
             y_with_fixed_kc, 'Saída $y(t)$',
             line_color = 'r-',
             ['$y(t)$ c/ $C_{at}$ e' fixed_kc_str])
                                    
hold on
plot_signal (sim.time, '',
             y_cl, '',
             line_color = 'c-',
             ['$y(t)$  s/ $C_{at}$ e' unit_kc_str])
plot_signal (sim.time, 'Tempo(jw)',
             y_cl_with_fixed_kc, 'Saída $y(t)$',
             line_color = 'm-',
             ['$y(t)$ s/ $C_{at}$ e' fixed_kc_str])

subplot(2,1,2)
plot_signal (sim.time, '',
             err, '',
             line_color = 'b-',
             ['$e(t)$ c/ $C_{at}$ e' unit_kc_str])
hold on                                    
plot_signal (sim.time, 'Tempo(jw)',
             err_with_fixed_kc, 'Erro $e(t)$',
             line_color = 'r-',
             ['$e(t)$ c/ $C_{at}$ e' fixed_kc_str])
plot_signal (sim.time, '',
             err_cl, '',
             line_color = 'c-',
             ['$e(t)$ s/ $C_{at}$ e' unit_kc_str])
plot_signal (sim.time, 'Tempo(jw)',
             err_cl_with_fixed_kc, 'Erro $e(t)$',
             line_color = 'm-',
             ['$e(t)$ s/ $C_{at}$ e' fixed_kc_str])

subplot(2,1,1)
set(legend, 'location', 'northeastoutside')
subplot(2,1,2)
set(legend, 'location', 'northeastoutside')
