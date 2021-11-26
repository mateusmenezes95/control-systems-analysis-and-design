run base_script.m
run specific_state_space_graphs_funtions.m

simulate_regulator = true;
simulate_reference_follower = true;

for i=1:4
    desirable_ctrb_eigenvalues(i) = e^(-sampling_period/4);
    desirable_obsv_eigenvalues(i) = e^(-sampling_period/2);
endfor

kd = place(gz_ss.A, gz_ss.B, desirable_ctrb_eigenvalues);
ld = place(gz_ss.A', gz_ss.C', desirable_obsv_eigenvalues)';

zero_bar = zeros(4,2);
idty = eye(2,2);
a_bar = [gz_ss.a zero_bar; gz_ss.c*gz_ss.a idty];
b_bar = [gz_ss.b; gz_ss.c*gz_ss.b];
desirable_ctrb_eigenvalues = [desirable_ctrb_eigenvalues desirable_ctrb_eigenvalues(1:2)];

k_bar = place(a_bar , b_bar , desirable_ctrb_eigenvalues);
kx = k_bar(:,1:4);
ki = k_bar(:,5:6);

% Initial contidions
xk0 = zeros(size(xeq_continue));
xk0_est = xk0 .+ 0.01;

a_euler = eye(size(gs_ss.A))+(integration_step_size*gs_ss.A);
b_euler = integration_step_size*gs_ss.B;

gz_a_minus_lc = gz_ss.a - ld*gz_ss.c;

gs_ss_for_euler = ss(a_euler, b_euler, gs_ss.c, gs_ss.d);

% Regulator action
ueq = ueq.*unit_step(1:2,:);
yr = yeq.*unit_step(1:2,:);
ref = [yr; ueq];

% Initial conditions
xt(:,1) = xk0;
xt_with_dist(:,1) = xk0;
xt_with_obsv(:,1) = xk0;
xt_with_int(:,1) = xk0;

xd(:,1) = xk0;
xd_with_integrator(:,1) = xk0;
xd_est(:,1:2) = zeros(4,2);

u0 = [0; 0];
ut(:,1) = u0;
ut_with_dist(:,1) = u0 + qu(:,1);
ut_with_obsv(:,1) = u0;
ut_with_int(:,1) = u0;

ud_with_int(:,1) = [0; 0];

k=2;

if simulate_regulator
    for i=1:sim_time_length
        xt(:,i+1) = a_euler*xt(:,i) + b_euler*ut(:,i);
        xt_with_dist(:,i+1) = a_euler*xt_with_dist(:,i) + b_euler*ut_with_dist(:,i);
        xt_with_int(:,i+1) = a_euler*xt_with_int(:,i) + b_euler*ut_with_int(:,i);
        xt_with_obsv(:,i+1) = a_euler*xt_with_obsv(:,i) + b_euler*ut_with_obsv(:,i);

        yt(:,i) = gs_ss.c*xt(:,i);
        yt_with_dist(:,i) = gs_ss.c*xt_with_dist(:,i) + qy(:,i);
        yt_with_int(:,i) = gs_ss.c*xt_with_int(:,i) + qy(:,i);
        yt_with_obsv(:,i) = gs_ss.c*xt_with_obsv(:,i);
 
        if (mod(i, integration_step_ratio) == 1) || i == 1
            t1 = (i-1)*integration_step_size;
            t2 = (k-1)*sampling_period;
            
            xd(:,k) = xt(:,i);
            xd_with_integrator(:,k) = xt_with_int(:,i);

            yd(:,k) = yt(:,i);
            yd_with_int(:,k) = yt_with_int(:,i);
            yrd(:,k) = yr(:,i);

            ud(:,k) = -kd*(xd(:,k) - xeq_discrete.*unit_step(:,i)) + ueq(:,i);
            ud_with_obsv(:,k) = -kd*(xd_est(:,k) - xeq_discrete.*unit_step(:,i)) + ueq(:,i);

            ud_with_int(:,k) = (ud_with_int(:,k-1)
                             - kx*(xd_with_integrator(:,k) - xd_with_integrator(:,k-1))
                             + ki*(yrd(:,k) - yd_with_int(:,k)));

            if i != sim_time_length
                xd_est(:,k+1) = gz_a_minus_lc*xd_est(:,k) + gz_ss.b*ud(:,k) + ld*yd(:,k);
            endif

            ed_est(:,k) = xd(:,k) - xd_est(:,k);

            k = k+1;
        endif

        if i != sim_time_length
            ut(:,i+1) = ud(:,k-1);
            ut_with_obsv(:,i+1) = ud_with_obsv(:,k-1);
            ut_with_dist(:,i+1) = ut(:,i) + qu(:,i);
            ut_with_int(:,i+1) = ud_with_int(:,k-1) + qu(:,i);
            yheld(:, i) = yd(:,k-1);
        endif
    endfor

    y_est_err = yt - yt_with_obsv;

    plot_mimo_response_and_control_signals(sim.time, figure_num = 1,
                                           yt, ut, ref,
                                           'b');

    fig_handle = figure(1);
    set(fig_handle, 'units', 'points')
    fig_pos_vec = get(fig_handle, 'position');
    fig_pos_vec(3) = 400;
    fig_pos_vec(4) = 250;
    set(fig_handle, 'position', fig_pos_vec)

    subplot(2,2,2)
    legend('show')
    set(legend, 'location', 'northoutside')
    set(legend, 'orientation', 'horizontal')
    legend({'equilibrio'})
    legend_pos = get(legend, 'position');
    legend_pos(4) = 0.044055;
    set(legend, 'position', legend_pos);

    plot_mimo_response_and_control_signals(sim.time, figure_num = 2,
                                           yt_with_dist, ut_with_dist, ref, 'b');

    plot_mimo_response_and_control_signals(sim.time, figure_num = 2,
                                           yt_with_int, ut_with_int, ref,
                                           'r', plot_ref = false);

    set_mimo_signals_graph_legend(fig = 2, x_pos = 0.25,
                                  {'referencia', 'sem acao integral', 'com acao integral'})

    plot_mimo_response_and_control_signals(sim.time, figure_num = 3,
                                           yt_with_obsv, ut_with_obsv, ref, 'b');

    plot_mimo_response_and_control_signals(sim.time, figure_num = 3,
                                           yt, ut, ref, 'r', plot_ref = false);

    set_mimo_signals_graph_legend(fig = 3, x_pos = 0.25,
                                  {'referencia', 'sem observador', 'com observador'})
                                  
    plot_mimo_states_estimation_error(size(ed_est)(2), ed_est, 'b', 4)

    plot_mimo_output_estimation_error(sim.time, y_est_err, 'b', 5)
endif
