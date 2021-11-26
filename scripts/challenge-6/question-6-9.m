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
xt_with_obsv = xt;

xd_est(:,1) = xk0;

ut(:,1) = [0;0];
ut_with_dist(:,1) = ut(:,1) + qu(:,1);
ut_with_obsv = ut;

k=1;

if simulate_regulator
    for i=1:sim_time_length
        xt(:,i+1) = a_euler*xt(:,i) + b_euler*ut(:,i);
        xt_with_dist(:,i+1) = a_euler*xt_with_dist(:,i) + b_euler*ut_with_dist(:,i);
        xt_with_obsv(:,i+1) = a_euler*xt_with_obsv(:,i) + b_euler*ut_with_obsv(:,i);

        yt(:,i) = gs_ss.c*xt(:,i);
        yt_with_dist(:,i) = gs_ss.c*xt_with_dist(:,i) + qy(:,i);
        yt_with_obsv(:,i) = gs_ss.c*xt_with_obsv(:,i);
 
        if (mod(i, integration_step_ratio) == 1) || i == 1
            t1 = (i-1)*integration_step_size;
            t2 = (k-1)*sampling_period;
            
            xs(:,k) = xt(:,i);
            ys(:,k) = yt(:,i);

            us(:,k) = -kd*(xs(:,k) - xeq_discrete.*unit_step(:,i)) + ueq(:,i);
            us_with_obsv(:,k) = -kd*(xd_est(:,k) - xeq_discrete.*unit_step(:,i)) + ueq(:,i);

            if i != sim_time_length
                xd_est(:,k+1) = gz_a_minus_lc*xd_est(:,k) + gz_ss.b*us(:,k) + ld*ys(:,k);
            endif

            ed_est(:,k) = xs(:,k) - xd_est(:,k);

            k = k+1;
        endif

        if i != sim_time_length
            ut(:,i+1) = us(:,k-1);
            ut_with_obsv(:,i+1) = us_with_obsv(:,k-1);
            ut_with_dist(:,i+1) = ut(:,i) + qu(:,i);
            yheld(:, i) = ys(:,k-1);
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

    plot_mimo_response_and_control_signals(sim.time, figure_num = 3,
                                           yt_with_obsv, ut_with_obsv, ref, 'b');

    plot_mimo_response_and_control_signals(sim.time, figure_num = 3,
                                           yt, ut, ref, 'r', plot_ref = false);

    set_mimo_signals_graph_legend(fig = 3, x_pos = 0.25,
                                  {'referencia', 'sem observador', 'com observador'})
                                  
    plot_mimo_states_estimation_error(size(ed_est)(2), ed_est, 'b', 4)

    plot_mimo_output_estimation_error(sim.time, y_est_err, 'b', 5)
endif

zero_bar = zeros(4,2);
idty = eye(2,2);
a_bar = [gz_ss.a zero_bar; gz_ss.c*gz_ss.a idty];
b_bar = [gz_ss.b; gz_ss.c*gz_ss.b];
desirable_ctrb_eigenvalues = [desirable_ctrb_eigenvalues desirable_ctrb_eigenvalues(1:2)];

k_bar = place(a_bar , b_bar , desirable_ctrb_eigenvalues);
kx = k_bar(:,1:4);
ki = k_bar(:,5:6);

% Reference following with disturbance rejection

us(:,1:3) = zeros(2,3);
xs(:,1:3) = [xk0 xk0 xk0];
ys(:,1:3) = zeros(2,3);
x(:,1) = xk0;
y(:,1) = [0;0];
k = 2;
if simulate_reference_follower
    for i=1:length(sim.time)
        if (mod(i, integration_step_ratio) == 1);
            k=k+1;
            xs(:,k) = x(:,i);
            ys(:,k) = y(:,i);
            yrs(:,k) = yr(:,i);
            us(:,k) = us(:,k-1) - kx*(xs(:,k) - xs(:,k-1)) + ki*(yrs(:,k) - ys(:,k));
        endif

        u(:,i) = us(:,k);
        [x(:,i+1), y(:,i+1)] = get_ss_output(x(:,i), gs_ss_for_euler, us(:,k) + qu(:,i));
        y(:,i+1) = y(:,i+1) + qy(:,i);
    endfor

    plot_mimo_response_and_control_signals(sim.time, figure_num = 2,
                                           y(:,2:end), u, ref,
                                           'r', plot_ref = false);

    set_mimo_signals_graph_legend(fig = 2, x_pos = 0.25,
                                  {'referencia', 'sem acao integral', 'com acao integral'})
endif
