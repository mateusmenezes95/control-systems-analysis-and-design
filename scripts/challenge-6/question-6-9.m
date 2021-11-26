run base_script.m
run specific_state_space_graphs_funtions.m

simulate_regulator = true;
simulate_reference_follower = true;

for i=1:4
    desirable_eigenvalues(i) = e^(-sampling_period/4);
endfor

kd = place(gz_ss.A, gz_ss.B, desirable_eigenvalues)

% Initial contidions
xk0 = zeros(size(xeq_continue));

a_euler = eye(size(gs_ss.A))+(integration_step_size*gs_ss.A);
b_euler = integration_step_size*gs_ss.B;

gs_ss_for_euler = ss(a_euler, b_euler, gs_ss.c, gs_ss.d)

% Regulator action
ueq = ueq.*unit_step(1:2,:);
yr = yeq.*unit_step(1:2,:);
ref = [yr; ueq];
if simulate_regulator
    for k=1:length(sim.time)
        if (mod(k, integration_step_ratio) == 1)
            t = (k-1)*integration_step_size;
            if k == 1
                uk = -kd*(xk0 - xeq_discrete.*unit_step(:,k)) + ueq(:,k);
            else
                uk = -kd*(xk - xeq_discrete.*unit_step(:,k)) + ueq(:,k);
            endif
        endif
        ukd = uk + qu(:,k);

        if k == 1
            [xk, yk] = get_ss_output(xk0, gs_ss_for_euler, uk);
            [xkd, ykd] = get_ss_output(xk0, gs_ss_for_euler, ukd);
            Y = yk;
            Yd = ykd + qy(:,k);
            U = uk;
            Ud = ukd;
            X = xk;
        else
            [xk, yk] = get_ss_output(xk, gs_ss_for_euler, uk);
            [xkd, ykd] = get_ss_output(xkd, gs_ss_for_euler, ukd);
            Y = [Y yk];
            ykd = ykd + qy(:,k); 
            Yd = [Yd ykd];
            U = [U uk];
            Ud = [Ud ukd];
            X = [X xk];
        endif
    endfor

    plot_mimo_response_and_control_signals(sim.time, figure_num = 1,
                                           Y, U, ref, 'b', {'sem integrador'});

    fig_handle = figure(1);
    set(fig_handle, 'units', 'points')
    fig_pos_vec = get(fig_handle, 'position')
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
                                           Yd, Ud, ref, 'b', {'sem integrador'});
endif

zero_bar = zeros(4,2);
idty = eye(2,2);
a_bar = [gz_ss.a zero_bar; gz_ss.c*gz_ss.a idty];
b_bar = [gz_ss.b; gz_ss.c*gz_ss.b];
desirable_eigenvalues = [desirable_eigenvalues desirable_eigenvalues(1:2)];

k_bar = place(a_bar , b_bar , desirable_eigenvalues);
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
                                           'r', {'sem integrador','com integrador'}, false);
    
    fig_handle = figure(2);
    set(fig_handle, 'units', 'points')
    fig_pos_vec = get(fig_handle, 'position')
    fig_pos_vec(3) = 400;
    fig_pos_vec(4) = 250;
    set(fig_handle, 'position', fig_pos_vec)

    subplot(2,2,1)
    legend('show')
    set(legend, 'location', 'northoutside')
    set(legend, 'orientation', 'horizontal')
    legend({'referencia', 'sem acao integral', 'com acao integral'})
    legend_pos = get(legend, 'position');
    legend_pos(1) = 0.25;
    set(legend, 'position', legend_pos);
endif
