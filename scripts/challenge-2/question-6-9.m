addpath(genpath("./")) % Add lib path to Octave script file search paths
run common_parameters_scripts.m

% =============================================================================
% Transfer Functions Definions
% =============================================================================

simulate_question(6:9) = [true true true true]
plot_saturation_analysis = true

disp('=======================================================================')
disp('Simulating plant response with input reference filter')
disp('=======================================================================')

model_num = 1
delay = l(model_num)
c = kc * ((s*ti + 1) / s*ti)
ci = kc / (s*ti)
g = g = k(model_num) / (tau(model_num) * s + 1)
f = 1 / ((ti*s) + 1)
f_unit = tf(1)

% =============================================================================
% Main of the script
% =============================================================================

if simulate_question(6)
    [u, y, err] = simulate_sys(sim.time, dt,
                               g, delay, c, f,
                               reference.signal,
                               step_output_disturbance.signal,
                               step_input_disturbance.signal);

    reference_filtered = step(sim.time, f);
    time_delay_size = length(0:dt:2);
    reference_filtered = reference_filtered(1:(length(reference_filtered)-time_delay_size), 1);
    reference_filtered = [zeros(time_delay_size, 1)', reference_filtered(1:end)'];
    subplot(3,1,1)
    plot(sim.time, reference_filtered, '-r;$r_{f}(t)$;');
    hold on

    plot_response_and_control_signals (sim.time, figure_num = 1,
                                      reference.signal,
                                      y, '$y_{m1}(t)$ com filtro de referência',
                                      err, '$e_{m1}(t)$ com filtro de referência',
                                      u, '$u_{m1}(t)$ com filtro de referência', 'b-')
end

if simulate_question(7)
    disp('=======================================================================')
    disp('Simulating plant response with I+P controller')
    disp('=======================================================================')

    [u, y, err] = simulate_sys(sim.time, dt,
                            g, delay, ci, f_unit,
                            reference.signal,
                            step_output_disturbance.signal,
                            step_input_disturbance.signal,
                            controller_type = 'I+P');

    hold on

    plot_response_and_control_signals (sim.time, figure_num = 1, 'none',
                                    y, '$y_{m1}(t)$ com controlador I+P',
                                    err, '$e_{m1}(t)$ com controlador I+P',
                                    u, '$u_{m1}(t)$ com controlador I+P', 'm-')
end

if simulate_question(8)
    disp('============================================================================')
    disp('Simulating plant response with I+P controller and control outuput saturation')
    disp('============================================================================')

    [u, y, err] = simulate_sys(sim.time, dt,
                            g, delay, ci, f_unit,
                            reference.signal,
                            step_output_disturbance.signal,
                            step_input_disturbance.signal,
                            controller_type = 'I+P',
                            saturation = [0 0.8]);

    hold on

    plot_response_and_control_signals (sim.time, figure_num = 1, 'none',
                                    y, '$y_{m1}(t)$ com saturação de $u(t)$',
                                    err, '$e_{m1}(t)$ com saturação',
                                    u, '$u_{m1}(t)$ com saturação', 'g-')
    if plot_saturation_analysis
        sat_idxs(1) = find(u >= 0.8)(1);
        sat_idxs(2)= find(u(1, ++sat_idxs(1):end) < 0.8)(1) + sat_idxs(1);
        sat_times = [(--sat_idxs(1) * dt) (--sat_idxs(2) * dt)]; % Time starts from zero!!!
        err_cross_idx(1) = find(err < 0)(1);
        err_cross_idx(2) = find(err(1, ++err_cross_idx(1):end) >= 0)(1) + err_cross_idx(1);
        neg_e_time = [(--err_cross_idx(1) * dt) (--err_cross_idx(2) * dt)];

        plot_coordinates(sat_times(1), y(1, sat_idxs(1)), [3 1 1], '$t_{1}$')
        plot_coordinates(neg_e_time(1), y(1, err_cross_idx(1)), [3 1 1], '$t_{2}$')
        plot_coordinates(sat_times(2), y(1, sat_idxs(2)), [3 1 1], '$t_{3}$')
        plot_coordinates(neg_e_time(2), y(1, err_cross_idx(1)), [3 1 1], '$t_{4}$')
        set(gca, 'ylim', [-0.2 1.2]); % This is called in Brazil "armengue!"
        set(gca, 'ytick', [-0.2:.2:1.2]); % This is called in Brazil "armengue!"
        set(gca, 'box', 'on'); % This is called in Brazil "armengue!"
        
        plot_coordinates(sat_times(1), err(1, sat_idxs(1)), [3 1 2], '$t_{1}$')
        plot_coordinates(neg_e_time(1), err(1, err_cross_idx(1)), [3 1 2], '$t_{2}$')
        plot_coordinates(sat_times(2), err(1, sat_idxs(2)), [3 1 2], '$t_{3}$')
        plot_coordinates(neg_e_time(2), err(1, err_cross_idx(2)), [3 1 2], '$t_{4}$')
        
        plot_coordinates(sat_times(1), u(1, sat_idxs(1)), [3 1 3], '$t_{1}$')
        plot_coordinates(neg_e_time(1), u(1, err_cross_idx(1)), [3 1 3], '$t_{2}$')
        plot_coordinates(sat_times(2), u(1, sat_idxs(2)), [3 1 3], '$t_{3}$')
        plot_coordinates(neg_e_time(2), u(1, err_cross_idx(2)), [3 1 3], '$t_{4}$')
    end
end

if simulate_question(9)
    disp('====================================================================')
    disp('Simulating plant response with I+P , saturation and anti-windup')
    disp('====================================================================')

    [u, y, err] = simulate_sys(sim.time, dt,
                            g, delay, ci, f_unit,
                            reference.signal,
                            step_output_disturbance.signal,
                            step_input_disturbance.signal,
                            controller_type = 'I+P',
                            saturation = [0 0.8],
                            anti_windup_tol = 10e-5);

    hold on

    plot_response_and_control_signals (sim.time, figure_num = 1, 'none',
                                    y, '$y_{m1}(t)$ com anti-windup',
                                    err, '$e_{m1}(t)$ com anti-windup',
                                    u, '$u_{m1}(t)$ com anti-windup', 'k-')

    set(gca, 'ytick', [0:.2:1.2]); % This is called in Brazil "armengue!"
end

if all(simulate_question(6:9))
    for i=1:3
        subplot(3,1,i);
        set(legend, 'location', 'eastoutside');
    end
end