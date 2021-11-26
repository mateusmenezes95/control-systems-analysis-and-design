printf("Loading graphs plot functions...\n");

% =============================================================================
% Defaul parameters
% =============================================================================

global font_size = 10;
global line_thickness = 1;
global y_axis_limits_offset = 0.2;

% =============================================================================
% Functions
% =============================================================================

function [y_min, y_max] = get_y_axis_limits (signal)
    global y_axis_limits_offset;
    y_min = min(signal) - y_axis_limits_offset;
    y_max = max(signal) + y_axis_limits_offset;
end

function plot_responses_of_disturbances_signals (sim_time, figure_num = 1,
                                                 reference,
                                                 control_loop_response,
                                                 controller_output,
                                                 input_disturbance,
                                                 output_disturbance)
    global font_size;
    global line_thickness;
    
    figure(figure_num)

    subplot(2, 2, 1)
    plot(sim_time, reference, '--r', 'linewidth', (line_thickness-0.5))
    hold on
    plot(sim_time, control_loop_response, 'b', 'linewidth', line_thickness)
    grid on
    [y_min, y_max] = get_y_axis_limits(control_loop_response);
    axis([0 sim_time(end) y_min y_max])
    set(gca,'fontsize',font_size)

    hx = xlabel('Tempo (s)');
    hy = ylabel('Saída $y(t)$');
    hl = legend('$r(t)$', '$y(t)$');
    set(hl, 'fontsize', font_size)
    set(hl, 'location', 'southeast')
    set(hx, 'fontsize', font_size)
    set(hy, 'fontsize', font_size)

    subplot(2, 2, 3)
    plot(sim_time, controller_output, 'b', 'linewidth', line_thickness)
    grid on
    [y_min, y_max] = get_y_axis_limits(controller_output);
    axis([0 sim_time(end) y_min y_max])
    set(gca,'fontsize',font_size)

    hx = xlabel('Tempo (s)');
    hy = ylabel('Sinal de Controle $u(t)$');
    set(hx,'fontsize', font_size)
    set(hy,'fontsize', font_size)

    subplot(2, 2, 2)
    plot(sim_time, output_disturbance, 'b', 'linewidth', line_thickness)
    grid on
    [y_min, y_max] = get_y_axis_limits(output_disturbance);
    axis([0 sim_time(end) y_min  y_max])
    set(gca,'fontsize', font_size)

    hx = xlabel('Tempo (s)');
    hy = ylabel('Perturbação na Saída $q_{y}(t)$');
    set(hx, 'fontsize', font_size)
    set(hy, 'fontsize', font_size)

    subplot(2, 2, 4)
    plot(sim_time, input_disturbance,'b','linewidth',line_thickness)
    grid on
    [y_min, y_max] = get_y_axis_limits(input_disturbance);
    axis([0 sim_time(end) y_min y_max])
    set(gca, 'fontsize', font_size)

    hx = xlabel('Tempo (s)');
    hy = ylabel('Perturbação na Entrada $q_{u}(t)$');
    set(hx, 'fontsize', font_size)
    set(hy, 'fontsize', font_size)
end

function plot_signal (x, x_name,
                      y, y_name,
                      line_color,
                      signal_legend)
    global font_size;
    global line_thickness;
    
    plot(x, y, [line_color ';' signal_legend ';'], 'linewidth', (line_thickness - 0.5))
    grid on
    [y_min, y_max] = get_y_axis_limits(y);

    axis([0 x(end) min(ylim()(1), y_min) max(ylim()(2), y_max)]);

    set(gca, 'fontsize', font_size)

    xlabel(x_name, 'fontsize', font_size);
    ylabel(y_name, 'fontsize', font_size);
end

function plot_response_and_control_signals (sim_time, figure_num = 1,
                                            reference = 'none',
                                            y, y_legend,
                                            err, err_legend, 
                                            u, u_legend,
                                            line_color)
    global font_size;
    global line_thickness;
    
    figure(figure_num)

    subplot(3,1,1)
    if ((!ischar(reference)) & (length(get(gca, 'children')) == 0))
        plot_signal(sim_time, '', reference, '', '--r', '$r(t)$')
    endif

    hold on

    plot_signal(sim_time, 'Tempo (s)',
                y, 'Saida $y(t)$',
                line_color, y_legend)

    subplot(3,1,2)
    plot_signal(sim_time, 'Tempo (s)',
                err, 'Erro $e(t)$',
                line_color, err_legend)

    hold on

    subplot(3,1,3)
    plot_signal(sim_time, 'Tempo (s)',
                u, 'Sinal de Controle $u(t)$',
                line_color, u_legend)
end

function generate_subplot (which_subplot = [1 1 1])
    subplot(which_subplot(1), which_subplot(2), which_subplot(3))
endfunction

function plot_coordinates(x, y, which_subplot = [1 1 1], point_description = '',
                          color = 'k', linestyle = ':', point_symbol = '.')
    global line_thickness;
    generate_subplot(which_subplot)
    set(legend, 'autoupdate', 'off')
    line([x x], [ylim()(1) y],
        'color', color, 'linewidth', (line_thickness - 0.5), 'linestyle', linestyle)
    line([xlim()(1) x], [y y],
        'color', color, 'linewidth', (line_thickness - 0.5), 'linestyle', linestyle)
    hold on
    plot(x, y, point_symbol, 'color', color)

    if !isnull(point_description)
        text(x-1.5, ylim()(1) + 0.08, point_description)
    end
end

function plot_semilogx(x, y, color = 'b', linestyle = '-', legend_str = '',
                       xtitle = '', ytitle = '')
    global line_thickness;
    semilogx(x, y, [color ';' legend_str ';'], 'linestyle', linestyle,
             'linewidth', (line_thickness - 0.5));
    xlabel(xtitle);
    ylabel(ytitle);
    grid on
endfunction

function plot_bode(sys, w, color, linestyle, legend_str,
                   gain_subplot = [2 1 1], phase_subplot = [2 1 2])
    [gain, phase, w] = bode(sys, w);

    generate_subplot(gain_subplot)
    hold on
    plot_semilogx(w, 20*log10(gain), color, linestyle, legend_str,
                  'Frequência [rad/s]', 'Ganho [db]')
    generate_subplot(phase_subplot)
    hold on
    plot_semilogx(w, phase, color, linestyle, legend_str,
                  'Frequência [rad/s]', 'Fase [graus]')
endfunction

function plot_phase_margin(sys, gain_subplot = [2 1 1], phase_subplot = [2 2 2],
                           color = 'k');
    global line_thickness;
    [gain_margin, phase_margin, w_gain_margin, w_phase_margin] = margin(sys);

    generate_subplot(gain_subplot)
    hold on
    set(legend, 'autoupdate', 'off')

    line(xlim(), [0 0], 'color', 'k', 'linewidth', line_thickness/2, 'linestyle', ':')
    line([w_phase_margin w_phase_margin], [ylim()(1) 0],
         'color', 'r', 'linewidth', line_thickness/2, 'linestyle', '-')
    set(legend, 'location', 'southwest')
    plot(w_phase_margin, 0, 'marker', '.', 'color', 'r')
    
    generate_subplot(phase_subplot)
    set(legend, 'autoupdate', 'off')

    line(xlim(), [-180 -180], 'color', 'k', 'linewidth', 0.5, 'linestyle', ':')

    y = (-180 + phase_margin);

    line([w_phase_margin w_phase_margin], [-180 y],
         'color', 'r', 'linewidth', line_thickness/2, 'linestyle', '-')
    % line([xlim()(1) w_phase_margin], [y y],
    % 'color', color, 'linewidth', (line_thickness - 0.5), 'linestyle', ':')
    set(legend, 'location', 'southwest')
    hold on
    plot(w_phase_margin, y, 'marker', '.', 'color', 'r')
endfunction

function plot_phase_peak(sys, w, gain_subplot = [2 1 1], phase_subplot = [2 2 2],
                         color = 'k');
    global line_thickness;
    [gain, phase, w] = bode(sys, w);

    [phi_max, phi_max_idx] = max(phase);

    plot_coordinates(w(phi_max_idx), 20*log10(gain(phi_max_idx)), gain_subplot,
                     '', 'k', '--', '.')
    plot_coordinates(w(phi_max_idx), phi_max, phase_subplot,
                     '', 'k', '--', '.')
endfunction

printf("Loaded successfuly graphs plot functions \n");
