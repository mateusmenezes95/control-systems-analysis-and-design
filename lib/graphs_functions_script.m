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

    hx = xlabel(x_name);
    hy = ylabel(y_name);
    legend('fontsize', font_size, 'location', 'southeast');
    set(hx, 'fontsize', font_size)
    set(hy, 'fontsize', font_size)
end

function plot_response_and_control_signals (sim_time, figure_num = 1,
                                            reference = 'none',
                                            y, y_legend, 
                                            u, u_legend,
                                            line_color)
    global font_size;
    global line_thickness;
    
    figure(figure_num)

    subplot(2, 1, 1)

    if !ischar(reference)
        plot_signal(sim_time, '', reference, '', '--r', 'r(t)')
    end

    hold on

    plot_signal(sim_time, 'Tempo (s)',
                y, 'Saida $y(t)$',
                line_color, y_legend)

    subplot(2,1,2)
    plot_signal(sim_time, 'Tempo (s)',
                u, 'Sinal de Controle $u(t)$',
                line_color, u_legend)
end
