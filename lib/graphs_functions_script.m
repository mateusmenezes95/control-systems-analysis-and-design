printf("Loading graphs plot functions...\n");

% =============================================================================
% Defaul parameters
% =============================================================================

global font_size = 15
global line_thickness = 2;
global y_axis_limits_offset = 0.2;

% =============================================================================
% Functions
% =============================================================================

function [y_min, y_max] = get_y_axis_limits (signal)
    global y_axis_limits_offset;
    y_min = min(signal) - y_axis_limits_offset;
    y_max = max(signal) + y_axis_limits_offset;
end

function plot_responses_of_disturbances_signals (sim_time, reference,
                                                 control_loop_response,
                                                 controller_output,
                                                 input_disturbance,
                                                 output_disturbance)
    global font_size;
    global line_thickness;
    
    figure(1)

    subplot(2, 2, 1)
    plot(sim_time, reference, '--r', 'linewidth', (line_thickness-0.5))
    hold on
    plot(sim_time, control_loop_response, 'b', 'linewidth', line_thickness)
    grid on
    [y_min, y_max] = get_y_axis_limits(control_loop_response);
    axis([0 sim_time(end) y_min y_max])
    set(gca,'fontsize',font_size)

    hx = xlabel('Tempo (s)');
    hy = ylabel('Saída');
    hl = legend('r(t)', 'y(t)');
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
    hy = ylabel('Sinal de Controle');
    set(hx,'fontsize', font_size)
    set(hy,'fontsize', font_size)

    subplot(2, 2, 2)
    plot(sim_time, output_disturbance, 'b', 'linewidth', line_thickness)
    grid on
    [y_min, y_max] = get_y_axis_limits(output_disturbance);
    axis([0 sim_time(end) y_min  y_max])
    set(gca,'fontsize', font_size)

    hx = xlabel('Tempo (s)');
    hy = ylabel('Perturbação - Saída');
    set(hx, 'fontsize', font_size)
    set(hy, 'fontsize', font_size)

    subplot(2, 2, 4)
    plot(sim_time, input_disturbance,'b','linewidth',line_thickness)
    grid on
    [y_min, y_max] = get_y_axis_limits(input_disturbance);
    axis([0 sim_time(end) y_min y_max])
    set(gca, 'fontsize', font_size)

    hx = xlabel('Tempo (s)');
    hy = ylabel('Perturbação - Entrada');
    set(hx, 'fontsize', font_size)
    set(hy, 'fontsize', font_size)
end
