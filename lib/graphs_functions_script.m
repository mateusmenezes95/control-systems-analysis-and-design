printf("Loading graphs plot functions...\n");

% ===============================================================================
% Defaul parameters
% ===============================================================================

global font_size = 15
global line_thickness = 2;

% ===============================================================================
% Functions
% ===============================================================================

function plot_responses_of_disturbances_signals (sim_time, reference, control_loop_response,
                                                 controller_output, input_disturbance, output_disturbance)
    global font_size;
    global line_thickness;
    
    figure(1)

    subplot(2, 2, 1)
    plot(sim_time, reference, '--r', 'linewidth', (line_thickness-0.5))
    hold on
    plot(sim_time, control_loop_response, 'b', 'linewidth', line_thickness)
    grid on
    axis([0 sim_time(end) (min(control_loop_response) - 0.2)  (max(control_loop_response) + 0.2)])
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
    axis([0 sim_time(end) (min(controller_output) - 0.2) (max(controller_output) + 0.2)])
    set(gca,'fontsize',font_size)

    hx = xlabel('Tempo (s)');
    hy = ylabel('Sinal de Controle');
    set(hx,'fontsize', font_size)
    set(hy,'fontsize', font_size)

    subplot(2, 2, 2)
    plot(sim_time, output_disturbance, 'b', 'linewidth', line_thickness)
    grid on
    axis([0 sim_time(end) min(output_disturbance)-0.2  max(output_disturbance)+0.2])
    set(gca,'fontsize', font_size)

    hx = xlabel('Tempo (s)');
    hy = ylabel('Perturbação - Saída');
    set(hx, 'fontsize', font_size)
    set(hy, 'fontsize', font_size)

    subplot(2, 2, 4)
    plot(sim_time, input_disturbance,'b','linewidth',line_thickness)
    grid on
    axis([0 sim_time(end) min(input_disturbance)-0.2 max(input_disturbance)+0.2])
    set(gca, 'fontsize', font_size)

    hx = xlabel('Tempo (s)');
    hy = ylabel('Perturbação - Entrada');
    set(hx, 'fontsize', font_size)
    set(hy, 'fontsize', font_size)
end
