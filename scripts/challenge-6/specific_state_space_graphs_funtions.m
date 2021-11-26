run graphs_functions_script

function plot_mimo_response_and_control_signals (sim_time, figure_num = 1,
                                                 y, u, reference, 
                                                 line_color, legend_str, plot_ref = true)
    global font_size;
    global line_thickness;
    
    figure(figure_num)

    subplot(2,2,1);
    if plot_ref
        plot_signal(sim_time, '', reference(1,:), '', '--k', '')
    endif
    hold on
    plot_signal(sim_time, 'Tempo (s)',
                y(1,:), 'Saida $y_1(t)$',
                line_color,'')
    legend('off')

    subplot(2,2,2);
    if plot_ref
        plot_signal(sim_time, '', reference(2,:), '', '--k', '')
    endif
    hold on
    plot_signal(sim_time, 'Tempo (s)',
                y(2,:), 'Saida $y_2(t)$',
                line_color, '')
    legend('off')

    subplot(2,2,3);
    if plot_ref
        plot_signal(sim_time, '', reference(3,:), '', '--k', '')
    endif
    hold on
    plot_signal(sim_time, 'Tempo (s)',
                u(1,:), 'Sinal de Controle $u_1(t)$',
                line_color, '')
    legend('off');

    subplot(2,2,4);
    if plot_ref
        plot_signal(sim_time, '', reference(4,:), '', '--k', '')
    endif
    hold on
    plot_signal(sim_time, 'Tempo (s)',
                u(2,:), 'Sinal de Controle $u_2(t)$',
                line_color, '')
    legend('off');
end
